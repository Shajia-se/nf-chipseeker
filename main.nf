#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


params.peaks         = params.peaks ?: "peaks/*"
params.outdir        = params.outdir ?: "chipseeker_output"
params.do_enrich     = (params.do_enrich ?: true) as boolean

def detectType(File f) {
    if (f.name.endsWith(".narrowPeak") || f.name.endsWith(".bed")) {
        return 'narrowPeak'
    }
    if (f.name.startsWith("all_peaks.") && f.name.endsWith(".tsv")) {
        return 'allpeaks'
    }
    if (f.name == "consensus_peaks.tsv") {
        return 'consensus'
    }
    return 'unknown'
}

workflow {

    Channel
        .fromPath(params.peaks)
        .map { f -> tuple(f, detectType(f)) }
        .filter { it[1] != 'unknown' }
        .tap { if (it.empty) error "No valid peak files found" }
        .set { peak_ch }


    narrowPeak_ch = peak_ch.filter { it[1] == 'narrowPeak' }
                           .map { tuple(it[0].baseName, it[0]) }

    allPeaks_ch   = peak_ch.filter { it[1] == 'allpeaks' }
                           .map { tuple(it[0].name.replaceFirst("all_peaks.","").replace(".tsv",""), it[0]) }

    consensus_ch  = peak_ch.filter { it[1] == 'consensus' }
                           .map { tuple("consensus", it[0]) }


    chipseeker_annotate(narrowPeak_ch)
    chipseeker_annotate(allPeaks_ch)
    chipseeker_annotate(consensus_ch)


    if (!consensus_ch.empty) {
        build_master_table(consensus_ch)
    }
}


process chipseeker_annotate {

    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(peak_file)

    output:
    path "*"

    script:
    """
    R --vanilla << 'EOF'
    suppressPackageStartupMessages({
        library(ChIPseeker)
        library(GenomicRanges)
        library(TxDb.Mmusculus.UCSC.mm39.knownGene)
        library(org.Mm.eg.db)
        library(clusterProfiler)
        library(ReactomePA)
        library(openxlsx)
        library(ggplot2)
    })

    txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
    OrgDb <- org.Mm.eg.db

    message("Reading peak file: ${peak_file}")
    fname <- "${peak_file}"

    # 自动识别 narrowPeak / all_peaks / consensus
    if (grepl(".narrowPeak$|.bed$", fname)) {
        df <- read.delim(fname, header=FALSE)
        colnames(df)[1:3] <- c("seqnames","start","end")
        df$strand <- "*"
        gr <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
        comp_name <- "${sample}"
    } else {
        df <- read.delim(fname, as.is=TRUE)
        comp_name <- "${sample}"
        gr <- GRanges(
          seqnames = df[, 'seqnames'],
          ranges   = IRanges(df[, 'start'], df[, 'end']),
          strand   = df[, 'strand'],
          peak_width  = df[, 'width'],
          Conc        = df[, 'Conc'],
          Fold        = df[, 'Fold'],
          FDR         = df[, 'FDR']
        )
        seqlevelsStyle(gr) <- "UCSC"
    }

    # ---------------------------
    # promoter heatmap + avg profile
    # ---------------------------
    promoters <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
    tagMatrix <- try(getTagMatrix(gr, windows=promoters, weightCol='Conc'), silent=TRUE)

    if (!inherits(tagMatrix, "try-error") && !is.null(dim(tagMatrix))) {
        pdf("tagMatrix_${sample}.pdf")
        tagHeatmap(tagMatrix, xlim=c(-2000,2000), color="red")
        plotAvgProf(tagMatrix, xlim=c(-2000,2000),
                    xlab="Genomic Region", ylab="Read Count Frequency")
        plotAvgProf(tagMatrix, xlim=c(-2000,2000), conf=0.95, resample=1000)
        dev.off()
    }

    # ---------------------------
    # ChIPseeker annotation
    # ---------------------------
    anno <- annotatePeak(gr, TxDb=txdb, annoDb=OrgDb,
                         tssRegion=c(-2000,2000))

    anno_df <- as.data.frame(anno)
    write.table(anno_df, "${sample}_annotated.tsv", sep="\\t",
                quote=FALSE, row.names=FALSE)
    write.xlsx(anno_df, "${sample}_annotated.xlsx", row.names=FALSE)

    # Pie + distToTSS
    pdf("annotation_${sample}.pdf", width=10, height=7)
        plotAnnoPie(anno)
        plotDistToTSS(anno)
    dev.off()

    # Barplot
    pdf("annotation_bar_${sample}.pdf")
        print(plotAnnoBar(anno))
    dev.off()

    # Violin plot (Fold by feature)
    if ("Fold" %in% colnames(anno_df)) {
        anno_df\$peak.group <- gsub("Intron.*","Intron",
                               gsub("Exon.*","Exon", anno_df\$annotation))
        pdf("feature_violin_${sample}.pdf", width=14, height=7)
        p <- ggplot(anno_df, aes(x=peak.group, y=Fold, fill=peak.group)) +
              geom_violin() +
              geom_hline(yintercept=0, lty=2) +
              ggtitle("${sample}") +
              stat_summary(geom='crossbar', size=0.2, fun=median)
        print(p)
        dev.off()
    }

    # ---------------------------
    # Enrichment (GO + KEGG + Reactome)
    # ---------------------------
    if (${params.do_enrich}) {
        genes <- unique(anno_df\$geneId)
        genes <- genes[!is.na(genes)]

        if (length(genes) >= 10) {

            # GO
            ego <- try(enrichGO(gene=genes, OrgDb=OrgDb, readable=TRUE), silent=TRUE)
            if (!inherits(ego, "try-error") && nrow(as.data.frame(ego))>0) {
                pdf("GO_dotplot_${sample}.pdf")
                print(dotplot(ego))
                dev.off()
            }

            # KEGG
            ekegg <- try(enrichKEGG(gene=genes, organism="mmu"), silent=TRUE)
            if (!inherits(ekegg, "try-error") && nrow(as.data.frame(ekegg))>0) {
                pdf("KEGG_dotplot_${sample}.pdf")
                print(dotplot(ekegg))
                dev.off()
            }

            # Reactome
            ereact <- try(enrichPathway(gene=genes, organism="mouse"), silent=TRUE)
            if (!inherits(ereact, "try-error") && nrow(as.data.frame(ereact))>0) {
                pdf("Reactome_dotplot_${sample}.pdf")
                print(dotplot(ereact))
                dev.off()
            }
        }
    }

    save.image("${sample}_session.Rdata")

    EOF
    """
}

process build_master_table {

    publishDir "${params.outdir}/master", mode:'copy', overwrite:true

    input:
    tuple val(s), path(consensus)

    script:
    """
    R --vanilla << 'EOF'
    suppressPackageStartupMessages({
        library(ChIPseeker)
        library(TxDb.Mmusculus.UCSC.mm39.knownGene)
        library(org.Mm.eg.db)
        library(openxlsx)
    })

    txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
    OrgDb <- org.Mm.eg.db

    D <- read.delim("${consensus}", as.is=TRUE)
    gr <- makeGRangesFromDataFrame(D, keep.extra.columns=TRUE,
              seqnames.field='chrom_names',
              start.field='START', end.field='END')
    seqlevelsStyle(gr) <- "UCSC"

    anno <- annotatePeak(gr, TxDb=txdb, annoDb=OrgDb, tssRegion=c(-2000,2000))
    df <- as.data.frame(anno)

    write.xlsx(df, "master_table.xlsx", row.names=FALSE)

    EOF
    """
}
