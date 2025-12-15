#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process chipseeker_annotate {
  tag { comp }

  publishDir "${params.chipseeker_output}/${comp}", mode: 'copy', overwrite: true

  input:
    tuple val(comp), path(all_peaks_tsv)
    path(gtf)

  output:
    tuple val(comp), path("stats.${comp}.tsv"), emit: stats_tsv

  script:
    """
    set -euo pipefail

    cat > run.R << 'RS'
    suppressPackageStartupMessages({
      library(ChIPseeker)
      library(GenomicRanges)
      library(GenomicFeatures)
      library(openxlsx)
      library(ggplot2)
    })

    has_reactome <- requireNamespace("ReactomePA", quietly=TRUE)
    has_enrichplot <- requireNamespace("enrichplot", quietly=TRUE)

    comp_name    <- "${comp}"
    peaks_file   <- "${all_peaks_tsv}"
    gtf_file     <- "${gtf}"

    annoDb       <- "${params.annoDb}"          # e.g., "org.Mm.eg.db"
    reactome_org <- "${params.reactome_org}"    # e.g., "mouse"
    do_enrich    <- as.logical("${params.do_enrich}")

    tss_up     <- as.integer("${params.tss_up}")
    tss_down   <- as.integer("${params.tss_down}")
    fdr_cutoff <- as.numeric("${params.fdr_cutoff}")

    # Build TxDb from GTF (mm39 + GENCODE vM34)
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format="gtf")

    message("Reading peaks: ", peaks_file)
    D <- read.delim(peaks_file, as.is=TRUE)

    # Required columns
    req <- c("seqnames","start","end")
    miss <- setdiff(req, colnames(D))
    if (length(miss) > 0) stop("Missing columns in ", peaks_file, ": ", paste(miss, collapse=", "))

    # Optional columns
    if (!("strand" %in% colnames(D))) D$strand <- "*"
    if (!("width"  %in% colnames(D))) D$width  <- D$end - D$start + 1
    for (nm in c("Conc","Fold","FDR")) if (!(nm %in% colnames(D))) D[[nm]] <- NA_real_

    gr <- GRanges(
      seqnames = D$seqnames,
      ranges   = IRanges(D$start, D$end),
      strand   = D$strand
    )
    mcols(gr)$Conc <- D$Conc
    mcols(gr)$Fold <- D$Fold
    mcols(gr)$FDR  <- D$FDR
    mcols(gr)$width <- D$width

    # Annotation
    peakAnno <- annotatePeak(gr, tssRegion=c(-tss_up, tss_down), TxDb=txdb, annoDb=annoDb)
    peak_df  <- as.data.frame(peakAnno)

    # Export annotated peaks
    write.table(peak_df, file=paste0("annotated_peaks.", comp_name, ".tsv"),
                sep="\\t", quote=FALSE, row.names=FALSE)
    openxlsx::write.xlsx(peak_df, file=paste0("annotated_peaks.", comp_name, ".xlsx"), rowNames=FALSE)

    # Annotation plots
    pdf(paste0("annotated.", comp_name, ".pdf"), width=10, height=7)
      plotAnnoPie(peakAnno)
      plotDistToTSS(peakAnno, title="Distribution relative to TSS")
    dev.off()

    # Promoter ± tss_up heatmap/profile (weighted by Conc if available)
    promoter <- getPromoters(TxDb=txdb, upstream=tss_up, downstream=tss_up)
    if (!all(is.na(mcols(gr)$Conc))) {
      tagMatrix <- try(getTagMatrix(gr, windows=promoter, weightCol="Conc"), silent=TRUE)
      if (!inherits(tagMatrix, "try-error") && !is.null(dim(tagMatrix))) {
        pdf(paste0("tagMatrix.", comp_name, ".pdf"))
          tagHeatmap(tagMatrix, xlim=c(-tss_up, tss_up))
          plotAvgProf(tagMatrix, xlim=c(-tss_up, tss_up),
                      xlab="Genomic Region (5'->3')", ylab="Weighted signal (Conc)")
        dev.off()
      }
    }

    # Violin: Fold by feature (only if Fold exists)
    if (!all(is.na(peak_df$Fold))) {
      peaks <- peak_df
      peaks$peak.group <- peaks$annotation
      peaks$peak.group[grepl("Exon",   peaks$peak.group)] <- "Exon"
      peaks$peak.group[grepl("Intron", peaks$peak.group)] <- "Intron"

      pdf(paste0("feature_foldChange_violin.", comp_name, ".pdf"), height=7, width=14)
        p <- ggplot(peaks, aes(x=peak.group, y=Fold, fill=peak.group)) +
          geom_violin() + ggtitle(comp_name) +
          geom_hline(yintercept=0, lwd=0.7, lty=2) +
          stat_summary(geom="crossbar", size=0.1, fun="median")
        print(p)
      dev.off()
    }

    # Enrichment (Reactome only; all + significant)
    if (do_enrich && has_reactome && ("geneId" %in% colnames(peak_df))) {
      suppressPackageStartupMessages(library(ReactomePA))

      genes_all <- unique(na.omit(peak_df$geneId))
      if (length(genes_all) > 0) {
        pw_all <- ReactomePA::enrichPathway(genes_all, organism=reactome_org)
        write.table(as.data.frame(pw_all), file=paste0("enriched_reactome_all.", comp_name, ".tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)
        openxlsx::write.xlsx(as.data.frame(pw_all), file=paste0("enriched_reactome_all.", comp_name, ".xlsx"), rowNames=FALSE)
        if (has_enrichplot && nrow(as.data.frame(pw_all)) > 0) {
          pdf(paste0("enrichment_all.", comp_name, ".pdf"), width=12, height=6)
            print(enrichplot::dotplot(pw_all))
          dev.off()
        }
      }

      if ("FDR" %in% colnames(peak_df)) {
        genes_sig <- unique(na.omit(subset(peak_df, FDR <= fdr_cutoff)$geneId))
        if (length(genes_sig) > 0) {
          pw_sig <- ReactomePA::enrichPathway(genes_sig, organism=reactome_org)
          write.table(as.data.frame(pw_sig), file=paste0("enriched_reactome_sig.", comp_name, ".tsv"),
                      sep="\\t", quote=FALSE, row.names=FALSE)
          openxlsx::write.xlsx(as.data.frame(pw_sig), file=paste0("enriched_reactome_sig.", comp_name, ".xlsx"), rowNames=FALSE)
          if (has_enrichplot && nrow(as.data.frame(pw_sig)) > 0) {
            pdf(paste0("enrichment_sig.", comp_name, ".pdf"), width=12, height=6)
              print(enrichplot::dotplot(pw_sig))
            dev.off()
          }
        }
      }
    }

    # stats file for master merge (merge_id + Conc/Fold/FDR)
    peak_df$merge_id <- paste(peak_df$seqnames, peak_df$start, peak_df$end, sep=":")
    out <- peak_df[, intersect(c("merge_id","Conc","Fold","FDR"), colnames(peak_df))]
    colnames(out)[colnames(out)=="Conc"] <- paste0("Conc|", comp_name)
    colnames(out)[colnames(out)=="Fold"] <- paste0("Fold|", comp_name)
    colnames(out)[colnames(out)=="FDR"]  <- paste0("FDR|",  comp_name)

    write.table(out, file=paste0("stats.", comp_name, ".tsv"),
                sep="\\t", quote=FALSE, row.names=FALSE)
    RS

    Rscript run.R
    """
}

process chipseeker_master {
  tag "master"

  publishDir "${params.chipseeker_output}", mode: 'copy', overwrite: true

  input:
    path(consensus_tsv)
    path(stats_files)
    path(gtf)

  output:
    path("annotated_master_table.xlsx")

  script:
    """
    set -euo pipefail

    cat > master.R << 'RS'
    suppressPackageStartupMessages({
      library(ChIPseeker)
      library(GenomicFeatures)
      library(openxlsx)
    })

    gtf_file   <- "${gtf}"
    cons_file  <- "${consensus_tsv}"
    stats_list <- strsplit("${stats_files.join(' ')}", " ")[[1]]

    annoDb   <- "${params.annoDb}"
    tss_up   <- as.integer("${params.tss_up}")
    tss_down <- as.integer("${params.tss_down}")

    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format="gtf")

    C <- read.delim(cons_file, as.is=TRUE)

    # Auto-detect common DiffBind consensus formats
    if (all(c("chrom_names","START","END") %in% colnames(C))) {
      seqnames <- C$chrom_names; start <- C$START; end <- C$END
    } else if (all(c("seqnames","start","end") %in% colnames(C))) {
      seqnames <- C$seqnames; start <- C$start; end <- C$end
    } else {
      stop("consensus_peaks.tsv columns not recognized. Found: ", paste(colnames(C), collapse=", "))
    }

    gr <- GRanges(seqnames=seqnames, ranges=IRanges(start, end))
    peakAnno <- annotatePeak(gr, tssRegion=c(-tss_up, tss_down), TxDb=txdb, annoDb=annoDb)
    D <- as.data.frame(peakAnno)
    D$merge_id <- paste(D$seqnames, D$start, D$end, sep=":")

    # Merge in each comparison's stats by merge_id
    for (f in stats_list) {
      S <- read.delim(f, as.is=TRUE)
      D <- merge(D, S, by="merge_id", all.x=TRUE)
    }

    D <- D[, c("merge_id", setdiff(colnames(D), "merge_id"))]
    write.xlsx(D, "annotated_master_table.xlsx", rowNames=FALSE)
    RS

    Rscript master.R
    """
}


workflow.onStart {
  if (!params.diffbind_output) error "Missing --diffbind_output"
  if (!params.gtf)         error "Missing --gtf (mm39 GTF)"
}

Channel
  .fromPath("${params.diffbind_output}/all_peaks.*.tsv")
  .ifEmpty { error "No all_peaks.*.tsv found under: ${params.diffbind_output}" }
  .map { f ->
    def comp = f.baseName.replaceFirst(/^all_peaks\./,'')
    tuple(comp, f)
  }
  .set { ch_all_peaks }

Channel
  .fromPath("${params.diffbind_output}/consensus_peaks.tsv")
  .ifEmpty { error "Missing consensus_peaks.tsv under: ${params.diffbind_output}" }
  .set { ch_consensus }

workflow {

  annotated = chipseeker_annotate(ch_all_peaks, file(params.gtf))

  chipseeker_master(ch_consensus, annotated.stats_tsv.collect(), file(params.gtf))
}

