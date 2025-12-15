#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * nf-chipseeker (DiffBind export -> ChIPseeker)
 * Inputs under --diffbind_output:
 *   - all_peaks.<comp>.tsv
 *   - consensus_peaks.tsv
 */

process chipseeker_annotate {
  tag { sample }

  publishDir "${params.chipseeker_output}/${sample}", mode: 'copy', overwrite: true

  input:
    tuple val(sample), path(all_peaks_tsv)
    path(gtf)

  output:
    tuple val(sample), path("stats.${sample}.tsv"), emit: stats_tsv

    path("annotated_peaks.${sample}.tsv"),  emit: anno_tsv
    path("annotated_peaks.${sample}.xlsx"), emit: anno_xlsx

    path("annotated.${sample}.pdf"),        emit: anno_pdf
    path("tagMatrix.${sample}.pdf"),        optional: true, emit: tag_pdf
    path("feature_foldChange_violin.${sample}.pdf"), optional: true, emit: violin_pdf

    path("enriched_reactome_all.${sample}.tsv"), optional: true, emit: react_all_tsv
    path("enriched_reactome_sig.${sample}.tsv"), optional: true, emit: react_sig_tsv

  script:
  """
  set -euo pipefail

  SAMPLE="${sample}"
  PEAKS="${all_peaks_tsv}"
  GTF="${gtf}"

  ANNODB="${params.annoDb}"
  REACTOME_ORG="${params.reactome_org}"
  DO_ENRICH="${params.do_enrich}"

  TSS_UP="${params.tss_up}"
  TSS_DOWN="${params.tss_down}"
  FDR_CUTOFF="${params.fdr_cutoff}"

  cat > run.R <<'RS'
  suppressPackageStartupMessages({
    library(ChIPseeker)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(GenomeInfoDb)
    library(openxlsx)
    library(ggplot2)
  })

  sample     <- "${sample}"
  peaks_file <- "${all_peaks_tsv}"
  gtf_file   <- "${gtf}"

  annoDb      <- "${params.annoDb}"
  reactome_org<- "${params.reactome_org}"
  do_enrich   <- as.logical("${params.do_enrich}")

  tss_up      <- as.integer("${params.tss_up}")
  tss_down    <- as.integer("${params.tss_down}")
  fdr_cutoff  <- as.numeric("${params.fdr_cutoff}")

  has_reactome <- requireNamespace("ReactomePA", quietly=TRUE)

  message("GTF: ", gtf_file)
  if (!file.exists(gtf_file)) stop("GTF not found: ", gtf_file)

  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format="gtf")

  message("Reading peaks: ", peaks_file)
  D <- read.delim(peaks_file, as.is=TRUE)

  req <- c("seqnames","start","end")
  miss <- setdiff(req, colnames(D))
  if (length(miss) > 0) stop("Missing columns in peaks file: ", paste(miss, collapse=", "))

  if (!("strand" %in% colnames(D))) D\$strand <- "*"
  if (!("width"  %in% colnames(D))) D\$width  <- D\$end - D\$start + 1
  for (nm in c("Conc","Fold","FDR")) if (!(nm %in% colnames(D))) D[[nm]] <- NA_real_

  gr <- GRanges(
    seqnames = D\$seqnames,
    ranges   = IRanges(D\$start, D\$end),
    strand   = D\$strand
  )
  mcols(gr)\$Conc  <- D\$Conc
  mcols(gr)\$Fold  <- D\$Fold
  mcols(gr)\$FDR   <- D\$FDR
  mcols(gr)\$width <- D\$width

  ## Harmonize seqlevels style
  tx_seqs <- seqlevels(txdb)
  tx_has_chr <- any(grepl("^chr", tx_seqs))
  gr_has_chr <- any(grepl("^chr", as.character(seqnames(gr))))

  if (tx_has_chr && !gr_has_chr) {
    suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC")
  } else if (!tx_has_chr && gr_has_chr) {
    suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI")
  }

  # keep only shared seqlevels (THIS is the fix)
  common_seqs <- intersect(seqlevels(gr), tx_seqs)
  gr <- GenomeInfoDb::keepSeqlevels(gr, common_seqs, pruning.mode="coarse")

  if (length(gr) == 0) {
    stop(
      "No peaks left after harmonizing seqlevels. ",
      "Example TxDb seqlevels: ", paste(head(tx_seqs, 10), collapse=", "),
      " | Example peak seqlevels: ", paste(head(seqlevels(gr), 10), collapse=", ")
    )
  }

  ## Annotation
  peakAnno <- annotatePeak(
    gr,
    tssRegion = c(-tss_up, tss_down),
    TxDb = txdb,
    annoDb = annoDb
  )

  peak_df <- as.data.frame(peakAnno)

  ## Outputs
  write.table(peak_df, file=paste0("annotated_peaks.", sample, ".tsv"),
              sep="\\t", quote=FALSE, row.names=FALSE)
  openxlsx::write.xlsx(peak_df, file=paste0("annotated_peaks.", sample, ".xlsx"),
                       rowNames=FALSE)

  pdf(paste0("annotated.", sample, ".pdf"), width=10, height=7)
    plotAnnoPie(peakAnno)
    plotDistToTSS(peakAnno)
  dev.off()

  ## Tag matrix + profile
  promoter <- getPromoters(TxDb=txdb, upstream=tss_up, downstream=tss_up)

  tagMatrix <- tryCatch({
    if (!all(is.na(mcols(gr)\$Conc))) {
      getTagMatrix(gr, windows=promoter, weightCol="Conc")
    } else {
      getTagMatrix(gr, windows=promoter)
    }
  }, error=function(e) NULL)

  if (!is.null(tagMatrix) && !is.null(dim(tagMatrix))) {
    pdf(paste0("tagMatrix.", sample, ".pdf"))
      ## tagHeatmap signature differs across versions; do not pass xlim
      tagHeatmap(tagMatrix)
      plotAvgProf(tagMatrix, xlim=c(-tss_up, tss_up))
    dev.off()
  }

  ## Violin by feature
  if ("Fold" %in% colnames(peak_df) && !all(is.na(peak_df\$Fold))) {
    peaks <- peak_df
    peaks\$peak.group <- peaks\$annotation
    peaks\$peak.group[grepl("Exon", peaks\$peak.group)]   <- "Exon"
    peaks\$peak.group[grepl("Intron", peaks\$peak.group)] <- "Intron"

    pdf(paste0("feature_foldChange_violin.", sample, ".pdf"), width=14, height=7)
      print(
        ggplot(peaks, aes(x=peak.group, y=Fold, fill=peak.group)) +
          geom_violin() +
          geom_hline(yintercept=0, lty=2) +
          stat_summary(fun=median, geom="crossbar", width=0.3) +
          ggtitle(sample)
      )
    dev.off()
  }

  ## Reactome enrichment (optional)
  if (do_enrich && has_reactome && ("geneId" %in% colnames(peak_df))) {
    suppressPackageStartupMessages(library(ReactomePA))
    genes_all <- unique(na.omit(peak_df\$geneId))

    if (length(genes_all) > 0) {
      pw_all <- enrichPathway(genes_all, organism=reactome_org)
      write.table(as.data.frame(pw_all),
        file=paste0("enriched_reactome_all.", sample, ".tsv"),
        sep="\\t", quote=FALSE, row.names=FALSE
      )
    }

    if ("FDR" %in% colnames(peak_df)) {
      genes_sig <- unique(na.omit(subset(peak_df, FDR <= fdr_cutoff)\$geneId))
      if (length(genes_sig) > 0) {
        pw_sig <- enrichPathway(genes_sig, organism=reactome_org)
        write.table(as.data.frame(pw_sig),
          file=paste0("enriched_reactome_sig.", sample, ".tsv"),
          sep="\\t", quote=FALSE, row.names=FALSE
        )
      }
    }
  }

  ## Stats for master merge
  peak_df\$merge_id <- paste(peak_df\$seqnames, peak_df\$start, peak_df\$end, sep=":")
  out <- peak_df[, intersect(c("merge_id","Conc","Fold","FDR"), colnames(peak_df)), drop=FALSE]

  if ("Conc" %in% colnames(out)) colnames(out)[colnames(out)=="Conc"] <- paste0("Conc|", sample)
  if ("Fold" %in% colnames(out)) colnames(out)[colnames(out)=="Fold"] <- paste0("Fold|", sample)
  if ("FDR"  %in% colnames(out)) colnames(out)[colnames(out)=="FDR"]  <- paste0("FDR|",  sample)

  write.table(out, file=paste0("stats.", sample, ".tsv"),
              sep="\\t", quote=FALSE, row.names=FALSE)
RS

  Rscript run.R

  ## Rename outputs to match declared names
  mv -f "annotated_peaks.${sample}.tsv"  "annotated_peaks.${sample}.tsv"  2>/dev/null || true
  mv -f "annotated_peaks.${sample}.xlsx" "annotated_peaks.${sample}.xlsx" 2>/dev/null || true
  mv -f "annotated.${sample}.pdf"        "annotated.${sample}.pdf"        2>/dev/null || true
  mv -f "tagMatrix.${sample}.pdf"        "tagMatrix.${sample}.pdf"        2>/dev/null || true
  mv -f "feature_foldChange_violin.${sample}.pdf" "feature_foldChange_violin.${sample}.pdf" 2>/dev/null || true
  mv -f "enriched_reactome_all.${sample}.tsv" "enriched_reactome_all.${sample}.tsv" 2>/dev/null || true
  mv -f "enriched_reactome_sig.${sample}.tsv" "enriched_reactome_sig.${sample}.tsv" 2>/dev/null || true

  ## Guarantee stats output exists
  test -s "stats.${sample}.tsv"
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
    path("annotated_master_table.tsv")

  script:
  """
  set -euo pipefail

  ## bring stats into cwd for simple list.files()
  cp -f ${stats_files.join(' ')} . || true

  cat > master.R <<'RS'
  suppressPackageStartupMessages({
    library(ChIPseeker)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(GenomeInfoDb)
    library(openxlsx)
  })

  cons_file <- "${consensus_tsv}"
  gtf_file  <- "${gtf}"
  annoDb    <- "${params.annoDb}"
  tss_up    <- as.integer("${params.tss_up}")
  tss_down  <- as.integer("${params.tss_down}")

  if (!file.exists(gtf_file)) stop("GTF not found: ", gtf_file)
  if (!file.exists(cons_file)) stop("Consensus not found: ", cons_file)

  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format="gtf")

  C <- read.delim(cons_file, as.is=TRUE)
  cn <- colnames(C)

  ## Support multiple common consensus formats
  if (all(c("CHR","START","END") %in% cn)) {
    seqnames <- C\$CHR; start <- C\$START; end <- C\$END
  } else if (all(c("chrom_names","START","END") %in% cn)) {
    seqnames <- C\$chrom_names; start <- C\$START; end <- C\$END
  } else if (all(c("seqnames","start","end") %in% cn)) {
    seqnames <- C\$seqnames; start <- C\$start; end <- C\$end
  } else if (all(c("seqnames","START","END") %in% cn)) {
    seqnames <- C\$seqnames; start <- C\$START; end <- C\$END
  } else if (all(c("chr","start","end") %in% cn)) {
    seqnames <- C\$chr; start <- C\$start; end <- C\$end
  } else {
    stop("Unrecognized consensus_peaks.tsv format. Columns: ", paste(cn, collapse=", "))
  }

  gr <- GRanges(seqnames=seqnames, ranges=IRanges(start, end))

  tx_seqs <- seqlevels(txdb)
  tx_has_chr <- any(grepl("^chr", tx_seqs))
  gr_has_chr <- any(grepl("^chr", as.character(seqnames(gr))))

  if (tx_has_chr && !gr_has_chr) {
    suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC")
  } else if (!tx_has_chr && gr_has_chr) {
    suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI")
  }

  # keep only shared seqlevels (THIS is the fix)
  common_seqs <- intersect(seqlevels(gr), tx_seqs)
  gr <- GenomeInfoDb::keepSeqlevels(gr, common_seqs, pruning.mode="coarse")

  if (length(gr) == 0) {
    stop(
      "No peaks left after harmonizing seqlevels. ",
      "Example TxDb seqlevels: ", paste(head(tx_seqs, 10), collapse=", "),
      " | Example peak seqlevels: ", paste(head(seqlevels(gr), 10), collapse=", ")
    )
  }

  peakAnno <- annotatePeak(gr, tssRegion=c(-tss_up, tss_down), TxDb=txdb, annoDb=annoDb)
  D <- as.data.frame(peakAnno)
  D\$merge_id <- paste(D\$seqnames, D\$start, D\$end, sep=":")

  stats_list <- list.files(full.names = TRUE)
  stats_list <- stats_list[ startsWith(basename(stats_list), "stats.") &
                           endsWith(basename(stats_list), ".tsv") ]
  if (length(stats_list) == 0) stop("No stats.*.tsv found in workdir.")

  for (f in stats_list) {
    S <- read.delim(f, as.is=TRUE)
    D <- merge(D, S, by="merge_id", all.x=TRUE)
  }

  write.table(D, "annotated_master_table.tsv", sep="\\t", quote=FALSE, row.names=FALSE)
  openxlsx::write.xlsx(D, "annotated_master_table.xlsx", rowNames=FALSE)
  
  RS
  Rscript master.R
  """
}

workflow {

  if (!params.diffbind_output) error "Missing --diffbind_output"
  if (!params.gtf)             error "Missing --gtf"
  if (!file(params.gtf).exists()) error "GTF not found: ${params.gtf}"

  ch_gtf = Channel.value( file(params.gtf) )

  Channel
    .fromPath("${params.diffbind_output}/all_peaks.*.tsv")
    .ifEmpty { error "No all_peaks.*.tsv found under: ${params.diffbind_output}" }
    .map { f -> tuple( f.baseName.replaceFirst(/^all_peaks\\./,''), f ) }
    .set { ch_all_peaks }

  Channel
    .fromPath("${params.diffbind_output}/consensus_peaks.tsv")
    .ifEmpty { error "Missing consensus_peaks.tsv under: ${params.diffbind_output}" }
    .set { ch_consensus }

  annotated = chipseeker_annotate(ch_all_peaks, ch_gtf)

  stats_paths = annotated.stats_tsv.map { s, f -> f }.collect()
  chipseeker_master(ch_consensus, stats_paths, ch_gtf)
}
