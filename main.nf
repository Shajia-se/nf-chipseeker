#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process chipseeker_annotate {
  tag { sample }

  publishDir "${params.project_folder}/${params.chipseeker_output}/${sample}", mode: 'copy', overwrite: true

  input:
    tuple val(sample), path(peak_file)
    path(gtf)

  output:
    tuple val(sample), path("stats.${sample}.tsv"), emit: stats_tsv

    path("annotated_peaks.${sample}.tsv"),  emit: anno_tsv
    path("annotated_peaks.${sample}.xlsx"), emit: anno_xlsx

    path("annotated.${sample}.pdf"),        emit: anno_pdf
    path("tagMatrix.${sample}.pdf"),        optional: true, emit: tag_pdf

    path("enriched_reactome_all.${sample}.tsv"), optional: true, emit: react_all_tsv
    path("enriched_reactome_sig.${sample}.tsv"), optional: true, emit: react_sig_tsv

  script:
  """
  set -euo pipefail
  mkdir -p tmp
  export TMPDIR=\$PWD/tmp
  export TEMP=\$PWD/tmp
  export TMP=\$PWD/tmp

  cat > run.R <<'RS'
  suppressPackageStartupMessages({
    library(ChIPseeker)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(GenomeInfoDb)
    library(openxlsx)
  })

  sample      <- "${sample}"
  peaks_file  <- "${peak_file}"
  gtf_file    <- "${gtf}"

  annoDb      <- "${params.annoDb}"
  reactome_org<- "${params.reactome_org}"
  do_enrich   <- as.logical("${params.do_enrich}")

  tss_up      <- as.integer("${params.tss_up}")
  tss_down    <- as.integer("${params.tss_down}")
  fdr_cutoff  <- as.numeric("${params.fdr_cutoff}")

  has_reactome <- requireNamespace("ReactomePA", quietly=TRUE)

  if (!file.exists(gtf_file)) stop("GTF not found: ", gtf_file)
  if (!file.exists(peaks_file)) stop("Peak file not found: ", peaks_file)

  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format="gtf")

  normalize_peak_df <- function(D) {
    cn <- colnames(D)
    if (all(c("seqnames","start","end") %in% cn)) return(D)
    if (all(c("chr","start","end") %in% cn)) {
      colnames(D)[colnames(D)=="chr"] <- "seqnames"
      return(D)
    }
    if (all(c("CHR","START","END") %in% cn)) {
      colnames(D)[colnames(D)=="CHR"] <- "seqnames"
      colnames(D)[colnames(D)=="START"] <- "start"
      colnames(D)[colnames(D)=="END"] <- "end"
      return(D)
    }
    NULL
  }

  D <- tryCatch(read.delim(peaks_file, as.is=TRUE, check.names=FALSE, comment.char=""), error=function(e) NULL)
  D <- if (!is.null(D)) normalize_peak_df(D) else NULL

  if (is.null(D)) {
    D0 <- read.table(peaks_file, sep="\t", header=FALSE, as.is=TRUE, comment.char="", quote="")
    if (ncol(D0) < 3) stop("Peak file must have at least 3 columns: ", peaks_file)

    colnames(D0)[1:3] <- c("seqnames","start","end")
    if (ncol(D0) >= 5) colnames(D0)[5] <- "score"
    if (ncol(D0) >= 7) colnames(D0)[7] <- "signalValue"
    if (ncol(D0) >= 8) colnames(D0)[8] <- "pValue"
    if (ncol(D0) >= 9) colnames(D0)[9] <- "qValue"
    D <- D0
  }

  if (!("strand" %in% colnames(D))) D\$strand <- "*"

  gr <- GRanges(
    seqnames = D\$seqnames,
    ranges   = IRanges(as.integer(D\$start), as.integer(D\$end)),
    strand   = D\$strand
  )

  if ("score" %in% colnames(D))       mcols(gr)\$score <- suppressWarnings(as.numeric(D\$score))
  if ("signalValue" %in% colnames(D)) mcols(gr)\$signalValue <- suppressWarnings(as.numeric(D\$signalValue))
  if ("pValue" %in% colnames(D))      mcols(gr)\$pValue <- suppressWarnings(as.numeric(D\$pValue))
  if ("qValue" %in% colnames(D))      mcols(gr)\$qValue <- suppressWarnings(as.numeric(D\$qValue))

  tx_seqs <- seqlevels(txdb)
  tx_has_chr <- any(grepl("^chr", tx_seqs))
  gr_has_chr <- any(grepl("^chr", as.character(seqnames(gr))))

  if (tx_has_chr && !gr_has_chr) {
    suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC")
  } else if (!tx_has_chr && gr_has_chr) {
    suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI")
  }

  common_seqs <- intersect(seqlevels(gr), tx_seqs)
  gr <- GenomeInfoDb::keepSeqlevels(gr, common_seqs, pruning.mode="coarse")

  if (length(gr) == 0) {
    stop("No peaks left after harmonizing seqlevels")
  }

  peakAnno <- annotatePeak(
    gr,
    tssRegion = c(-tss_up, tss_down),
    TxDb = txdb,
    annoDb = annoDb
  )

  peak_df <- as.data.frame(peakAnno)

  write.table(peak_df, file=paste0("annotated_peaks.", sample, ".tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  openxlsx::write.xlsx(peak_df, file=paste0("annotated_peaks.", sample, ".xlsx"), rowNames=FALSE)

  pdf(paste0("annotated.", sample, ".pdf"), width=9, height=5)
    par(mfrow=c(1,2))
    plotAnnoBar(peakAnno)
    plotDistToTSS(peakAnno, title=paste("Distance to TSS -", sample))
  dev.off()

  promoter <- getPromoters(TxDb=txdb, upstream=tss_up, downstream=tss_up)
  tagMatrix <- tryCatch({ getTagMatrix(gr, windows=promoter) }, error=function(e) NULL)
  if (!is.null(tagMatrix) && !is.null(dim(tagMatrix))) {
    pdf(paste0("tagMatrix.", sample, ".pdf"))
      tagHeatmap(tagMatrix)
      plotAvgProf(tagMatrix, xlim=c(-tss_up, tss_up))
    dev.off()
  }

  if (do_enrich && has_reactome && ("geneId" %in% colnames(peak_df))) {
    suppressPackageStartupMessages(library(ReactomePA))
    genes_all <- unique(na.omit(peak_df\$geneId))
    if (length(genes_all) > 0) {
      pw_all <- enrichPathway(genes_all, organism=reactome_org)
      write.table(as.data.frame(pw_all),
        file=paste0("enriched_reactome_all.", sample, ".tsv"),
        sep="\t", quote=FALSE, row.names=FALSE
      )
    }

    if ("qValue" %in% colnames(peak_df)) {
      genes_sig <- unique(na.omit(subset(peak_df, qValue <= fdr_cutoff)\$geneId))
      if (length(genes_sig) > 0) {
        pw_sig <- enrichPathway(genes_sig, organism=reactome_org)
        write.table(as.data.frame(pw_sig),
          file=paste0("enriched_reactome_sig.", sample, ".tsv"),
          sep="\t", quote=FALSE, row.names=FALSE
        )
      }
    }
  }

  peak_df\$merge_id <- paste(peak_df\$seqnames, peak_df\$start, peak_df\$end, sep=":")

  out <- peak_df[, intersect(c("merge_id","annotation","distanceToTSS","score","signalValue","pValue","qValue"), colnames(peak_df)), drop=FALSE]

  for (nm in setdiff(colnames(out), "merge_id")) {
    colnames(out)[colnames(out)==nm] <- paste0(nm, "|", sample)
  }

  write.table(out, file=paste0("stats.", sample, ".tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  RS

  Rscript run.R
  test -s "stats.${sample}.tsv"
  """
}

process chipseeker_master {
  tag "master"

  publishDir "${params.project_folder}/${params.chipseeker_output}", mode: 'copy', overwrite: true

  input:
    path(stats_files)

  output:
    path("annotated_master_table.xlsx"), emit: master_xlsx
    path("annotated_master_table.tsv"), emit: master_tsv

  script:
  '''
  set -euo pipefail
  mkdir -p tmp
  export TMPDIR=$PWD/tmp
  export TEMP=$PWD/tmp
  export TMP=$PWD/tmp

  cat > master.R <<'RS'
  suppressPackageStartupMessages({
    library(openxlsx)
  })

  files <- list.files(pattern="^stats[.].*[.]tsv$", full.names=TRUE)
  if (length(files) == 0) stop("No stats files found")

  M <- read.delim(files[1], as.is=TRUE)
  for (f in files[-1]) {
    D <- read.delim(f, as.is=TRUE)
    M <- merge(M, D, by="merge_id", all=TRUE)
  }

  write.table(M, "annotated_master_table.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  openxlsx::write.xlsx(M, "annotated_master_table.xlsx", rowNames=FALSE)
  RS

  Rscript master.R
  '''
}

process chipseeker_summary {
  tag "summary"

  publishDir "${params.project_folder}/${params.chipseeker_output}", mode: 'copy', overwrite: true

  input:
    path(master_tsv)

  output:
    path("annotation_summary.by_sample.tsv")
    path("annotation_summary.by_sample.pdf")

  script:
  '''
  set -euo pipefail
  mkdir -p tmp
  export TMPDIR=$PWD/tmp
  export TEMP=$PWD/tmp
  export TMP=$PWD/tmp

  cat > summary.R <<'RS'
  D <- read.delim("annotated_master_table.tsv", as.is=TRUE, check.names=FALSE)

  anno_cols <- grep("annotation|", colnames(D), value=TRUE, fixed=TRUE)
  if (length(anno_cols) == 0) {
    write.table(data.frame(sample=character(), annotation=character(), n=integer(), fraction=numeric()),
                "annotation_summary.by_sample.tsv", sep="\\t", quote=FALSE, row.names=FALSE)
    pdf("annotation_summary.by_sample.pdf", width=8, height=5)
    plot.new(); title("No annotation|sample columns found")
    dev.off()
    quit(save="no")
  }

  out <- data.frame(sample=character(), annotation=character(), n=integer(), fraction=numeric(), stringsAsFactors=FALSE)
  for (cn in anno_cols) {
    sample <- sub("annotation|", "", cn, fixed=TRUE)
    vals <- D[[cn]]
    vals <- vals[!is.na(vals) & vals != ""]
    if (length(vals) == 0) next

    # collapse detailed labels for clearer overview
    vals2 <- vals
    vals2[grepl("Promoter", vals2)] <- "Promoter"
    vals2[grepl("Exon", vals2)] <- "Exon"
    vals2[grepl("Intron", vals2)] <- "Intron"
    vals2[grepl("Intergenic", vals2)] <- "Intergenic"
    vals2[grepl("UTR", vals2)] <- "UTR"

    T <- sort(table(vals2), decreasing=TRUE)
    tmp <- data.frame(
      sample = sample,
      annotation = names(T),
      n = as.integer(T),
      fraction = as.numeric(T) / sum(T),
      stringsAsFactors = FALSE
    )
    out <- rbind(out, tmp)
  }

  write.table(out, "annotation_summary.by_sample.tsv", sep="\\t", quote=FALSE, row.names=FALSE)

  pdf("annotation_summary.by_sample.pdf", width=10, height=max(2.5, 1.2 + 0.8 * length(unique(out$sample))))
  if (nrow(out) == 0) {
    plot.new(); title("No annotation summary available")
  } else {
    samples <- unique(out$sample)
    annos <- unique(out$annotation)
    M <- matrix(0, nrow=length(annos), ncol=length(samples), dimnames=list(annos, samples))
    for (i in seq_len(nrow(out))) {
      M[out$annotation[i], out$sample[i]] <- out$fraction[i]
    }
    cols <- c(
      Promoter   = "#9dbfd3",
      Exon       = "#e15759",
      Intron     = "#f28e2b",
      Intergenic = "#6b4c9a",
      UTR        = "#59a14f"
    )
    use_cols <- cols[rownames(M)]
    use_cols[is.na(use_cols)] <- "#bdbdbd"
    par(mar=c(4, 8, 2, 8), xpd=NA)
    barplot(
      t(M),
      horiz=TRUE,
      beside=FALSE,
      col=use_cols,
      border=NA,
      xlim=c(0, 1),
      las=1,
      xlab="Fraction",
      main="Peak Annotation Composition"
    )
    legend("topright", inset=c(-0.18, 0), legend=rownames(M), fill=use_cols, bty="n", cex=0.8)
  }
  dev.off()
  RS

  Rscript summary.R
  '''
}

workflow {
  if (!params.gtf) error "Missing --gtf"
  if (!file(params.gtf).exists()) error "GTF not found: ${params.gtf}"

  def idrPattern = params.idr_peak_pattern ?: "*_idr.sorted.chr.narrowPeak"
  def consensusPattern = params.consensus_peak_pattern ?: "*_consensus.bed"
  def diffbindPattern = params.diffbind_peak_pattern ?: "*.bed"
  def peakSources = (params.chipseeker_peak_sources ?: 'idr,consensus,diffbind')
    .toString()
    .split(',')
    *.trim()
    .findAll { it }
    .unique()
  def selectedPairs = null as Set

  if (params.idr_pairs_csv && file(params.idr_pairs_csv).exists()) {
    selectedPairs = [] as Set
    file(params.idr_pairs_csv).eachLine { line, n ->
      if (n == 1 || !line?.trim()) return
      def cols = line.split(',', -1)*.trim()
      if (cols.size() > 0 && cols[0]) {
        selectedPairs << cols[0]
      }
    }
  }

  ch_gtf = Channel.value(file(params.gtf))

  def peakRows = []

  if (peakSources.contains('idr')) {
    def idrDir = file(params.idr_output)
    assert idrDir.exists() : "idr_output not found: ${params.idr_output}"
    def idrFiles = file("${params.idr_output}").listFiles()?.findAll { f ->
      f.isFile() && f.name ==~ globToRegex(idrPattern)
    } ?: []
    idrFiles.each { f ->
      def sample = f.baseName
        .replaceFirst(/_idr\.sorted\.chr$/, '')
        .replaceFirst(/_idr\.sorted$/, '')
        .replaceFirst(/\.narrowPeak$/, '')
      if (selectedPairs == null || selectedPairs.contains(sample)) {
        peakRows << tuple("idr__${sample}", file(f.toString()))
      }
    }
  }

  if (peakSources.contains('consensus')) {
    def consensusDir = file(params.peak_consensus_output)
    assert consensusDir.exists() : "peak_consensus_output not found: ${params.peak_consensus_output}"
    def consensusFiles = file("${params.peak_consensus_output}").listFiles()?.findAll { f ->
      f.isFile() && f.name ==~ globToRegex(consensusPattern)
    } ?: []
    consensusFiles.each { f ->
      peakRows << tuple("consensus__${f.baseName}", file(f.toString()))
    }
  }

  if (peakSources.contains('diffbind')) {
    def diffbindDir = file(params.diffbind_output)
    assert diffbindDir.exists() : "diffbind_output not found: ${params.diffbind_output}"
    def diffbindFiles = file("${params.diffbind_output}").listFiles()?.findAll { f ->
      f.isFile() && f.name ==~ globToRegex(diffbindPattern)
    } ?: []
    diffbindFiles.each { f ->
      peakRows << tuple("diffbind__${f.baseName}", file(f.toString()))
    }
  }

  Channel
    .fromList(peakRows)
    .ifEmpty { error "No peak files found for ChIPseeker. Check configured source directories and patterns." }
    .set { ch_peak_sets }

  annotated = chipseeker_annotate(ch_peak_sets, ch_gtf)
  stats_paths = annotated.stats_tsv.map { s, f -> f }.collect()
  master = chipseeker_master(stats_paths)
  chipseeker_summary(master.master_tsv)
}

def globToRegex(pattern) {
  '^' + pattern
    .replace('.', '\\.')
    .replace('*', '.*')
    .replace('?', '.') + '$'
}
