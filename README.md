# nf-chipseeker

Nextflow DSL2 module for ChIPseeker-based peak annotation using IDR-filtered peaks.

## What This Module Does

This module takes IDR peak files (`*_idr.sorted.chr.narrowPeak`) and performs:

- genomic annotation with ChIPseeker (`annotatePeak`)
- annotation plots (pie + distance to TSS)
- optional tag matrix plot
- optional Reactome pathway enrichment
- per-sample annotation tables (`tsv` and `xlsx`)
- merged master table across all samples

This is intended to run directly after `nf-idr`.

## Inputs

Required:

- `--idr_output`: directory containing IDR peak files
- `--gtf`: GTF annotation file

Optional:

- `--idr_peak_pattern` (default: `*_idr.sorted.chr.narrowPeak`)
- `--annoDb` (default: `org.Mm.eg.db`)
- `--reactome_org` (default: `mouse`)
- `--do_enrich` (default: `true`)
- `--tss_up` / `--tss_down` (default: `2000` / `2000`)
- `--fdr_cutoff` (default: `0.05`)

## Outputs

Output root: `${project_folder}/${chipseeker_output}`

Per sample (`${chipseeker_output}/<sample>/`):

- `annotated_peaks.<sample>.tsv`
- `annotated_peaks.<sample>.xlsx`
- `annotated.<sample>.pdf`
- `tagMatrix.<sample>.pdf` (optional)
- `enriched_reactome_all.<sample>.tsv` (optional)
- `enriched_reactome_sig.<sample>.tsv` (optional)
- `stats.<sample>.tsv`

Merged output (`${chipseeker_output}/`):

- `annotated_master_table.tsv`
- `annotated_master_table.xlsx`
- `annotation_summary.by_sample.tsv`
- `annotation_summary.by_sample.pdf`

## Run

HPC:

```bash
nextflow run main.nf -profile hpc
```

Example with explicit input path:

```bash
nextflow run main.nf -profile hpc \
  --idr_output /ictstr01/groups/idc/projects/uhlenhaut/jiang/pipelines/nf-idr/idr_output \
  --gtf /ictstr01/groups/idc/projects/uhlenhaut/jiang/reference/gtf/mm39.gencode.vM34.annotation.gtf
```

Resume:

```bash
nextflow run main.nf -profile hpc -resume
```

## Notes

- Sample names are derived from peak filenames by stripping suffixes like `_idr.sorted.chr.narrowPeak`.
- Ensure chromosome naming style between peaks and GTF is compatible (`chr1` vs `1`). The module auto-harmonizes styles and keeps shared seqlevels.
