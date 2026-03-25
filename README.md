# nf-chipseeker

Nextflow DSL2 module for ChIPseeker-based peak annotation across IDR, consensus, universe, and DiffBind peak sets.

## What This Module Does

This module can take peak files from:

- `nf-idr`
- `nf-peak-consensus/strict_q0.01`
- `nf-peak-consensus/consensus_q0.05`
- optional universe peak beds from each consensus profile, including:
  - `universe_peaks.bed`
  - `consensus_first_universe_peaks.bed`
  - `union_first_universe_peaks.bed`
- `nf-diffbind`

and performs:

- genomic annotation with ChIPseeker (`annotatePeak`)
- annotation plots (pie + distance to TSS)
- optional tag matrix plot
- optional Reactome pathway enrichment
- per-sample annotation tables (`tsv` and `xlsx`)
- merged master table across all samples

This is intended to run after peak-generation modules are complete.

## Inputs

Required:

- `--gtf`: GTF annotation file

Optional:

- `--idr_peak_pattern` (default: `*_idr.sorted.chr.narrowPeak`)
- `--idr_pairs_csv` (optional; if provided, only annotate pair names listed in `pair_name` column)
- `--annoDb` (default: `org.Mm.eg.db`)
- `--reactome_org` (default: `mouse`)
- `--do_enrich` (default: `true`)
- `--tss_up` / `--tss_down` (default: `3000` / `3000`)
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
  --gtf /ictstr01/groups/idc/projects/uhlenhaut/jiang/reference/gtf/gencode.vM27.primary_assembly.annotation.gtf
```

Restrict to selected IDR pairs:

```bash
nextflow run main.nf -profile hpc \
  --idr_output /ictstr01/groups/idc/projects/uhlenhaut/jiang/pipelines/nf-idr/idr_output \
  --idr_pairs_csv /ictstr01/groups/idc/projects/uhlenhaut/jiang/pipelines/nf-idr/idr_pairs.csv \
  --gtf /ictstr01/groups/idc/projects/uhlenhaut/jiang/reference/gtf/gencode.vM27.primary_assembly.annotation.gtf
```

Resume:

```bash
nextflow run main.nf -profile hpc -resume
```

## Notes

- Sample names are derived from peak filenames by stripping suffixes like `_idr.sorted.chr.narrowPeak`.
- Default peak sources are `idr,consensus_q0.01,consensus_q0.05,diffbind`.
- Additional optional source names supported in `--chipseeker_peak_sources` are `consensus` (alias of `consensus_q0.01`), `universe_q0.01`, and `universe_q0.05`.
- When `universe_q0.05` is enabled, the module now annotates any available q0.05 universe files under `nf-peak-consensus/consensus_q0.05/`, including `universe_peaks.bed`, `consensus_first_universe_peaks.bed`, and `union_first_universe_peaks.bed`.
- Ensure chromosome naming style between peaks and GTF is compatible (`chr1` vs `1`). The module auto-harmonizes styles and keeps shared seqlevels.
