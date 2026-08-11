# nf-chipseeker

`nf-chipseeker` annotates peak files and creates peak annotation summary tables/plots.

This module is optional and can consume peaks from IDR, peak consensus, universe peaks, or DiffBind.

## Inputs

Required:

```bash
--gtf /path/to/annotation.gtf
```

Select peak sources:

```bash
--chipseeker_peak_sources idr,consensus_q0.01,consensus_q0.05,universe_q0.05,diffbind
```

Provide only directories required by selected sources:

- `idr` needs `--idr_output`
- consensus/universe sources need `--peak_consensus_output`
- `diffbind` needs `--diffbind_output`

## Outputs

```text
${project_folder}/${chipseeker_output}/
  annotated_master_table.tsv
  annotated_master_table.xlsx
  annotation_summary.by_sample.tsv
  annotation_summary.by_sample.pdf
```

Each peak set also gets its own output folder.

## Run

```bash
nextflow run main.nf -profile hpc \
  --gtf /path/to/annotation.gtf \
  --idr_output /path/to/idr_output \
  --peak_consensus_output /path/to/peak_consensus_output \
  --diffbind_output /path/to/diffbind_output \
  --project_folder /path/to/output_project
```

Actual execution should be tested where Nextflow is installed.
