# raccoon

<p align="center">
  <img src="docs/raccoon_logo.png" alt="raccoon logo" width="240" />
</p>

<p align="center"><strong>Rigorous Alignment Curation: Cleanup Of Outliers and Noise</strong></p>

Raccoon is a lightweight toolkit for alignment and phylogenetic QC workflows. It identifies problematic sites (e.g., clustered SNPs, SNPs near Ns/gaps, and frame‑breaking indels) and produces mask files and summaries for downstream analyses.

---

## Contents

- [Use cases](#use-cases)
- [Installation](#installation)
- [Quickstart](#quickstart)
- [CLI usage](#cli-usage)
- [Mask notes](#mask-notes)
- [Example data](#example-data)

## Use cases

- Flag clustered SNPs that may indicate contamination, recombination, or misalignment.
- Detect SNPs adjacent to low-coverage regions (Ns) or gaps.
- Identify frame-breaking indels in coding regions using a GenBank reference.
- Generate mask files to exclude suspect sites prior to phylogenetic or evolutionary analyses.

## Installation

From source:

```bash
pip install .
```

For development (editable install):

```bash
pip install -e .
```

## Quickstart

```bash
raccoon alignment examples/constructed_alignment.fasta -d outdir \
	--genbank examples/constructed_reference.gb --reference-id ref
```

Outputs:

- mask_sites.csv
- alignment_qc_summary.txt

## CLI usage

Show help:

```bash
raccoon --help
```

Alignment QC:

```bash
raccoon alignment <alignment.fasta> -d outdir
```

With a GenBank reference for frame‑break detection:

```bash
raccoon alignment <alignment.fasta> -d outdir \
  --genbank <reference.gb> --reference-id <ref_id>
```

Masking toggles (defaults are enabled):

```bash
raccoon alignment <alignment.fasta> -d outdir \
  --no-mask-n-adjacent --no-mask-gap-adjacent
```

Key alignment options:

- `--n-threshold`: fraction of Ns allowed per sequence before flagging.
- `--cluster-window`: window size (bp) for clustered SNP detection.
- `--cluster-count`: minimum SNPs within a window to flag as clustered.
- `--mask-clustered/--no-mask-clustered`: include/exclude clustered SNPs.
- `--mask-n-adjacent/--no-mask-n-adjacent`: include/exclude SNPs adjacent to Ns.
- `--mask-gap-adjacent/--no-mask-gap-adjacent`: include/exclude SNPs adjacent to gaps.
- `--mask-frame-break/--no-mask-frame-break`: include/exclude frame-breaking indels.

Combine FASTA files:

```bash
raccoon combine a.fasta b.fasta -o combined.fasta
```

With metadata-driven headers:

```bash
raccoon combine a.fasta b.fasta -o combined.fasta \
  --metadata metadata.csv other_metadata.csv --metadata-id-field id \
  --metadata-location-field location --metadata-date-field date \
  --header-separator '|'
```

See full CLI details in [docs/cli.md](docs/cli.md).

## Mask notes

Mask output uses the following note values:

| Note | Meaning |
| --- | --- |
| clustered_snps | Clustered SNPs within the configured window. |
| N_adjacent | SNPs adjacent to an N run within the configured window. |
| gap_adjacent | SNPs adjacent to a gap within the configured window. |
| frame_break | Gap sites that break the CDS frame length. |

## Example data

The [examples](examples) folder includes a constructed alignment and GenBank reference suitable for quick testing:

- [examples/constructed_alignment.fasta](examples/constructed_alignment.fasta)
- [examples/constructed_reference.gb](examples/constructed_reference.gb)