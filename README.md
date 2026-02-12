# raccoon

Raccoon is a lightweight toolkit for alignment and phylogenetic QC workflows. It focuses on identifying problematic sites (e.g., clustered SNPs, SNPs near Ns/gaps, and frame-breaking indels) and producing mask files and summaries for downstream analyses.

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

## CLI usage

Show help:

```bash
raccoon --help
```

Alignment QC:

```bash
raccoon alignment <alignment.fasta> -d outdir
```

With a GenBank reference for frame-break detection:

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

See full CLI details in [docs/cli.md](docs/cli.md).

## Example data

The [examples](examples) folder includes a constructed alignment and GenBank reference suitable for quick testing:

- [examples/constructed_alignment.fasta](examples/constructed_alignment.fasta)
- [examples/constructed_reference.gb](examples/constructed_reference.gb)