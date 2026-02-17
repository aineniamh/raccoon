# raccoon CLI

Use the `raccoon` top-level command with subcommands for different QC tasks.

Examples:

- Alignment QC
  ```bash
  raccoon alignment input_alignment.fasta -d outdir --flag-n
  ```

- Phylogenetic QC
  ```bash
  raccoon phylo --phylogeny mytree --assembly-refs refs.fasta -d outdir --run-apobec
  ```

- Combine FASTA files
  ```bash
  raccoon combine a.fasta b.fasta -o combined.fasta
  ```

Each subcommand has its own help available with `raccoon <subcommand> --help`.

# raccoon CLI

Use the `raccoon` top-level command with subcommands for different QC tasks.

Top-level usage

```bash
raccoon <subcommand> [options]
```

Run `raccoon <subcommand> --help` to see subcommand-specific options.

Alignment subcommand

Purpose: run alignment quality-control checks and produce a mask file and summary.

Basic usage:

```bash
raccoon alignment <alignment.fasta> -d outdir
```

Key options

- `alignment` (positional): path to the input alignment (FASTA)
- `-d, --output-dir`: directory to write outputs (created if missing)
- `--genbank`: optional GenBank file containing CDS/features for frame-breaking indel checks
- `--reference-id`: optional sequence id in the GenBank file used as reference mapping
- `--n-threshold`: fraction of N allowed per sequence before flagged (default: 0.2)
- `--cluster-window`: window size (bp) to search for clustered SNPs (default: 10)
- `--cluster-count`: minimum number of SNPs within the window to be considered clustered (default: 3)
- `--mask-clustered/--no-mask-clustered`: include/exclude clustered SNPs in mask output (default: include)
- `--mask-n-adjacent/--no-mask-n-adjacent`: include/exclude SNPs adjacent to Ns in mask output (default: include)
- `--mask-gap-adjacent/--no-mask-gap-adjacent`: include/exclude SNPs adjacent to gaps in mask output (default: include)
- `--mask-frame-break/--no-mask-frame-break`: include/exclude frame-breaking indels in mask output (default: include)

Behavior and exit codes

- The command will try to create `--output-dir` if it does not exist and will verify it is writable.
- Optional input files (`--genbank`, `--reference-id`) are validated if provided.
- Exit code `0` on success, `1` on validation/file errors, `2` on unexpected failures.

Example

```bash
raccoon alignment data/sequences.fasta -d results/alignment_qc --genbank refs/ref.gb --reference-id NC_000000
```

Example with masking toggles

```bash
raccoon alignment data/sequences.fasta -d results/alignment_qc \
  --genbank refs/ref.gb --reference-id NC_000000 \
  --no-mask-n-adjacent --no-mask-gap-adjacent
```

Phylogenetic subcommand

Purpose: run phylogenetic QC (SNP anomaly checks, apobec analyses, plotting helpers).

Basic usage:

```bash
raccoon phylo --phylogeny treefile.newick --assembly-refs refs.fasta -d outdir
```

Key options

- `--phylogeny`: path to Newick tree file
- `--assembly-refs`: path to assembly/reference FASTA used for mapping
- `-d, --output-dir`: directory to write outputs (created if missing)
- `--outgroup-ids`: comma-separated list of outgroup sequence ids
- `--mask-file`: optional mask CSV file with sites to ignore
- `--height`: plotting height parameter (optional)
- `--run-apobec`: flag to run APOBEC3-specific analyses

Behavior and exit codes

- The command validates `--assembly-refs` and optional mask/tree files, and ensures `--output-dir` is writable.
- Exit code `0` on success, `1` on validation/file errors, `2` on unexpected failures.

Example

```bash
raccoon phylo --phylogeny trees/rep.tree --assembly-refs data/refs.fasta -d results/phylo --run-apobec
```

Combine subcommand

Purpose: combine one or more FASTA files into a single upper-case, unwrapped FASTA, with optional metadata-driven header harmonisation.

Basic usage:

```bash
raccoon combine a.fasta b.fasta -o combined.fasta
```

Header harmonisation using metadata:

```bash
raccoon combine a.fasta b.fasta -o combined.fasta \
  --metadata metadata.csv \
  --metadata-id-field id \
  --metadata-location-field location \
  --metadata-date-field date \
  --header-separator '|'
```

Key options

- `inputs` (positional): one or more input FASTA files
- `-o, --output`: output FASTA file (use `-` for stdout)
- `--metadata`: one or more metadata CSVs used to harmonise headers
- `--metadata-delimiter`: metadata delimiter (default: `,`)
- `--metadata-id-field`: metadata id column (default: `id`)
- `--metadata-location-field`: metadata location column (default: `location`)
- `--metadata-date-field`: metadata date column (default: `date`)
- `--header-separator`: header separator (default: `|`)
