# UMI-tools-rs

Independent Rust reimplementation of the UMI-tools `dedup` command.

This is not the original UMI-tools package. The original project is
[CGATOxford/UMI-tools](https://github.com/CGATOxford/UMI-tools), which remains
the behavioral reference for this crate. This Rust crate currently supports
only `dedup`.

This reimplementation was produced by
[OpenAI Codex](https://openai.com/codex/) and
[Claude Code](https://www.anthropic.com/product/claude-code).

## Build

```bash
cargo build --release
```

The optimized binary is written to:

```bash
target/release/umi-tools-rs
```

## Usage

The command is shaped to match UMI-tools:

```bash
umi-tools-rs dedup \
  --stdin input.bam \
  --stdout dedup.bam \
  --method directional
```

SAM output:

```bash
umi-tools-rs dedup -I input.bam -S dedup.sam --out-sam
```

Tag-based UMI extraction:

```bash
umi-tools-rs dedup -I input.bam -S dedup.bam \
  --extract-umi-method=tag \
  --umi-tag=RX
```

Paired-end deduplication:

```bash
umi-tools-rs dedup -I input.bam -S dedup.bam --paired
```

By default, output is coordinate-sorted with `samtools sort`, matching UMI-tools behavior. Use `--no-sort-output` to write records directly.

## Implemented

- UMI grouping methods: `unique`, `percentile`, `cluster`, `adjacency`, `directional`
- UMI extraction from read names, BAM tags, and `UMI_`/`CELL_` read-name fields
- `--per-cell`, `--per-gene --gene-tag`, and `--per-contig`
- `--gene-transcript-map` metacontig mode for indexed BAM/CRAM inputs
- `--ignore-umi`, `--ignore-tlen`, `--read-length`, `--spliced-is-unique`
- `--mapping-quality`, `--chrom`, `--subset`, `--random-seed`
- `--output-stats` with `_per_umi_per_position.tsv`, `_per_umi.tsv`, and
  `_edit_distance.tsv` outputs
- Paired-end read1 deduplication with second-pass mate recovery
- Hidden UMI whitelist filtering used by UMI-tools tests
- SAM/BAM/CRAM reading and writing via `rust-htslib`

## Compatibility Notes

The observed stats tables match the Python fixtures for deterministic count
columns. The `_edit_distance.tsv` null-model columns are sampled with Rust's RNG
instead of NumPy's RNG, so they are reproducible with `--random-seed` but are
not bit-for-bit identical to UMI-tools.

## License

MIT.
