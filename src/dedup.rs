use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{Context, Result, bail};
use clap::{Args, ValueEnum};
use indexmap::IndexMap;
use rand::distributions::{Distribution, WeightedIndex};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use regex::Regex;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{self, CompressionLevel, Format, Header, Read, Record, Writer};
use tempfile::TempDir;

use crate::network::{ClusterMethod, cluster_umis};

#[derive(Debug, Args)]
pub struct DedupArgs {
    #[arg(short = 'I', long = "stdin", value_name = "FILE")]
    pub input: PathBuf,

    #[arg(short = 'S', long = "stdout", value_name = "FILE")]
    pub output: Option<PathBuf>,

    #[arg(short = 'L', long = "log", value_name = "FILE")]
    pub log: Option<PathBuf>,

    #[arg(long = "log2stderr", default_value_t = false)]
    pub log2stderr: bool,

    #[arg(long = "in-format", value_enum)]
    pub in_format: Option<AlignmentFormat>,

    #[arg(short = 'i', long = "in-sam", default_value_t = false)]
    pub in_sam: bool,

    #[arg(long = "out-format", value_enum)]
    pub out_format: Option<AlignmentFormat>,

    #[arg(short = 'o', long = "out-sam", default_value_t = false)]
    pub out_sam: bool,

    #[arg(long = "reference-filename", alias = "reference-file")]
    pub reference_filename: Option<PathBuf>,

    #[arg(long = "input-options")]
    pub input_options: Option<String>,

    #[arg(long = "output-options")]
    pub output_options: Option<String>,

    #[arg(long = "temp-dir")]
    pub temp_dir: Option<PathBuf>,

    #[arg(long = "compresslevel")]
    pub compresslevel: Option<u32>,

    #[arg(long = "threads", short = '@', default_value_t = 0)]
    pub threads: usize,

    #[arg(long = "method", value_enum, default_value_t = ClusterMethod::Directional)]
    pub method: ClusterMethod,

    #[arg(long = "edit-distance-threshold", default_value_t = 1)]
    pub threshold: usize,

    #[arg(long = "spliced-is-unique", default_value_t = false)]
    pub spliced: bool,

    #[arg(long = "soft-clip-threshold", default_value_t = 4.0)]
    pub soft_clip_threshold: f64,

    #[arg(long = "read-length", default_value_t = false)]
    pub read_length: bool,

    #[arg(long = "extract-umi-method", value_enum, default_value_t = ExtractUmiMethod::ReadId)]
    pub get_umi_method: ExtractUmiMethod,

    #[arg(long = "umi-separator", default_value = "_")]
    pub umi_sep: String,

    #[arg(long = "umi-tag", default_value = "RX")]
    pub umi_tag: String,

    #[arg(long = "umi-tag-split")]
    pub umi_tag_split: Option<String>,

    #[arg(long = "umi-tag-delimiter")]
    pub umi_tag_delim: Option<String>,

    #[arg(long = "cell-tag")]
    pub cell_tag: Option<String>,

    #[arg(long = "cell-tag-split", default_value = "-")]
    pub cell_tag_split: String,

    #[arg(long = "cell-tag-delimiter")]
    pub cell_tag_delim: Option<String>,

    #[arg(long = "filter-umi", default_value_t = false, hide = true)]
    pub filter_umi: bool,

    #[arg(long = "umi-whitelist", hide = true)]
    pub umi_whitelist: Option<PathBuf>,

    #[arg(long = "umi-whitelist-paired", hide = true)]
    pub umi_whitelist_paired: Option<PathBuf>,

    #[arg(long = "per-gene", default_value_t = false)]
    pub per_gene: bool,

    #[arg(long = "gene-tag")]
    pub gene_tag: Option<String>,

    #[arg(long = "assigned-status-tag")]
    pub assigned_tag: Option<String>,

    #[arg(long = "skip-tags-regex", default_value = "^(__|Unassigned)")]
    pub skip_regex: String,

    #[arg(long = "per-contig", default_value_t = false)]
    pub per_contig: bool,

    #[arg(long = "gene-transcript-map")]
    pub gene_transcript_map: Option<PathBuf>,

    #[arg(long = "per-cell", default_value_t = false)]
    pub per_cell: bool,

    #[arg(
        long = "buffer-whole-contig",
        alias = "whole-contig",
        default_value_t = false
    )]
    pub whole_contig: bool,

    #[arg(long = "multimapping-detection-method", value_enum)]
    pub detection_method: Option<DetectionMethod>,

    #[arg(long = "mapping-quality", default_value_t = 0)]
    pub mapping_quality: u8,

    #[arg(long = "output-unmapped", default_value_t = false, hide = true)]
    pub output_unmapped: bool,

    #[arg(long = "ignore-umi", default_value_t = false)]
    pub ignore_umi: bool,

    #[arg(long = "ignore-tlen", default_value_t = false)]
    pub ignore_tlen: bool,

    #[arg(long = "chrom")]
    pub chrom: Option<String>,

    #[arg(long = "subset")]
    pub subset: Option<f64>,

    #[arg(long = "paired", default_value_t = false)]
    pub paired: bool,

    #[arg(long = "no-sort-output", default_value_t = false)]
    pub no_sort_output: bool,

    #[arg(long = "unmapped-reads", value_enum, default_value_t = ReadHandling::Discard)]
    pub unmapped_reads: ReadHandling,

    #[arg(long = "chimeric-pairs", value_enum, default_value_t = ReadHandling::Use)]
    pub chimeric_pairs: ReadHandling,

    #[arg(long = "unpaired-reads", value_enum, default_value_t = ReadHandling::Use)]
    pub unpaired_reads: ReadHandling,

    #[arg(long = "random-seed")]
    pub random_seed: Option<u64>,

    #[arg(long = "output-stats")]
    pub output_stats: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum AlignmentFormat {
    Sam,
    Bam,
    Cram,
}

impl AlignmentFormat {
    fn htslib(self) -> Format {
        match self {
            AlignmentFormat::Sam => Format::Sam,
            AlignmentFormat::Bam => Format::Bam,
            AlignmentFormat::Cram => Format::Cram,
        }
    }

    fn samtools_name(self) -> &'static str {
        match self {
            AlignmentFormat::Sam => "SAM",
            AlignmentFormat::Bam => "BAM",
            AlignmentFormat::Cram => "CRAM",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum ExtractUmiMethod {
    #[value(name = "read_id", alias = "read-id")]
    ReadId,
    Tag,
    Umis,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum DetectionMethod {
    #[value(name = "NH", alias = "nh")]
    NH,
    #[value(name = "X0", alias = "x0")]
    X0,
    #[value(name = "XT", alias = "xt")]
    XT,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum ReadHandling {
    Discard,
    Use,
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
enum PositionKey {
    Coord { tid: i32, pos: i64 },
    Gene(String),
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
enum GroupInner {
    Coord {
        reverse: bool,
        spliced_offset: i64,
        tlen: i64,
        read_len: usize,
    },
    Gene(String),
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct GroupKey {
    inner: GroupInner,
    cell: Option<Vec<u8>>,
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct MateKey {
    qname: Vec<u8>,
    tid: i32,
    pos: i64,
}

#[derive(Clone)]
struct BundleEntry {
    read: Record,
    count: usize,
    ties_seen: u64,
}

type Bundle = IndexMap<Vec<u8>, BundleEntry>;

struct BundleBuilder {
    reads: BTreeMap<PositionKey, BTreeMap<GroupKey, Bundle>>,
    last_tid: Option<i32>,
    last_start: i64,
    events: HashMap<String, u64>,
    n_input: u64,
}

#[derive(Clone, Debug)]
struct SelectedUmi {
    parent: Vec<u8>,
    count: usize,
}

struct StatsCollector {
    prefix: PathBuf,
    pre_rows: Vec<(Vec<u8>, usize)>,
    post_rows: Vec<(Vec<u8>, usize)>,
    pre_distances: Vec<f64>,
    post_distances: Vec<f64>,
    pre_null_distances: Vec<f64>,
    post_null_distances: Vec<f64>,
    random_umis: RandomUmiSampler,
}

struct RandomUmiSampler {
    umis: Vec<Vec<u8>>,
    distribution: Option<WeightedIndex<usize>>,
    rng: StdRng,
}

struct Logger {
    file: Option<BufWriter<File>>,
}

impl Logger {
    fn new(path: Option<&Path>) -> Result<Self> {
        let file = match path {
            Some(path) => Some(BufWriter::new(
                File::create(path).with_context(|| format!("failed to create log {path:?}"))?,
            )),
            None => None,
        };

        Ok(Self { file })
    }

    fn info(&mut self, message: impl AsRef<str>) -> Result<()> {
        let message = message.as_ref();
        if let Some(file) = self.file.as_mut() {
            writeln!(file, "{message}")?;
        } else {
            eprintln!("{message}");
        }
        Ok(())
    }
}

impl StatsCollector {
    fn new(prefix: PathBuf, args: &DedupArgs, target_names: &[String]) -> Result<Self> {
        Ok(Self {
            prefix,
            pre_rows: Vec::new(),
            post_rows: Vec::new(),
            pre_distances: Vec::new(),
            post_distances: Vec::new(),
            pre_null_distances: Vec::new(),
            post_null_distances: Vec::new(),
            random_umis: RandomUmiSampler::from_input(args, target_names)?,
        })
    }

    fn record_bundle(&mut self, counts: &IndexMap<Vec<u8>, usize>, selected: &[SelectedUmi]) {
        let pre_umis: Vec<Vec<u8>> = counts.keys().cloned().collect();
        let post_umis: Vec<Vec<u8>> = selected.iter().map(|item| item.parent.clone()).collect();

        self.pre_distances.push(average_umi_distance(&pre_umis));
        self.post_distances.push(average_umi_distance(&post_umis));

        let random_pre = self.random_umis.sample(pre_umis.len());
        let random_post = self.random_umis.sample(post_umis.len());
        self.pre_null_distances
            .push(average_umi_distance(&random_pre));
        self.post_null_distances
            .push(average_umi_distance(&random_post));

        self.pre_rows
            .extend(counts.iter().map(|(umi, count)| (umi.clone(), *count)));
        self.post_rows.extend(
            selected
                .iter()
                .map(|item| (item.parent.clone(), item.count)),
        );
    }

    fn write(&self, method: ClusterMethod) -> Result<()> {
        self.write_per_umi_per_position()?;
        self.write_per_umi()?;
        self.write_edit_distance(method)?;
        Ok(())
    }

    fn write_per_umi_per_position(&self) -> Result<()> {
        let pre_counts = count_instances_by_count(&self.pre_rows);
        let post_counts = count_instances_by_count(&self.post_rows);
        let keys: BTreeSet<usize> = pre_counts
            .keys()
            .chain(post_counts.keys())
            .copied()
            .collect();

        let path = suffixed_path(&self.prefix, "_per_umi_per_position.tsv");
        let mut out = BufWriter::new(
            File::create(&path).with_context(|| format!("failed to create {}", path.display()))?,
        );
        writeln!(out, "counts\tinstances_pre\tinstances_post")?;
        for count in keys {
            writeln!(
                out,
                "{count}\t{}\t{}",
                pre_counts.get(&count).copied().unwrap_or(0),
                post_counts.get(&count).copied().unwrap_or(0)
            )?;
        }
        Ok(())
    }

    fn write_per_umi(&self) -> Result<()> {
        let pre = aggregate_umi_rows(&self.pre_rows);
        let post = aggregate_umi_rows(&self.post_rows);
        let path = suffixed_path(&self.prefix, "_per_umi.tsv");
        let mut out = BufWriter::new(
            File::create(&path).with_context(|| format!("failed to create {}", path.display()))?,
        );

        writeln!(
            out,
            "UMI\tmedian_counts_pre\ttimes_observed_pre\ttotal_counts_pre\tmedian_counts_post\ttimes_observed_post\ttotal_counts_post"
        )?;
        for (umi, pre_agg) in pre {
            let post_agg = post.get(&umi).copied().unwrap_or_default();
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                String::from_utf8_lossy(&umi),
                pre_agg.median_counts,
                pre_agg.times_observed,
                pre_agg.total_counts,
                post_agg.median_counts,
                post_agg.times_observed,
                post_agg.total_counts
            )?;
        }
        Ok(())
    }

    fn write_edit_distance(&self, method: ClusterMethod) -> Result<()> {
        let max_ed = self
            .pre_distances
            .iter()
            .chain(self.post_distances.iter())
            .chain(self.pre_null_distances.iter())
            .chain(self.post_null_distances.iter())
            .copied()
            .fold(-1.0_f64, f64::max)
            .floor() as i32;
        let bins: Vec<i32> = (-1..=(max_ed + 1)).collect();
        let unique = tally_digitized(&self.pre_distances, &bins);
        let unique_null = tally_digitized(&self.pre_null_distances, &bins);
        let post = tally_digitized(&self.post_distances, &bins);
        let post_null = tally_digitized(&self.post_null_distances, &bins);
        let method_name = cluster_method_name(method);

        let path = suffixed_path(&self.prefix, "_edit_distance.tsv");
        let mut out = BufWriter::new(
            File::create(&path).with_context(|| format!("failed to create {}", path.display()))?,
        );
        writeln!(
            out,
            "unique\tunique_null\t{method_name}\t{method_name}_null\tedit_distance"
        )?;
        for (idx, label_value) in bins.iter().enumerate() {
            let label = if idx == 0 {
                "Single_UMI".to_string()
            } else {
                label_value.to_string()
            };
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}",
                unique.get(idx).copied().unwrap_or(0),
                unique_null.get(idx).copied().unwrap_or(0),
                post.get(idx).copied().unwrap_or(0),
                post_null.get(idx).copied().unwrap_or(0),
                label
            )?;
        }
        Ok(())
    }
}

impl RandomUmiSampler {
    fn from_input(args: &DedupArgs, target_names: &[String]) -> Result<Self> {
        let mut reader = open_reader(args)?;
        let mut record = Record::new();
        let mut counts: IndexMap<Vec<u8>, usize> = IndexMap::new();

        while let Some(read_result) = reader.read(&mut record) {
            read_result.context("failed to read input while collecting UMI stats")?;
            if record.is_unmapped() || record.is_last_in_template() {
                continue;
            }
            if !should_visit_chrom(&record, target_names, args.chrom.as_deref()) {
                continue;
            }
            if let Ok((umi, _)) = extract_barcodes(&record, args) {
                *counts.entry(umi).or_insert(0) += 1;
            }
        }

        let umis: Vec<Vec<u8>> = counts.keys().cloned().collect();
        let weights: Vec<usize> = counts.values().copied().collect();
        let distribution = if weights.is_empty() {
            None
        } else {
            Some(WeightedIndex::new(weights).context("failed to build UMI sampler")?)
        };
        let rng = match args.random_seed {
            Some(seed) => StdRng::seed_from_u64(seed ^ 0x5eed_5eed_u64),
            None => StdRng::from_entropy(),
        };

        Ok(Self {
            umis,
            distribution,
            rng,
        })
    }

    fn sample(&mut self, n: usize) -> Vec<Vec<u8>> {
        let Some(distribution) = self.distribution.as_ref() else {
            return Vec::new();
        };

        (0..n)
            .map(|_| self.umis[distribution.sample(&mut self.rng)].clone())
            .collect()
    }
}

pub fn run(args: DedupArgs) -> Result<()> {
    validate_args(&args)?;

    let mut logger = Logger::new(args.log.as_deref())?;
    let started = std::time::Instant::now();
    logger.info(format!(
        "umi-tools-rs dedup started; cwd={}",
        std::env::current_dir()?.display()
    ))?;

    let output_format = output_format(&args);
    let whitelist = load_whitelist(&args)?;
    let temp_dir = make_temp_dir(args.temp_dir.as_deref())?;
    let intermediate_path = if args.no_sort_output {
        None
    } else {
        Some(temp_dir.path().join("dedup.unsorted.bam"))
    };

    let write_target = intermediate_path
        .as_deref()
        .map(WriteTarget::Path)
        .unwrap_or_else(|| match args.output.as_deref() {
            Some(path) if path != Path::new("-") => WriteTarget::Path(path),
            _ => WriteTarget::Stdout,
        });

    let mut reader = open_reader(&args)?;
    let header = Header::from_template(reader.header());
    let target_names = target_names(reader.header());
    let metacontigs = load_metacontig_map(&args, &target_names)?;
    let mut stats = if let Some(prefix) = args.output_stats.as_ref() {
        Some(StatsCollector::new(prefix.clone(), &args, &target_names)?)
    } else {
        None
    };

    let mut rng = match args.random_seed {
        Some(seed) => StdRng::seed_from_u64(seed),
        None => StdRng::from_entropy(),
    };
    let mut builder = BundleBuilder::default();
    let mut mate_keys = HashSet::new();
    let mut n_output = 0_u64;

    {
        let mut writer = open_writer(write_target, &header, output_format.htslib(), &args)?;
        if let Some(metacontigs) = metacontigs.as_ref() {
            process_metacontigs(
                metacontigs,
                &args,
                &target_names,
                &mut builder,
                &mut writer,
                &mut mate_keys,
                &whitelist,
                &mut stats,
                &mut rng,
                &mut n_output,
            )?;
        } else {
            process_stream(
                &mut reader,
                &args,
                &target_names,
                &mut builder,
                &mut writer,
                &mut mate_keys,
                &whitelist,
                &mut stats,
                &mut rng,
                &mut n_output,
            )?;
        }

        let remaining = builder.remaining_keys();
        flush_bundles(
            remaining,
            &mut builder,
            &mut writer,
            &mut mate_keys,
            &whitelist,
            &args,
            &mut stats,
            &mut n_output,
        )?;

        if args.paired && !mate_keys.is_empty() {
            write_mates(&args, &mut writer, &mut mate_keys, &mut n_output)?;
        }
    }

    if let Some(intermediate) = intermediate_path.as_deref() {
        sort_output(
            intermediate,
            args.output.as_deref(),
            output_format,
            &args,
            temp_dir.path(),
        )?;
    }
    if let Some(stats) = stats.as_ref() {
        stats.write(args.method)?;
    }

    logger.info(format!("Number of reads out: {n_output}"))?;
    logger.info(format!(
        "Number of input read1 records considered: {}",
        builder.n_input
    ))?;
    logger.info(format!("Run time: {:.3}s", started.elapsed().as_secs_f64()))?;
    for (event, count) in sorted_events(&builder.events) {
        logger.info(format!("{event}: {count}"))?;
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process_stream(
    reader: &mut bam::Reader,
    args: &DedupArgs,
    target_names: &[String],
    builder: &mut BundleBuilder,
    writer: &mut Writer,
    mate_keys: &mut HashSet<MateKey>,
    whitelist: &Option<HashSet<Vec<u8>>>,
    stats: &mut Option<StatsCollector>,
    rng: &mut StdRng,
    n_output: &mut u64,
) -> Result<()> {
    let mut record = Record::new();
    while let Some(read_result) = reader.read(&mut record) {
        read_result.context("failed to read input alignment")?;
        let read = record.clone();

        if !should_visit_chrom(&read, target_names, args.chrom.as_deref()) {
            continue;
        }

        process_read(
            read,
            None,
            args,
            target_names,
            builder,
            writer,
            mate_keys,
            whitelist,
            stats,
            rng,
            n_output,
        )?;
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process_metacontigs(
    metacontigs: &IndexMap<String, Vec<String>>,
    args: &DedupArgs,
    target_names: &[String],
    builder: &mut BundleBuilder,
    writer: &mut Writer,
    mate_keys: &mut HashSet<MateKey>,
    whitelist: &Option<HashSet<Vec<u8>>>,
    stats: &mut Option<StatsCollector>,
    rng: &mut StdRng,
    n_output: &mut u64,
) -> Result<()> {
    let mut reader = open_indexed_reader(args)?;
    let mut record = Record::new();

    for (gene, transcripts) in metacontigs {
        for transcript in transcripts {
            reader
                .fetch(transcript.as_str())
                .with_context(|| format!("failed to fetch transcript {transcript}"))?;
            while let Some(read_result) = reader.read(&mut record) {
                read_result.context("failed to read transcript alignment")?;
                let mut read = record.clone();
                set_string_tag(&mut read, b"MC", gene)?;
                process_read(
                    read,
                    Some(gene.as_str()),
                    args,
                    target_names,
                    builder,
                    writer,
                    mate_keys,
                    whitelist,
                    stats,
                    rng,
                    n_output,
                )?;
            }
        }

        let remaining = builder.remaining_keys();
        flush_bundles(
            remaining, builder, writer, mate_keys, whitelist, args, stats, n_output,
        )?;
        builder.last_tid = None;
        builder.last_start = 0;
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process_read(
    read: Record,
    meta_gene: Option<&str>,
    args: &DedupArgs,
    target_names: &[String],
    builder: &mut BundleBuilder,
    writer: &mut Writer,
    mate_keys: &mut HashSet<MateKey>,
    whitelist: &Option<HashSet<Vec<u8>>>,
    stats: &mut Option<StatsCollector>,
    rng: &mut StdRng,
    n_output: &mut u64,
) -> Result<()> {
    let flush_keys = builder.add_read(read, meta_gene, args, target_names, rng)?;
    flush_bundles(
        flush_keys, builder, writer, mate_keys, whitelist, args, stats, n_output,
    )
}

impl Default for BundleBuilder {
    fn default() -> Self {
        Self {
            reads: BTreeMap::new(),
            last_tid: None,
            last_start: 0,
            events: HashMap::new(),
            n_input: 0,
        }
    }
}

impl BundleBuilder {
    fn add_read(
        &mut self,
        read: Record,
        meta_gene: Option<&str>,
        args: &DedupArgs,
        target_names: &[String],
        rng: &mut StdRng,
    ) -> Result<Vec<PositionKey>> {
        if read.is_last_in_template() {
            return Ok(Vec::new());
        }

        self.n_input += 1;
        self.event("Input Reads");

        if args.paired {
            if read.is_paired() {
                self.event("Read pairs");
            } else {
                self.event("Unpaired reads");
                if args.unpaired_reads == ReadHandling::Discard {
                    return Ok(Vec::new());
                }
            }
        }

        if read.is_unmapped() {
            self.event(if args.paired {
                if read.is_mate_unmapped() {
                    "Both unmapped"
                } else {
                    "Read 1 unmapped"
                }
            } else {
                "Single end unmapped"
            });
            return Ok(Vec::new());
        }

        if args.paired && read.is_mate_unmapped() {
            self.event("Read 2 unmapped");
            if args.unmapped_reads != ReadHandling::Use {
                return Ok(Vec::new());
            }
        }

        if args.paired && read.is_paired() && read.tid() != read.mtid() {
            self.event("Chimeric read pair");
            if args.chimeric_pairs == ReadHandling::Discard {
                return Ok(Vec::new());
            }
        }

        if let Some(subset) = args.subset {
            if rng.gen_range(0.0..1.0) >= subset {
                self.event("Randomly excluded");
                return Ok(Vec::new());
            }
        }

        if read.mapq() < args.mapping_quality {
            self.event("< MAPQ threshold");
            return Ok(Vec::new());
        }

        let (umi, cell) = match extract_barcodes(&read, args) {
            Ok((umi, cell)) => {
                if args.ignore_umi {
                    (Vec::new(), cell)
                } else {
                    (umi, cell)
                }
            }
            Err(_) if args.ignore_umi && !args.per_cell => (Vec::new(), None),
            Err(_) => {
                self.event("Read skipped, missing umi and/or cell tag");
                return Ok(Vec::new());
            }
        };

        let tid = read.tid();
        let flush_keys =
            self.flush_keys_before_read(&read, args, target_names, meta_gene.is_some())?;
        let (position_key, group_key) = if args.per_gene {
            let gene = gene_for_read(&read, meta_gene, args, target_names)?;
            let gene = match gene {
                Some(gene) => gene,
                None => {
                    self.event("Read skipped, no gene");
                    return Ok(flush_keys);
                }
            };
            (
                PositionKey::Gene(gene.clone()),
                GroupKey {
                    inner: GroupInner::Gene(gene),
                    cell,
                },
            )
        } else {
            let (start, pos, splice_offset) = read_position(&read, args.soft_clip_threshold);
            let spliced_offset = if args.spliced {
                splice_offset.unwrap_or(0)
            } else {
                0
            };
            let tlen = if args.paired && !args.ignore_tlen {
                read.insert_size()
            } else {
                0
            };
            let read_len = if args.read_length { read.seq_len() } else { 0 };

            self.last_start = start;
            (
                PositionKey::Coord { tid, pos },
                GroupKey {
                    inner: GroupInner::Coord {
                        reverse: read.is_reverse(),
                        spliced_offset,
                        tlen,
                        read_len,
                    },
                    cell,
                },
            )
        };

        update_bundle(
            self.reads
                .entry(position_key)
                .or_default()
                .entry(group_key)
                .or_default(),
            umi,
            read,
            args.detection_method,
            rng,
        );

        self.last_tid = Some(tid);
        Ok(flush_keys)
    }

    fn flush_keys_before_read(
        &mut self,
        read: &Record,
        args: &DedupArgs,
        target_names: &[String],
        suppress_flush: bool,
    ) -> Result<Vec<PositionKey>> {
        if suppress_flush {
            return Ok(Vec::new());
        }

        let current_tid = read.tid();
        let Some(last_tid) = self.last_tid else {
            self.last_tid = Some(current_tid);
            return Ok(Vec::new());
        };

        if args.per_gene || args.whole_contig {
            if current_tid != last_tid {
                return Ok(self.remaining_keys());
            }
            return Ok(Vec::new());
        }

        let (start, _, _) = read_position(read, args.soft_clip_threshold);
        let mut flush_all = current_tid != last_tid;
        if start > self.last_start + 1000 {
            flush_all = true;
        }

        if !flush_all {
            return Ok(Vec::new());
        }

        if current_tid != last_tid {
            return Ok(self.remaining_keys());
        }

        let cutoff = start - 1000;
        let keys = self
            .reads
            .keys()
            .filter(|key| match key {
                PositionKey::Coord { tid, pos } => *tid == current_tid && *pos <= cutoff,
                PositionKey::Gene(gene) => target_names
                    .get(current_tid as usize)
                    .is_some_and(|name| name == gene),
            })
            .cloned()
            .collect();
        Ok(keys)
    }

    fn remaining_keys(&self) -> Vec<PositionKey> {
        self.reads.keys().cloned().collect()
    }

    fn event(&mut self, event: impl Into<String>) {
        *self.events.entry(event.into()).or_insert(0) += 1;
    }
}

fn update_bundle(
    bundle: &mut Bundle,
    umi: Vec<u8>,
    read: Record,
    detection_method: Option<DetectionMethod>,
    rng: &mut StdRng,
) {
    let Some(entry) = bundle.get_mut(&umi) else {
        bundle.insert(
            umi,
            BundleEntry {
                read,
                count: 1,
                ties_seen: 0,
            },
        );
        return;
    };

    entry.count += 1;
    match compare_representative(&entry.read, &read, detection_method) {
        Ordering::Greater => {
            entry.read = read;
            entry.ties_seen = 0;
        }
        Ordering::Equal => {
            entry.ties_seen += 1;
            if rng.gen_bool(1.0 / entry.ties_seen as f64) {
                entry.read = read;
            }
        }
        Ordering::Less => {}
    }
}

fn compare_representative(
    current: &Record,
    candidate: &Record,
    detection_method: Option<DetectionMethod>,
) -> Ordering {
    match candidate.mapq().cmp(&current.mapq()) {
        Ordering::Equal => {}
        ordering => return ordering,
    }

    match detection_method {
        Some(DetectionMethod::NH) => compare_lower_aux(current, candidate, b"NH"),
        Some(DetectionMethod::X0) => compare_lower_aux(current, candidate, b"X0"),
        Some(DetectionMethod::XT) => compare_xt(current, candidate),
        None => Ordering::Equal,
    }
}

fn compare_lower_aux(current: &Record, candidate: &Record, tag: &[u8]) -> Ordering {
    let current_value = aux_i64(current, tag);
    let candidate_value = aux_i64(candidate, tag);

    match (current_value, candidate_value) {
        (Some(current), Some(candidate)) => current.cmp(&candidate),
        _ => Ordering::Equal,
    }
}

fn compare_xt(current: &Record, candidate: &Record) -> Ordering {
    let current_unique = aux_string(current, b"XT").as_deref() == Some("U");
    let candidate_unique = aux_string(candidate, b"XT").as_deref() == Some("U");

    match (current_unique, candidate_unique) {
        (true, false) => Ordering::Less,
        (false, true) => Ordering::Greater,
        _ => Ordering::Equal,
    }
}

fn flush_bundles(
    keys: Vec<PositionKey>,
    builder: &mut BundleBuilder,
    writer: &mut Writer,
    mate_keys: &mut HashSet<MateKey>,
    whitelist: &Option<HashSet<Vec<u8>>>,
    args: &DedupArgs,
    stats: &mut Option<StatsCollector>,
    n_output: &mut u64,
) -> Result<()> {
    for key in keys {
        let Some(group_map) = builder.reads.remove(&key) else {
            continue;
        };

        for (_, bundle) in group_map {
            if args.ignore_umi {
                for entry in bundle.values() {
                    write_selected_read(&entry.read, writer, mate_keys, args, n_output)?;
                }
                continue;
            }

            let counts: IndexMap<Vec<u8>, usize> = bundle
                .iter()
                .map(|(umi, entry)| (umi.clone(), entry.count))
                .collect();
            let clusters = cluster_umis(&counts, args.method, args.threshold);
            let selected = selected_umis(&counts, clusters, whitelist);

            if selected.is_empty() {
                continue;
            }

            if let Some(stats) = stats.as_mut() {
                stats.record_bundle(&counts, &selected);
            }

            for selected_umi in selected {
                if let Some(entry) = bundle.get(&selected_umi.parent) {
                    write_selected_read(&entry.read, writer, mate_keys, args, n_output)?;
                }
            }
        }
    }

    Ok(())
}

fn selected_umis(
    counts: &IndexMap<Vec<u8>, usize>,
    clusters: Vec<Vec<Vec<u8>>>,
    whitelist: &Option<HashSet<Vec<u8>>>,
) -> Vec<SelectedUmi> {
    let mut selected = Vec::new();

    for cluster in clusters {
        let Some(parent) = cluster.first() else {
            continue;
        };
        if let Some(whitelist) = whitelist {
            if !whitelist.contains(parent) {
                continue;
            }
        }
        let count = cluster.iter().filter_map(|umi| counts.get(umi)).sum();
        selected.push(SelectedUmi {
            parent: parent.clone(),
            count,
        });
    }

    selected
}

fn write_selected_read(
    read: &Record,
    writer: &mut Writer,
    mate_keys: &mut HashSet<MateKey>,
    args: &DedupArgs,
    n_output: &mut u64,
) -> Result<()> {
    writer
        .write(read)
        .context("failed to write selected read")?;
    *n_output += 1;

    if args.paired {
        if let Some(key) = mate_key_for_read1(read) {
            mate_keys.insert(key);
        }
    }

    Ok(())
}

fn mate_key_for_read1(read: &Record) -> Option<MateKey> {
    if read.is_unmapped() || read.is_mate_unmapped() || read.mtid() < 0 {
        return None;
    }

    Some(MateKey {
        qname: read.qname().to_vec(),
        tid: read.mtid(),
        pos: read.mpos(),
    })
}

fn write_mates(
    args: &DedupArgs,
    writer: &mut Writer,
    mate_keys: &mut HashSet<MateKey>,
    n_output: &mut u64,
) -> Result<()> {
    let mut reader = open_reader(args)?;
    let mut record = Record::new();

    while let Some(read_result) = reader.read(&mut record) {
        read_result.context("failed to read input alignment while recovering mates")?;

        if record.is_unmapped()
            || record.is_mate_unmapped()
            || record.is_first_in_template()
            || record.tid() < 0
        {
            continue;
        }

        let key = MateKey {
            qname: record.qname().to_vec(),
            tid: record.tid(),
            pos: record.pos(),
        };

        if mate_keys.remove(&key) {
            writer.write(&record).context("failed to write mate read")?;
            *n_output += 1;
        }

        if mate_keys.is_empty() {
            break;
        }
    }

    Ok(())
}

fn extract_barcodes(read: &Record, args: &DedupArgs) -> Result<(Vec<u8>, Option<Vec<u8>>)> {
    match args.get_umi_method {
        ExtractUmiMethod::ReadId => extract_barcodes_from_read_id(read, args),
        ExtractUmiMethod::Tag => extract_barcodes_from_tags(read, args),
        ExtractUmiMethod::Umis => extract_barcodes_from_umis_name(read, args),
    }
}

fn extract_barcodes_from_read_id(
    read: &Record,
    args: &DedupArgs,
) -> Result<(Vec<u8>, Option<Vec<u8>>)> {
    let qname = String::from_utf8_lossy(read.qname());
    let parts: Vec<&str> = qname.split(&args.umi_sep).collect();
    let umi = parts
        .last()
        .context("could not extract UMI from read name")?
        .as_bytes()
        .to_vec();
    let cell = if args.per_cell {
        Some(
            parts
                .get(parts.len().saturating_sub(2))
                .context("could not extract cell barcode from read name")?
                .as_bytes()
                .to_vec(),
        )
    } else {
        None
    };

    Ok((umi, cell))
}

fn extract_barcodes_from_tags(
    read: &Record,
    args: &DedupArgs,
) -> Result<(Vec<u8>, Option<Vec<u8>>)> {
    let mut umi = aux_string(read, args.umi_tag.as_bytes())
        .with_context(|| format!("missing UMI tag {}", args.umi_tag))?;
    if let Some(split) = &args.umi_tag_split {
        umi = umi.split(split).next().unwrap_or_default().to_owned();
    }
    if let Some(delim) = &args.umi_tag_delim {
        umi = umi.split(delim).collect::<String>();
    }

    let cell = if args.per_cell {
        let cell_tag = args
            .cell_tag
            .as_deref()
            .context("--per-cell with --extract-umi-method=tag requires --cell-tag")?;
        let mut cell = aux_string(read, cell_tag.as_bytes())
            .with_context(|| format!("missing cell tag {cell_tag}"))?;
        if !args.cell_tag_split.is_empty() {
            cell = cell
                .split(&args.cell_tag_split)
                .next()
                .unwrap_or_default()
                .to_owned();
        }
        if let Some(delim) = &args.cell_tag_delim {
            cell = cell.split(delim).collect::<String>();
        }
        Some(cell.into_bytes())
    } else {
        None
    };

    Ok((umi.into_bytes(), cell))
}

fn extract_barcodes_from_umis_name(
    read: &Record,
    args: &DedupArgs,
) -> Result<(Vec<u8>, Option<Vec<u8>>)> {
    let qname = String::from_utf8_lossy(read.qname());
    let mut umi = None;
    let mut cell = None;

    for element in qname.split(':') {
        if let Some(value) = element.strip_prefix("UMI_") {
            umi = Some(value.as_bytes().to_vec());
        } else if args.per_cell {
            if let Some(value) = element.strip_prefix("CELL_") {
                cell = Some(value.as_bytes().to_vec());
            }
        }
    }

    Ok((
        umi.context("could not extract UMI_ field from read name")?,
        cell,
    ))
}

fn gene_for_read(
    read: &Record,
    meta_gene: Option<&str>,
    args: &DedupArgs,
    target_names: &[String],
) -> Result<Option<String>> {
    if let Some(gene) = meta_gene {
        return Ok(Some(gene.to_owned()));
    }

    if args.per_contig {
        return Ok(target_names.get(read.tid() as usize).cloned());
    }

    let Some(gene_tag) = args.gene_tag.as_deref() else {
        return Ok(None);
    };
    let assigned_tag = args.assigned_tag.as_deref().unwrap_or(gene_tag);
    let assigned = aux_string(read, assigned_tag.as_bytes());
    let gene = aux_string(read, gene_tag.as_bytes());
    let Some(gene) = gene else {
        return Ok(None);
    };
    if gene.is_empty() {
        return Ok(None);
    }

    let skip_re = Regex::new(&args.skip_regex)?;
    if assigned
        .as_deref()
        .is_some_and(|value| skip_re.is_match(value))
    {
        return Ok(None);
    }

    Ok(Some(gene))
}

fn read_position(read: &Record, soft_clip_threshold: f64) -> (i64, i64, Option<i64>) {
    let cigar = read.cigar();
    let cigar_ops: Vec<Cigar> = cigar.iter().copied().collect();
    if cigar_ops.is_empty() {
        return (read.pos(), read.pos(), None);
    }

    if read.is_reverse() {
        let mut pos = cigar.end_pos();
        if matches!(cigar_ops.last(), Some(Cigar::SoftClip(_))) {
            if let Some(Cigar::SoftClip(clip)) = cigar_ops.last() {
                pos += *clip as i64;
            }
        }
        let start = read.pos();
        let is_spliced = if has_refskip(&cigar_ops)
            || matches!(cigar_ops.first(), Some(Cigar::SoftClip(clip)) if *clip as f64 > soft_clip_threshold)
        {
            let reversed: Vec<Cigar> = cigar_ops.iter().rev().copied().collect();
            find_splice(&reversed)
        } else {
            None
        };
        (start, pos, is_spliced)
    } else {
        let mut pos = read.pos();
        if let Some(Cigar::SoftClip(clip)) = cigar_ops.first() {
            pos -= *clip as i64;
        }
        let start = pos;
        let is_spliced = if has_refskip(&cigar_ops)
            || matches!(cigar_ops.last(), Some(Cigar::SoftClip(clip)) if *clip as f64 > soft_clip_threshold)
        {
            find_splice(&cigar_ops)
        } else {
            None
        };
        (start, pos, is_spliced)
    }
}

fn find_splice(cigar: &[Cigar]) -> Option<i64> {
    let mut offset = 0_i64;
    let mut iter = cigar.iter();

    if let Some(Cigar::SoftClip(clip)) = cigar.first() {
        offset = *clip as i64;
        iter.next();
    }

    for op in iter {
        match op {
            Cigar::RefSkip(_) | Cigar::SoftClip(_) => return Some(offset),
            Cigar::Match(bases) | Cigar::Del(bases) | Cigar::Equal(bases) | Cigar::Diff(bases) => {
                offset += *bases as i64;
            }
            Cigar::Ins(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    None
}

fn has_refskip(cigar: &[Cigar]) -> bool {
    cigar.iter().any(|op| matches!(op, Cigar::RefSkip(_)))
}

fn aux_string(read: &Record, tag: &[u8]) -> Option<String> {
    match read.aux(tag).ok()? {
        Aux::String(value) | Aux::HexByteArray(value) => Some(value.to_owned()),
        Aux::Char(value) => Some((value as char).to_string()),
        Aux::I8(value) => Some(value.to_string()),
        Aux::U8(value) => Some(value.to_string()),
        Aux::I16(value) => Some(value.to_string()),
        Aux::U16(value) => Some(value.to_string()),
        Aux::I32(value) => Some(value.to_string()),
        Aux::U32(value) => Some(value.to_string()),
        Aux::Float(value) => Some(value.to_string()),
        Aux::Double(value) => Some(value.to_string()),
        _ => None,
    }
}

fn aux_i64(read: &Record, tag: &[u8]) -> Option<i64> {
    match read.aux(tag).ok()? {
        Aux::I8(value) => Some(value as i64),
        Aux::U8(value) => Some(value as i64),
        Aux::I16(value) => Some(value as i64),
        Aux::U16(value) => Some(value as i64),
        Aux::I32(value) => Some(value as i64),
        Aux::U32(value) => Some(value as i64),
        Aux::String(value) => value.parse().ok(),
        _ => None,
    }
}

fn validate_args(args: &DedupArgs) -> Result<()> {
    if args.output_stats.is_some() && args.ignore_umi {
        bail!("--output-stats and --ignore-umi cannot be used together");
    }
    if args.gene_tag.is_some() && args.per_contig {
        bail!("use either --per-contig or --gene-tag, not both");
    }
    if args.gene_tag.is_some() && !args.per_gene {
        bail!("--gene-tag requires --per-gene");
    }
    if args.per_contig && !args.per_gene {
        bail!("--per-contig requires --per-gene");
    }
    if args.per_gene && !args.per_contig && args.gene_tag.is_none() {
        bail!("--per-gene requires --gene-tag or --per-contig");
    }
    if args.gene_transcript_map.is_some() && !args.per_contig {
        bail!("--gene-transcript-map requires --per-contig and --per-gene");
    }
    if args.filter_umi && args.umi_whitelist.is_none() {
        bail!("--filter-umi requires --umi-whitelist");
    }
    if let Some(subset) = args.subset {
        if !(0.0..=1.0).contains(&subset) {
            bail!("--subset must be between 0 and 1");
        }
    }
    if let Some(level) = args.compresslevel {
        if level > 9 {
            bail!("--compresslevel must be between 0 and 9");
        }
    }
    Ok(())
}

fn load_whitelist(args: &DedupArgs) -> Result<Option<HashSet<Vec<u8>>>> {
    if !args.filter_umi {
        return Ok(None);
    }

    let first = read_whitelist_file(
        args.umi_whitelist
            .as_deref()
            .context("--filter-umi requires --umi-whitelist")?,
    )?;

    let whitelist = if let Some(second_path) = args.umi_whitelist_paired.as_deref() {
        let second = read_whitelist_file(second_path)?;
        let mut paired = HashSet::new();
        for left in &first {
            for right in &second {
                let mut umi = left.clone();
                umi.extend_from_slice(right);
                paired.insert(umi);
            }
        }
        paired
    } else {
        first
    };

    Ok(Some(whitelist))
}

fn read_whitelist_file(path: &Path) -> Result<HashSet<Vec<u8>>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read whitelist {}", path.display()))?;
    let mut values = HashSet::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        if let Some(first) = line.split('\t').next() {
            values.insert(first.as_bytes().to_vec());
        }
    }
    Ok(values)
}

fn load_metacontig_map(
    args: &DedupArgs,
    target_names: &[String],
) -> Result<Option<IndexMap<String, Vec<String>>>> {
    let Some(map_path) = args.gene_transcript_map.as_deref() else {
        return Ok(None);
    };
    if args.chrom.is_some() {
        return Ok(None);
    }

    let target_set: HashSet<&str> = target_names.iter().map(String::as_str).collect();
    let text = std::fs::read_to_string(map_path)
        .with_context(|| format!("failed to read {}", map_path.display()))?;
    let mut metacontigs: IndexMap<String, Vec<String>> = IndexMap::new();

    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        if line.trim().is_empty() {
            break;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 2 {
            bail!(
                "gene transcript map line must contain two tab-separated fields: {}",
                line
            );
        }
        let gene = fields[0].to_owned();
        let transcript = fields[1].to_owned();
        if target_set.contains(transcript.as_str()) {
            let transcripts = metacontigs.entry(gene).or_default();
            if !transcripts.contains(&transcript) {
                transcripts.push(transcript);
            }
        }
    }

    Ok(Some(metacontigs))
}

#[derive(Clone, Copy, Debug, Default)]
struct UmiAggregate {
    median_counts: usize,
    times_observed: usize,
    total_counts: usize,
}

fn count_instances_by_count(rows: &[(Vec<u8>, usize)]) -> BTreeMap<usize, usize> {
    let mut counts = BTreeMap::new();
    for (_, count) in rows {
        *counts.entry(*count).or_insert(0) += 1;
    }
    counts
}

fn aggregate_umi_rows(rows: &[(Vec<u8>, usize)]) -> BTreeMap<Vec<u8>, UmiAggregate> {
    let mut grouped: BTreeMap<Vec<u8>, Vec<usize>> = BTreeMap::new();
    for (umi, count) in rows {
        grouped.entry(umi.clone()).or_default().push(*count);
    }

    grouped
        .into_iter()
        .map(|(umi, mut counts)| {
            counts.sort_unstable();
            let total_counts = counts.iter().sum();
            let times_observed = counts.len();
            let median_counts = median_floor(&counts);
            (
                umi,
                UmiAggregate {
                    median_counts,
                    times_observed,
                    total_counts,
                },
            )
        })
        .collect()
}

fn median_floor(values: &[usize]) -> usize {
    if values.is_empty() {
        return 0;
    }
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[mid - 1] + values[mid]) / 2
    } else {
        values[mid]
    }
}

fn average_umi_distance(umis: &[Vec<u8>]) -> f64 {
    if umis.len() <= 1 {
        return -1.0;
    }

    let mut total = 0_usize;
    let mut pairs = 0_usize;
    for i in 0..umis.len() {
        for j in (i + 1)..umis.len() {
            total += hamming_distance_lenient(&umis[i], &umis[j]);
            pairs += 1;
        }
    }

    total as f64 / pairs as f64
}

fn hamming_distance_lenient(a: &[u8], b: &[u8]) -> usize {
    let mismatches = a
        .iter()
        .zip(b.iter())
        .filter(|(left, right)| left != right)
        .count();
    mismatches + a.len().abs_diff(b.len())
}

fn tally_digitized(values: &[f64], bins: &[i32]) -> Vec<usize> {
    let mut tally = vec![0_usize; bins.len()];
    for value in values {
        let idx = bins
            .iter()
            .position(|bin| *value <= *bin as f64)
            .unwrap_or_else(|| bins.len().saturating_sub(1));
        if let Some(slot) = tally.get_mut(idx) {
            *slot += 1;
        }
    }
    tally
}

fn suffixed_path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut value = prefix.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

fn cluster_method_name(method: ClusterMethod) -> &'static str {
    match method {
        ClusterMethod::Unique => "unique",
        ClusterMethod::Percentile => "percentile",
        ClusterMethod::Cluster => "cluster",
        ClusterMethod::Adjacency => "adjacency",
        ClusterMethod::Directional => "directional",
    }
}

enum WriteTarget<'a> {
    Path(&'a Path),
    Stdout,
}

fn open_reader(args: &DedupArgs) -> Result<bam::Reader> {
    let mut reader = bam::Reader::from_path(&args.input)
        .with_context(|| format!("failed to open input {}", args.input.display()))?;
    if args.threads > 0 {
        reader.set_threads(args.threads)?;
    }
    if let Some(reference) = args.reference_filename.as_deref() {
        reader.set_reference(reference)?;
    }
    Ok(reader)
}

fn open_indexed_reader(args: &DedupArgs) -> Result<bam::IndexedReader> {
    let mut reader = bam::IndexedReader::from_path(&args.input)
        .with_context(|| format!("failed to open indexed input {}", args.input.display()))?;
    if args.threads > 0 {
        reader.set_threads(args.threads)?;
    }
    if let Some(reference) = args.reference_filename.as_deref() {
        reader.set_reference(reference)?;
    }
    Ok(reader)
}

fn open_writer(
    target: WriteTarget<'_>,
    header: &Header,
    format: Format,
    args: &DedupArgs,
) -> Result<Writer> {
    let mut writer = match target {
        WriteTarget::Path(path) => Writer::from_path(path, header, format)
            .with_context(|| format!("failed to open output {}", path.display()))?,
        WriteTarget::Stdout => {
            Writer::from_stdout(header, format).context("failed to open stdout")?
        }
    };
    if args.threads > 0 {
        writer.set_threads(args.threads)?;
    }
    if let Some(reference) = args.reference_filename.as_deref() {
        writer.set_reference(reference)?;
    }
    if let Some(level) = args.compresslevel {
        writer.set_compression_level(CompressionLevel::Level(level))?;
    }
    Ok(writer)
}

fn set_string_tag(read: &mut Record, tag: &[u8], value: &str) -> Result<()> {
    if read.aux(tag).is_ok() {
        read.remove_aux(tag).with_context(|| {
            format!(
                "failed to remove existing {} tag",
                String::from_utf8_lossy(tag)
            )
        })?;
    }
    read.push_aux(tag, Aux::String(value))
        .with_context(|| format!("failed to set {} tag", String::from_utf8_lossy(tag)))
}

fn output_format(args: &DedupArgs) -> AlignmentFormat {
    if args.out_sam {
        return AlignmentFormat::Sam;
    }
    if let Some(format) = args.out_format {
        return format;
    }
    if let Some(output) = args.output.as_deref() {
        match output.extension().and_then(|ext| ext.to_str()) {
            Some("sam") => return AlignmentFormat::Sam,
            Some("cram") => return AlignmentFormat::Cram,
            _ => {}
        }
    }
    AlignmentFormat::Bam
}

fn make_temp_dir(temp_dir: Option<&Path>) -> Result<TempDir> {
    let mut builder = tempfile::Builder::new();
    builder.prefix("umi-tools-rs.");
    match temp_dir {
        Some(path) => builder.tempdir_in(path),
        None => builder.tempdir(),
    }
    .context("failed to create temporary directory")
}

fn sort_output(
    intermediate: &Path,
    output: Option<&Path>,
    output_format: AlignmentFormat,
    args: &DedupArgs,
    temp_dir: &Path,
) -> Result<()> {
    let mut command = Command::new("samtools");
    command
        .arg("sort")
        .arg("--no-PG")
        .arg("-O")
        .arg(output_format.samtools_name())
        .arg("-T")
        .arg(temp_dir.join("sort"));

    if args.threads > 0 {
        command.arg("-@").arg(args.threads.to_string());
    }
    if let Some(reference) = args.reference_filename.as_deref() {
        command.arg("--reference").arg(reference);
    }
    if let Some(output) = output {
        if output != Path::new("-") {
            command.arg("-o").arg(output);
        }
    }
    command.arg(intermediate);

    let status = command.status().context("failed to run samtools sort")?;
    if !status.success() {
        bail!("samtools sort failed with status {status}");
    }

    Ok(())
}

fn should_visit_chrom(read: &Record, target_names: &[String], chrom: Option<&str>) -> bool {
    let Some(chrom) = chrom else {
        return true;
    };
    if read.tid() < 0 {
        return false;
    }
    target_names
        .get(read.tid() as usize)
        .is_some_and(|name| name == chrom)
}

fn target_names(header: &bam::HeaderView) -> Vec<String> {
    header
        .target_names()
        .iter()
        .map(|name| String::from_utf8_lossy(name).into_owned())
        .collect()
}

fn sorted_events(events: &HashMap<String, u64>) -> Vec<(&String, &u64)> {
    let mut values: Vec<_> = events.iter().collect();
    values.sort_by(|a, b| a.0.cmp(b.0));
    values
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record_with_xt(value: &str) -> Record {
        let mut record = Record::new();
        record.push_aux(b"XT", Aux::String(value)).unwrap();
        record
    }

    #[test]
    fn xt_tie_break_prefers_unique_alignment() {
        let unique = record_with_xt("U");
        let repeat = record_with_xt("R");

        assert_eq!(compare_xt(&repeat, &unique), Ordering::Greater);
        assert_eq!(compare_xt(&unique, &repeat), Ordering::Less);
        assert_eq!(compare_xt(&unique, &unique), Ordering::Equal);
    }
}
