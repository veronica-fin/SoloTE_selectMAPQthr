# SoloTE
[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs42003--022--04020--5-yellow)](https://www.nature.com/articles/s42003-022-04020-5)

SoloTE quantifies transposable element (TE) expression in single-cell RNA-seq data alongside gene expression. This fork reports pipeline version **1.10** and adds a configurable MAPQ threshold, fractional TE overlap, strand-aware annotation, paired-read reconciliation, and combined gene–TE features.

## 0. DEPENDENCIES
SoloTE requires the following tools to be installed, and available in your PATH environment variable:
- Samtools v1.16 or higher (http://www.htslib.org/download/)
- BEDTools v2.29.2 or higher (https://github.com/arq5x/bedtools2/releases)
- R v4 or higher (https://www.r-project.org/)

Python 3.9.5 or higher should also be available, with the following modules for the pipeline:
- Pysam (https://pysam.readthedocs.io/en/latest/installation.html)
- Pandas v1.5.0 or higher (https://pandas.pydata.org/)
- NumPy (https://numpy.org/)
- intervaltree (https://pypi.org/project/intervaltree/)

The RepeatMasker helper additionally requires `requests`, `colorama`, and `tqdm`.

Both `python` and `Rscript` must be available in `PATH`: the pipeline invokes its Python annotation helper and R matrix writer using those commands. `generate_mtx.R` uses base R and requires no additional R packages. Keep the pipeline and its helper scripts together in the repository directory. The pipeline uses the multiprocessing `fork` start method, so run it on a system that supports it, such as Linux or macOS.


## 1. SETTING UP NECESSARY FILES

RepeatMasker files (\*rm.out) for several genomes, can be found at the UCSC webpage. The helper utility, `SoloTE_RepeatMasker_to_BED.py`, is packed with SoloTE, and it streamlines the download of RepeatMasker file from UCSC, and conversion to BED format. It can be run like this:
```
python SoloTE_RepeatMasker_to_BED.py -g GenomeVersion
```
where
- GenomeVersion: genome version identifier. For example, `hg38` for human, and `mm10` for mouse.

For example, `-g hg38` produces `hg38_rmsk.bed`. To convert a local gzip-compressed RepeatMasker `.out` file instead:

```bash
python SoloTE_RepeatMasker_to_BED.py -i hg38.fa.out.gz
```

Additionally, the utility can be called with the `-l` option, and it will list identifiers that can be supplied as the `-g` parameter, of all available genomes at UCSC (first 5 lines shown next):
```
python SoloTE_RepeatMasker_to_BED.py -l
ailMel1	| Panda [Ailuropoda melanoleuca, Dec. 2009 (BGI-Shenzhen 1.0/ailMel1)]
allMis1	| American alligator [Alligator mississippiensis, Aug. 2012 (allMis0.2/allMis1)]
anoCar1	| Lizard [Anolis carolinensis, Feb. 2007 (Broad/anoCar1)]
anoCar2	| Lizard [Anolis carolinensis, May 2010 (Broad AnoCar2.0/anoCar2)]
anoGam1	| A. gambiae [Anopheles gambiae, Feb. 2003 (IAGEC MOZ2/anoGam1)]
```

If you have your own RepeatMasker file and/or a file corresponding to Transposable Elements obtained from another tool, make sure to adapt it to the following format for SoloTE:
```
sequenceName	startPosition	endPosition	sequenceName|startPosition|endPosition|TE_Subfamily:TE_Family:TE_Class|score|strand	score	strand
```

Use a tab-separated, headerless BED file with six columns. Column 4 is the TE identifier, combining its locus, Subfamily:Family:Class annotation, score, and strand. Column 6 must contain `+` or `-` for stranded analysis. Chromosome names must match the BAM file. The file should look like this:
```
chr1	11505	11675	chr1|11505|11675|L1MC5a:L1:LINE|25.1|-	25.1	-
chr1	11678	11780	chr1|11678|11780|MER5B:hAT-Charlie:DNA|29.4|-	29.4	-
chr1	15265	15355	chr1|15265|15355|MIR3:MIR:SINE|23.0|-	23.0	-
chr1	18907	19048	chr1|18907|19048|L2a:L2:LINE|33.8|+	33.8	+
chr1	19972	20405	chr1|19972|20405|L3:CR1:LINE|31.2|+	31.2	+
```


## 2. RUNNING SOLOTE

The input must be a coordinate-sorted BAM file. Reads used for TE quantification must have cell-barcode (`CB`) and UMI (`UB`) tags with values other than `-`. Existing gene annotations are read from the `GN` tag; reads with missing `GN` or `GN=-` are treated as non-genic. The pipeline indexes the input BAM, so its directory must be writable.

Run from the repository directory, or use the absolute path to `SoloTE_pipeline.py`:

```bash
python SoloTE_pipeline.py \
  --threads 8 \
  --bam /path/to/coordinate_sorted.bam \
  --teannotation /path/to/hg38_rmsk.bed \
  --outputprefix sample \
  --outputdir /path/to/results
```

Use paths and an output prefix without spaces or shell metacharacters, because the current pipeline constructs shell commands from these values.

### Required arguments

| Argument | Description |
| --- | --- |
| `-t`, `--threads` | Number of threads/processes to use. |
| `-b`, `--bam` | Input BAM file with `CB` and `UB` tags and, for gene counts, `GN` tags. |
| `-a`, `--teannotation` | TE annotation in the six-column BED format above. |
| `-o`, `--outputprefix` | Prefix for output files and directories. |
| `-d`, `--outputdir` | Results directory; created if it does not exist. |

### Optional arguments

| Argument | Default | Description |
| --- | --- | --- |
| `--dual` | Off | Include gene-annotated reads when assigning TE expression. Without this flag, only reads with missing `GN` or `GN=-` are selected for TE annotation. |
| `-p`, `--overlap_pct` | `0.5` | Minimum fraction of the aligned read length that must overlap a TE. Supply a value from `0.0` to `1.0`; the current parser does not validate this range. |
| `--locusMAPQthr` | `255` | Minimum MAPQ for a locus-specific TE assignment. Reads below the threshold receive a Subfamily:Family:Class assignment instead; they are not discarded by this threshold. |
| `-s`, `--strand_mode` | `Unstranded` | One of `Unstranded`, `Forward`, or `Reverse` (case-sensitive). |

For example, to include gene-annotated reads, require 80% TE overlap, use MAPQ 2 for locus assignment, and apply reverse-strand logic:

```bash
python SoloTE_pipeline.py \
  --threads 8 \
  --bam /path/to/coordinate_sorted.bam \
  --teannotation /path/to/hg38_rmsk.bed \
  --outputprefix sample \
  --outputdir /path/to/results \
  --dual \
  --overlap_pct 0.8 \
  --locusMAPQthr 2 \
  --strand_mode Reverse
```

Choose the MAPQ threshold and strand mode to match the aligner and library preparation used for your BAM file.

## 3. TE ASSIGNMENT AND COUNTING

The pipeline calls `annotateBAM_strandSpecific_paired2.py`. The other `annotateBAM*.py` files are not used by the current pipeline.

- **Overlap:** For each read, the annotator sums TE overlap across aligned blocks, excluding skipped reference regions such as introns. It selects the TE with the greatest overlap and requires that overlap to be at least `query_alignment_length × overlap_pct`. Equal-overlap ties have no explicit tie-breaking rule.
- **MAPQ:** Each read is assigned the full BED locus identifier when its MAPQ is at least `--locusMAPQthr`; otherwise it is assigned the TE's Subfamily:Family:Class identifier. TE tags receive the prefix `SoloTE|`.
- **Strand:** `Unstranded` ignores TE strand. In `Forward` mode, the inferred RNA strand matches the alignment strand of read 1 (or a single-end read), and is opposite to read 2. `Reverse` inverts that relationship. Only TEs matching the inferred RNA strand are considered.
- **Pairs:** Reads marked as paired are buffered by query name. If both assignments agree, that assignment is used for both records. If only one has a TE assignment, it is applied to both. If both have different nonempty assignments, both receive `SoloTE|AmbiguousTE`. This comparison uses the MAPQ-dependent assignment strings, so a locus assignment and a subfamily assignment can also conflict.
- **Gene–TE features:** With `--dual`, an existing gene name is preserved and the TE assignment is appended, for example `GeneA;SoloTE|chr1|...`. These combined identifiers remain single features in all five matrices; they are not split between gene and TE rows or included in the pure-TE family/class/subfamily aggregations.
- **UMIs:** Counts are the number of distinct `UB` values per feature and `CB` within each chromosome, followed by summation across chromosomes. UMI deduplication is not performed globally across chromosomes or across different features. Original gene-annotated BAM records are merged with the annotated TE records, so gene and combined gene–TE rows can represent the same underlying molecule in dual mode.

### Current paired-read limitations

Pair reconciliation only sees records retained by the initial TE-region and tag filters. The selection command uses `samtools view -M -L`; it does not explicitly fetch mates outside the TE regions. The annotator does not flush unmatched reads from its pair buffer at the end, so paired records without a retained matching record are omitted from the TE-annotated BAM. It matches records by query name without additional checks for secondary or supplementary alignments. Conflicting assignments produce `SoloTE|AmbiguousTE`, for which downstream aggregation has no dedicated handling and may fail if the expected TE annotation columns are absent.

## 4. OUTPUTS

For `--outputprefix sample --outputdir /path/to/results`, successfully generated matrices are written to:

```text
/path/to/results/sample_SoloTE_output/
├── sample_legacytes_MATRIX/
├── sample_locustes_MATRIX/
├── sample_classtes_MATRIX/
├── sample_familytes_MATRIX/
└── sample_subfamilytes_MATRIX/
```

Each matrix directory contains uncompressed `matrix.mtx`, `barcodes.tsv`, and `features.tsv` files. Matrix rows are features and columns are cell barcodes. `features.tsv` has two identical columns containing the feature identifier/name.

All matrices include gene features and combined gene–TE features, when present. Their pure-TE content differs:

| Matrix | Pure-TE features |
| --- | --- |
| `legacytes` | Original TE assignments: a mixture of full locus identifiers and below-threshold Subfamily:Family:Class identifiers. |
| `locustes` | Locus-specific TE assignments only. The current implementation recognizes these by the pattern `SoloTE\|chr`, so locus features on chromosomes without a `chr` prefix are excluded from this matrix. |
| `classtes` | TE counts summed by class. |
| `familytes` | TE counts summed by family; empty/missing family annotations are excluded. |
| `subfamilytes` | TE counts summed by subfamily. |

Temporary and intermediate files remain in `/path/to/results/sample_SoloTE_temp/`, including:

- `sample_final.bam` and its index: merged gene and TE records, retaining `GN`, `CB`, and `UB` tags.
- `sample_allcounts.txt`: headerless per-chromosome counts with feature, barcode, and UMI count columns.
- `sample.chimeric_counts.tsv`: combined gene–TE counts with the header `feature`, `barcode`, `umi`. Here, “chimeric” means a feature identifier containing both a semicolon and `SoloTE`; it is not a separate structural-chimera detection result.
- `sample_legacytes.txt`, `sample_locustes.txt`, `sample_classtes.txt`, `sample_familytes.txt`, and `sample_subfamilytes.txt`: matrix input tables with the header `Gene_id`, `barcode`, `umi`.

`sample_SoloTE.stats` is written to the directory from which the pipeline was launched, rather than to `--outputdir`. Its totals include pure genes and pure TEs but exclude combined gene–TE features, even though those features are present in the matrices. Its locus-specific statistics use the same `SoloTE|chr` pattern as the locus matrix.

Existing intermediates are reused when their files are found. When changing the input, annotation, MAPQ threshold, overlap fraction, strand mode, or dual mode, use a new output prefix or a clean results directory to avoid reusing counts from a previous run. If an R matrix-generation step fails, the pipeline logs an error and continues with the remaining matrices; check that all five output directories were produced.

See [Seurat_analysis_with_SoloTE.md](Seurat_analysis_with_SoloTE.md) for a downstream analysis example. For this pipeline, select one of the matrix subdirectories above when loading the data; the bundled example uses compressed files, while newly generated matrices are uncompressed.
