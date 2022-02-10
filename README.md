<img src="https://github.com/chasewnelson/VCFgenie/blob/main/VCFgenie_logo.png?raw=true" title="VCFgenie logo by Mitch Lin" alt="VCFgenie logo by Mitch Lin" align="middle">

# VCFgenie
**VCFgenie** is a Python program for dynamic (site-by-site) processing of within-sample (pooled-sequencing) variants in a VCF file. The program implements arbitrarily complex filtering rules using any keys in the `INFO` or `sample` columns, and controls for a user-provided false-positive count.

The program is executed at the Unix command line or Mac Terminal as follows:

	VCFgenie.py --seq_len=7857 --error_per_site=0.01103 --num_samples=1 --FP_cutoff=0.05 --VCF_files example.vcf > example.out

Find specific [examples](#examples) below. 

## <a name="contents"></a>Contents

* [Description](#description)
* [How it Works](#how-it-works)
* [False-Positive Count Method](#false-positive-count-method)
* [Allele Processing Rules](#allele-processing-rules)
* [Options](#options)
* [Examples](#examples)
* [Output](#output)
* [TODO](#todo)
* [Troubleshooting](#troubleshooting)
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)
* [Contact](#contact)

## <a name="description"></a>Description

<a target="_blank" rel="noopener noreferrer" href="http://samtools.github.io/hts-specs/VCFv4.2.pdf">**Variant Call Format (VCF)** files</a> are widely used to store genetic polymorphism data. These data may summarize either **(1) individual-based sequencing data**, wherein one sequencing reaction corresponds to one individual (any ploidy); or **(2) pooled-sequencing data**, wherein one sequencing reaction corresponds to genetic material from multiple individuals (any ploidy) representing a population under study. Common examples of pooled-sequencing data include somatic tissues containing populations of tumor cells; samples combining known quantities of genetic material from multiple indivudals; and within-host populations of pathogens such as viruses.

VCF files are the result of single nucleotide polymorphism (SNP) 'calling', for which numerous tools exist, including Freebayes, Genome Analysis Tool Kit (GATK), Ion Proton Variant Caller (TVC), Intrahost Variant Analysis of Replicates (iVar), LowFreq, Samtools mpileup, and VarScan 2.

After SNP calling, various tools and libraries also exist for downstream analysis and filtering of VCF files. However, current options are designed for filtering but not modifying VCF records using a predefined set of relevant parameters; treat reference (REF) and alternative (ALT) alleles inconsistently; or do not provide options for handling sites that are multiallelic (i.e., more than one ALT allele). Thus, there is a need for a tool that implements arbitrary user-provided rules to both filter and modify records, handle failing REF and multiple ALT alleles consistently, and control for a false-positive count among multiple samples, specifically in a way that is reproducible (i.e., no custom code needed). VCFgenie was written to meet this need, and allows the aforementioned steps to be performed before downstream (e.g., evolutionary and population genetic) analyses that rely on high-confidence polymorphism data.

## <a name="how-it-works"></a>How it Works
VCFgenie is written in Python 3.8. Its dependencies include the following libraries, which the user may need to install (e.g., using `pip` or `conda`):

* `argparse`
* `collections`
* `numpy`
* `operator`
* `os`
* `re`
* `scipy`
* `sys`
* `typing`

The program examines one or more user-provided VCF file(s) to filter, process, and modify records based on an arbitrarily complex set of user-provided rules. One record corresponds to one line of the VCF file, which reports a miniumum of two alleles: REF and ALT. Arbitrarily large numbers of ALT alleles (i.e., MULTIALLELIC sites) are allowed. 

At present, VCFgenie explicitly examines the first sample reported in the file, i.e., the sample reported in column 10, with identifying keys provided in column 9 (FORMAT). Because VCF file conventions differ widely, VCFgenie allows the user to specify which FORMAT data key should be used to determine allele frequencies:

1. `allele count`: the number of reads matching the ALT allele (**DEFAULT**: `AC`)
2. `allele frequency`: the frequency of the ALT allele (**DEFAULT**: `AF`)
3. `read depth (coverage)`: the number of mapped reads overlapping the site (**DEFAULT**: `DP`)

Thus, to determine allele frequencies and coverage, VCFgenie currently analyzes `AC`, `AF`, and `DP` (coverage) based on the `sample` data (column 10), not the `INFO` data (column 8). Many VCF files contain only one pooled sample, in which case this will be sufficient. However, if the user wishes to analyze multiple sample columns present in the same file, it is currently necessary to create a separate VCF file for each sample to be analyzed, with that sample present in column 10. Future instantiations of VCFgenie will allow (1) analysis based on `AC`/`AF`/`DP` data from the INFO column; and (2) automatic analysis of multiple samples, i.e., columns 10 and above.

Additional user-provided filtering criteria may be specified using any key/value pairs present in the `INFO` or `FORMAT`/`sample` columns, described in [Options](#options). Both REF and ALT alleles can fail if they do not meet a criterion, and can fail multiple criteria. Summary statistics reported in the output may therefore include redundant data, e.g., the same allele may be classified as `failAF` and `failDP` (see [Allele Processing Rules](allele-processing-rules)).

## <a name="false-positive-count-method"></a>False-Positive Count Method

In addition to simple filtering rules, VCFgenie implements a dynamic (site-by-site) binomial false-positive (FP) count filter based on the <a target="_blank" rel="noopener noreferrer" href="https://github.com/chasewnelson/SARS-CoV-2-ORF3d#figure-8">filter_VCF.py</a> script described by <a target="_blank" rel="noopener noreferrer" href="https://elifesciences.org/articles/59633">Nelson et al. (2020)</a>.

The method is **dynamic** because each site is evaluated separately, with results depending on that site's `DP` and each allele's `AC`. The method is **binomial** because, given a sequencing error rate per base (`p`) and coverage (`DP`), it determines the probability that each allele's `AC` (`x`) could have resulted from sequencing error alone. Finally, the method controls for a **false-positive count** defined as a user-provided maximum acceptable number of false variant calls, rejecting any variants that are too likely the result of error alone. This false-positive count is a function of four parameters, each of which must or can be supplied by the user (see [Options](#options)):

1. `--error_per_site`: error rate per site expected as the result of sample preparation, library preparation, and sequencing (currently assumes all nucleotide errors equally probable; see [TODO](#todo))
2. `--DP_key`: the read depth (coverage) at a site (i.e., each **read** is another opportunity for error)
3. `--seq_len`: length of reference sequence (e.g., contig, chromosome, or genome) in nucleotides (i.e., each **site** is another opportunity for error)
4. `--num_samples`: number of samples examined (i.e., each **sample** is another opportunity for error).

The binomial calculation uses the `scipy.stats` function `binom.cdf`:

	binom.cdf(x, n, p)

where:

	n = read depth (coverage)
	p = error rate / 3
	x = allele count (number of reads)
	FP cutoff = (1 - binom.cdf(x - 1, n, p)) * seq_len * num_samples
	
Thus, the method assumes a single nucleotide variant (SNV) model in which all three possible nucleotide errors occur at the same rate. Future instantiations of VCFgenie will allow more complex error models to be specified.

To implement the above method, the user must supply a **FP cutoff** (`--FP_cutoff`). For example, one might wish to use a FP count cutoff of 0.05 for applications requiring conservative SNP calling. For other contexts, an FP count cutoff of 1 (a maximum acceptable mean of 1 false variant across the entire analysis) or more may be acceptable.

Note that the advantage of using the FP method is that it precludes the need for a fixed, arbitrary `AF` cutoff (although one may still be applied using `--min_AF`).

## <a name="allele-processing-rules"></a>Allele Processing Rules

Resolution of failing alleles requires careful consideration. VCFgenie always maintains the same `DP` (coverage) at a site before and after processing and filtering. However, allele counts (`AC`) may be reassigned to different alleles, depending on which alleles pass the filtering criteria.

The following filtering rules are tested:

1. `failZeroAC`: allele fails when 0 reads support it (REF or ALT)
2. `failDP`: read depth (coverage) fails `--min_DP`
3. `failAC`: allele count (REF or ALT) fails `--min_AC`
4. `failMinAF`: allele frequency (REF or ALT) fails `--min_AF`
5. `failMaxAF`: allele frequency (REF or ALT) fails `--max_AF`
6. `failINFO`: fails one or more of the user-provided `--INFO_rules`
7. `failsample`: fails one or more of the user-provided `--sample_rules`
8. `failFP`: fails `--FP_cutoff`

When appropriate, the following data will be added to the `INFO` column:

1. `DECISION`=(`fixedREF`|`fixedALT`|`fail`|`failZeroAC`|`failDP`|`failAC`|`failAF`|`failINFO`|`failsample`|`failFP`|`pass`): arbitrarily many flags may be added, because each allele (including REF) may fail for multiple reasons
2. `STATUS`=(`PASS`|`FAIL`): one flag per ALT allele will be added
3. `fixedREF` **[Flag]**: site fixed for the REF allele (100% reference)
4. `fixedALT` **[Flag]**: site fixed for a particular ALT allele (100% non-reference)
5. `MULTIALLELIC` **[Flag]**: the site is multiallelic (more than one ALT allele)
6. `REF_FAIL` **[Flag]**: the REF allele failed to meet the criteria

For a site with one or more failing alleles, **one flag per allele** (including REF) will be reported in the `FILTER` column (semicolon-separated). This will include `PASS` for any passing alleles. Distinct values will be reported only once, in alphabetical order (i.e., order is not meaningful). Thus, a site that has one ALT that fails will have a minimum of two flags: one for the REF and one or more for the ALT (e.g., `PASS;failAC`).

The reads corresponding to a failing allele (REF or ALT) are assumed to represent sequencing error. Thus, when an allele fails, its reads are given to the allele(s) that passes.

**In most contexts,** the majority of sites are **biallelic** (one REF and one ALT allele) and the REF is the major allele. In such cases, if the ALT allele fails, the ALT reads will be assumed to represent error, and will be re-assigned to the REF (i.e., the REF allele will be fixed at 100%). Complementarily, at a biallelic site where the ALT is the major allele and the REF fails, the REF reads will be reassigned to the ALT (i.e., the ALT allele will be fixed at 100%).

**Other rare instances** may occur, especially with pooled-sequencing data involving poorly-supported REF alleles and multiallelic sites. Methods of application and resolution are listed below:

1. If an entire site fails (e.g., due to insufficient coverage), or if all alleles (including REF) fail at that site, all reads are given to the **major** (most common) allele, whether REF or ALT.
2. After all alleles are evaluated, reads from failing allele(s) are distributed to the passing allele(s) in proportion to each passing allele's frequency among the read total attributed to all passing alleles.
3. If a key supplied in `--INFO_rules` or `--sample_rules` has `Number=A`, then there is one value per ALT allele, and the rule is evaluated separately for each ALT (REF cannot be evaluated).
4. If a key supplied in `--INFO_rules` or `--sample_rules` has `Number=1`, then there is one value for the site, and the rule is applied equally to all alleles (REF and ALT).
5. Any flags other than `PASS` that were present in the `FILTER` column of the VCF input will be mainained and carried over from the original file.
6. `NOCALL`: `PASS` is converted to `NOCALL` at sites where a major allele fails. In other words, this flag is added to any site at which a MAJOR ALLELE fails, while any `PASS` flags from previous passing (minor) alleles are removed.
7. It is possible for a minor allele to pass while a major allele fails, e.g., when a `--sample_rule` specifies a criterion using a key with `Number=A` (only ALT alleles have a value) and an ALT allele is the major allele. In such cases, **no** `AF` **corrections are performed**; `AC` and `AF` are left as they were prior to evaluation, and passing alleles are converted to `NOCALL`.

## <a name="options"></a>Options
Call VCFgenie using the following options: 

**REQUIRED:**

* `-i`, `--VCF_files` **[FILE(S)]**: input variant call format (VCF) file(s)
* `-e`, `--error_per_site` **[float]**: error rate per site expected as the result of sample preparation, library preparation, and sequencing (currently assumes all nucleotide errors equally probable
* `-L`, `--seq_len` **[int]**: length of reference sequence (e.g., contig, chromosome, or genome) in nucleotides
* `-n`, `--num_samples` **[int]**: number of samples (VCF files) in full analysis
* `-f`, `--FP_cutoff` **[float]**: analysis-wide false-positive (FP) count cutoff

**OPTIONAL:**

* `-o`, `--out_dir` **[DIR]**: output directory name (**DEFAULT**: `VCFgenie_output`)
* `-c`, `--AC_key` **[str]**: FORMAT key to use for obtaining ALT allele count (**DEFAULT**: `AC`)
* `-C`, `--AC_key_new` **[str]**: FORMAT key to use for new (filtered) ALT allele count (**DEFAULT**: `NAC`)
* `-a`, `--AF_key` **[str]**: FORMAT key to use for ALT allele frequency (**DEFAULT**: `AF`)
* `-A`, `--AF_key_new` **[str]**: FORMAT key to use for new (filtered) ALT allele frequency (**DEFAULT**: `NAF`)
* `-d`, `--DP_key` **[str]**: FORMAT key to use for read depth (coverage) (**DEFAULT**: `DP`)
* `-t`, `--min_AC` **[int]**: allele count cutoff (min reads allowed to support allele, REF or ALT) (**DEFAULT**: `0`)
* `-q`, `--min_AF` **[float]**: minor allele frequency cutoff (min allowed, REF or ALT) (**DEFAULT**: `0`)
* `-Q`, `--max_AF` **[float]**: major allele frequency cutoff (max allowed, REF or ALT) (**DEFAULT**: `1`)
* `-D`, `--min_DP` **[int]**: read depth (coverage) cutoff (min allowed, REF or ALT) (**DEFAULT**: `1`)
* `-I`, `--INFO_rules` **[str]**: custom rules for filtering based on the INFO column, formatted as a comma-separated string in the case of multiple rules; must use keys present in the INFO column, and must use the comparison operators `==`/`!=`/`<`/`<=`/`>=`/`>` (**DEFAULT**: `None`)
	* Example: `--INFO_rules="STB>0.5,STB<0.9"`
* `-s`, `--sample_rules` **[str]**: custom rules for filtering based on the sample column, formatted as a comma-separated string in the case of multiple rules; must use keys present in the FORMAT column, and must use the comparison operators `==`/`!=`/`<`/`<=`/`>=`/`>` (**DEFAULT**: `None`)
	* Example: `--sample_rules="FSRF>=5,FSRR>=5,FSAF>=5,FSAR>=5"`

## <a name="examples"></a>EXAMPLES

### EXAMPLE 1: A Simple Run with Custom AC, AF, and DP Keys

At a minimum, five options are required: `--seq_len`, `--error_per_site`, `--num_samples`, `--FP_cutoff`, and `--VCF_files`.

By default, VCFgenie looks for allele count, allele frequency, and read depth data in the `FORMAT`/`sample` columns using the default data keys `AC`, `AF`, and `DP`, respectivley. If your analysis will **not** use these default data keys, you must provide both the required options and the keys to use for obtaining allele count (`--AC_key`), allele frequency (`--AF_key`), and read depth (`--DP_key`).

Our file `example_A.vcf` uses `FAO`, `AF`, and `FDP` (quality-corrected values) for these keys, so we will specify those as follows:

	VCFgenie.py --seq_len=7857 --error_per_site=0.01103 --num_samples=1 --FP_cutoff=0.05 --out_dir=ex1_out --AC_key=FAO --AF_key=AF --DP_key=FDP --VCF_files example_A.vcf > ex1.out

### EXAMPLE 2: Miniumum AC, AF, and DP Values

If, in addition to the specified FP cutoff (`--FP_cutoff`), you wish to specify a minimum acceptable allele count (`--min_AC`), minimum acceptable allele frequency (`--min_AF`), or minimum acceptable read depth (`--min_DP`), they may be specified as follows:

	VCFgenie.py --seq_len=7857 --error_per_site=0.01103 --num_samples=1 --FP_cutoff=0.05 --out_dir=ex2_out --AC_key=FAO --AF_key=AF --DP_key=FDP --min_AC=10 --min_AF=0.01 --min_DP=100 --VCF_files example_A.vcf > ex2.out

### EXAMPLE 3: INFO and Sample-Based Rules

To implement rules based on the `INFO` column, here specifying the value of `STB` (strand bias) must be both greater than 0.5 and less than 0.9, run as follows:

	VCFgenie.py --seq_len=7857 --error_per_site=0.01103 --num_samples=1 --FP_cutoff=0.05 --out_dir=ex3a_out --AC_key=FAO --AF_key=AF --DP_key=FDP --INFO_rules="STB>0.5,STB<0.9" --VCF_files example_A.vcf > ex3a.out
	
To implement rules based on the `<sample>` column, here specifying the allele counts for  REF and ALT alleles must be at least 5 on the forward (`FSRF` and `FSAF`) and reverse (`FSRR` and `FSAR`) strands, run as follows:

	VCFgenie.py --seq_len=7857 --error_per_site=0.01103 --num_samples=1 --FP_cutoff=0.05 --out_dir=ex3b_out --AC_key=FAO --AF_key=AF --DP_key=FDP --sample_rules="FSRF>=5,FSRR>=5,FSAF>=5,FSAR>=5" --VCF_files example_A.vcf > ex3b.out

### EXAMPLE 4: Multiple Input Files

To provide a specific list of more than one VCF files as input, run as follows:

	VCFgenie.py --seq_len=7857 --error_per_site=0.01103 --num_samples=1 --FP_cutoff=0.05 --out_dir=ex4a_out --AC_key=FAO --AF_key=AF --DP_key=FDP --VCF_files example_A.vcf example_B.vcf > ex4a.out

To use **all** VCF files in the working directory as input, run as follows:

	VCFgenie.py --seq_len=7857 --error_per_site=0.01103 --num_samples=1 --FP_cutoff=0.05 --out_dir=ex4b_out --AC_key=FAO --AF_key=AF --DP_key=FDP --VCF_files *.vcf > ex4b.out

### EXAMPLE 5: All Options

To use all currently available options, run as follows:

	VCFgenie.py --seq_len=7857 --error_per_site=0.01103 --num_samples=93 --FP_cutoff=1 --out_dir=ex5_out --AF_key=AF --AF_key_new==AFN --AC_key=FAO --AC_key_new=FAON --DP_key=FDP --min_DP=100 --min_AC=10 --min_AF=0 --max_AF=1 --INFO_rules="STB>0.5,STB<0.9" --sample_rules="FSRF>=5,FSRR>=5,FSAF>=5,FSAR>=5" --VCF_files *.vcf > ex5.out

## <a name="output"></a>Output

**VCFgenie** outputs the following data:

### <a name="standard-output"></a>Output 1: Standard Output

At the command line (Terminal), VCFgenie will report:

1. a log of the input parameters
2. the file(s) to be processed
3. data regarding each failing allele and MULTIALLELIC site
4. PER-FILE summary statistics
5. TOTAL summary statistics related to `AC`, `AF`, and `DP` of passing and failing alleles (original values, before correction)

Results will be reported separately for the following categories:

* `ALLELES`: all alleles
* `REF`: REF alleles
* `ALT`: ALT alleles
* `MAJOR`: MAJOR (most common) alleles
* `MINOR`: MINOR (not most common) alleles

Statistics include:

* `n`: number of alleles
* `min`: minimum
* `Q1`: first quartile (25th percentile)
* `mean`: mean
* `std`: standard deviation
* `median`: median (50th percentile)
* `Q3`: third quartile (75th percentile)
* `max`: maximum

Note that:

1. The number of MINOR alleles will always equal the number of ALT alleles, because there can only be one MAJOR allele at a site and there can also only be one REF allele at a site. However, although equal in number, the alleles in each group are not the same.
2. In rare instances, the `AF` of a MAJOR (most common) allele may be <50% at multiallelic sites. For example, suppose a site has the following alleles and allele counts: A (10 reads), C (700 reads), G (530 reads), and T (450 reads), for a total read depth of 1,690. Even though C is the major allele, its `AF` is only 41.4%.

### <a name="filtered-VCF-files"></a>Output 2: Filtered VCF Output File(s)

VCFgenie will produce one `*_filtered.vcf` file for each `*.vcf` file in the working directory, to be placed in `--out_dir` (DEFAULT: `VCFgenie_output`)

## <a name="todo"></a>TODO

Planned additions to VCFgenie include:

1. Implement a flag for overwriting `--AC_key_new` and `--AF_key_new` values in the `INFO` column
2. Allow an arbitrarily large number of samples (columns >10) to be automatically processed
3. Allow analyses to be based on INFO column
4. Implement more complex sequencing error models in the form of user-provided 4x4 (nucleotide) and 64x4 (trinucleotide-context) matrices
5. Graphical user interface (GUI)

## <a name="troubleshooting"></a>Troubleshooting

If you have questions about **VCFgenie**, please click on the <a target="_blank" href="https://github.com/chasewnelson/VCFgenie/issues">Issues</a> tab at the top of this page and begin a new thread, so that others might benefit from the discussion. Common questions will be addressed in this section.

## <a name="acknowledgments"></a>Acknowledgments
**VCFgenie** was written and is maintained with support from a Research Fellowship from the National Cancer Institute (NCI), National Institutes of Health (NIH) to C.W.N. (2021-present), Lisa Mirabello group. The logo image was designed by Mitch Lin (2022); copyright-free DNA helix obtained from Pixabay. This product is the result of work by Laurie Burdette, Lisa Mirabello, Sambit Mishra, Chase W. Nelson, Maisa Pinheiro, and Meredith Yeager.

## <a name="citation"></a>Citation

When using this software, please refer to and cite this page:

>https://github.com/chasewnelson/VCFgenie

## <a name="contact"></a>Contact
If you have questions about **VCFgenie**, please click on the <a target="_blank" href="https://github.com/chasewnelson/VCFgenie/issues">Issues</a> tab at the top of this page and begin a new thread, so that others might benefit from the discussion.

Other correspondence should be addressed to Chase W. Nelson: 

* chase.nelson <**AT**> nih <**DOT**> gov