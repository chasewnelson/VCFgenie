<img src="https://github.com/chasewnelson/VCFgenie/blob/main/VCFgenie_logo.png?raw=true" title="VCFgenie logo by Mitch Lin" alt="VCFgenie logo by Mitch Lin" align="middle">

# VCFgenie
**VCFgenie** is a Python program for dynamically processing within-sample (pooled-sequencing) variants in a VCF file by implementing arbitrarily complex filtering rules and controlling for a user-provided false-discovery rate.

To test the software with [example data](#examples), execute the program at the Unix command line or Mac Terminal as follows:

### FORMAT:

	VCFgenie.py --seq_len=7857 --error_per_site=0.01103 --num_samples=1 --FDR=0.05  --VCF_files example.vcf > example.out

Find more [examples](#examples) below. 

## <a name="contents"></a>Contents

* [Description](#description)
* [How it Works](#how-it-works)
* [Allele Filtering Rules](#allele-filtering-rules)
* [Options](#options)
* [Examples](#examples) - PENDING
* [Output](#output) - PENDING
* [Troubleshooting](#troubleshooting) - PENDING
* [Acknowledgments](#acknowledgments) - PENDING
* [Citation](#citation) - PENDING
* [Contact](#contact) - PENDING
* [References](#references) - PENDING

## <a name="description"></a>Description

<a target="_blank" rel="noopener noreferrer" href="http://samtools.github.io/hts-specs/VCFv4.2.pdf">**Variant Call Format (VCF)** files</a> are currently the primary means of storing genetic polymorphism data. These data may summarize either **(1) individual-based sequencing data**, wherein one sequencing reaction corresponds to one individual (any ploidy); or **(2) pooled-sequencing data**, wherein one sequencing reaction corresponds to genetic material from multiple individuals from a population under study (any ploidy). Common uses of pooled-sequencing include analysis of somatic tissue samples containing populations of tumor cells; samples containing known quantities of genetic material from multiple indivudals; and within-host populations of pathogens such as viruses.

VCF files are the result of single nucleotide polymorphism (SNP) calling, for which numerous tools exist, including Freebayes, Genome Analysis Tool Kit (GATK), Ion Proton Variant Caller (TVC), Intrahost Variant Analysis of Replicates (iVar), LowFreq, Samtools mpileup, and VarScan 2.

After SNP calling, various tools and libraries also exist for downstream analysis and filtering of VCF files. However, current options are designed for filtering but not modifying VCF records; treat reference (REF) alleles inconsistently; or do not provide options for handling sites that are multiallelic, i.e., more than one alternate (ALT) allele. Thus, there is a need for a tool that implements arbitrary user-provided filtering rules to both filter and modify records, handle failing REF and multiple ALT alleles consistently, and control for a false-discovery rate among multiple samples, specifically in a way that is reproducible (i.e., no custom code needed). VCFgenie was written to meet this need, and allows the aforementioned steps to be performed before downstream (e.g., evolutionary and population genetic) analyses that rely on the called polymorphisms.

## <a name="how-it-works"></a>How it Works
VCFgenie is written in Python. Its dependencies include the following Python libraries, which may need to be installed:

* `argparse`
* `collections`
* `numpy`
* `operator`
* `os`
* `re`
* `scipy`
* `sys`
* `typing`

The program examines one or more user-provided VCF file(s) to filter and modify records based on an arbitrarily complex set of user-provided rules. One record corresponds to one line of the VCF file, which reports a miniumum of two alleles: REF and ALT. Arbitrarily larger numbers of ALT alleles (i.e., MULTIALLELIC sites) are allowed. 

At present, VCFgenie explicitly examines the first sample reported in the file, i.e., the sample reported in column 10, with identifying keys provided in column 9 (FORMAT). Because VCF file conventions differ widely, VCFgenie allows the user to specify which FORMAT data key should be used for:

1. `allele count`: the number of reads matching to the ALT allele (**DEFAULT**: "AC")
2. `allele frequency`: the frequency of the ALT allele (**DEFAULT**: "AF")
3. `read depth (coverage)`: the number of mapping reads overlapping the site (**DEFAULT**: "DP")

Many VCF files contain only one (pooled) sample, in which case the above analysis will be sufficient. However, if the user wishes to analyze multiple sample columns from the same file, it is currently necessary to create a separate VCF file for each sample to be analyzed, with that sample present in column 10. Future instantiations of VCFgenie may allow (1) analysis based on AC/AF/DP data from the INFO column; and (2) automatic analysis of multiple samples, i.e., columns 10 and above.

In addition to AC/AF/DP data, additional user-provided filtering criteria may be based on any keys/values present in the INFO or FORMAT/sample columns. Criteria for accepting or rejecting each allele are evaluated in the order described in [Options](#options). VCFgenie stops analyzing an allele at its first point of failure. Thus, an allele that fails because the site has insufficient DP may also have failed other criteria if DP were sufficient. Thus, summary statistics reported in the output refer explicity to the first (highest-level) criterion that failed.

VCFgenie implements a dynamic (site-by-site) binomial false-discovery rate filter based on the <a target="_blank" rel="noopener noreferrer" href="https://github.com/chasewnelson/SARS-CoV-2-ORF3d#figure-8">filter_VCF.py</a> script described by <a target="_blank" rel="noopener noreferrer" href="https://elifesciences.org/articles/59633">Nelson et al. (2020)</a>.

The method is **dynamic** because each site is evaluated separately, with results depending on that site's `DP` and each allele's `AC`. The method is **binomial** because it determines the probability, given a sequencing error rate, that each allele's AC could have resulted from sequencing error alone. Finally, the method controls for a **false-discovery rate (FDR)** by utilize a user-defined maximum acceptable number of false variant calls, rejecting any variants that were too likely the result of error alone. This false discovery rate is a function of four parameters, each of which must or can be supplied by the user (see [Options](#options)):

1. `--error_per_site`: the given sequencing technology's error rate per base
2. `--DP_key`, corresponding to read depth (coverage) at a site (i.e., each **read** is another opportunity for error)
3. `--seq_len`: length of reference sequence (e.g., contig, chromosome, or genome) in nucleotides (i.e., each **site** is another opportunity for error)
4. `--num_samples`: number of samples examined (i.e., each **sample** is another opportunity for error).

The binomial calculation uses the `scipy.stats` function `binom.cdf`:

	binom.cdf(x, n, p)

where:

	n = total reads
	p = error rate / 3
	x = number of reads of minor allele
	FDR cutoff = (1 - binom.cdf(x - 1, n, p)) * seq_len * num_samples
	
Thus, the method assumes a single nucleotide variant (SNV) model in which all three possible nucleotide errors occur at the same rate. Future instantiations of VCFgenie will allow more complex error models to be specified.

To implement the above method, the user must supply a **FDR cutoff** (`--FDR_cutoff`). For example, if one sample is being analyzed, one might wish to use an FDR cutoff of 0.05. If multiple samples (e.g., 100) are being analyzed, one might wish to specify that a maximum mean of 1 false variant across the entire analysis is acceptable.

Note the advantage of using the FDR is that it precludes the need for a fixed arbitrary `AF` cutoff (although one may still be applied using `--min_MAF`).

## <a name="allele-filtering-rules"></a>Allele Filtering Rules

Resolution of failing alleles requires careful consideration. VCFgenie always maintains the same `DP` (coverage) at a site before and after filtering. However, allele counts (`AC`) may be reassigned to different alleles, depending on which alleles pass the filtering criteria.

Filtering rules are tested in the following order:

1. `failZeroAC`: the number of reads corresponding to the tested allele (REF or ALT) is already 0 (i.e., `AF` is 0%)
2. `failDP`: read depth (coverage) is insufficient at this site
3. `failAC`: minor allele count (whether REF or ALT) is insufficient at this site
4. `failMAF`: minor allele frequency (whether REF or ALT) fails `--min_MAF`
5. `failINFO`: fails one of the user-provided `--INFO_rules`
6. `failsample`: fails one of the user-provided `--sample_rules`
7. `failFDR`: fails because it implies too high a false-discovery rate

An allele also be given the flag `fixedALT`, `fixedREF`, or `pass` if apropriate.

**Usually**, most sites are biallelic (one REF and one ALT allele) and the REF is the major allele. In such cases, if the ALT allele fails, the ALT reads will be assumed to represent error, and will simply be re-assigned to the REF allele (i.e., the REF allele will be fixed at 100%). Complementarily, at a biallelic site where the ALT is the major allele and the REF fails, the REF reads will be reassigned to the ALT (i.e., the ALT allele will be fixed at 100%).

**However, in rare instances,** other unusual cases may occur, especially with pooled-sequencing data involving poorly-supported REF alleles and multiallelic sites. Some rules, peculiar instances, and methods of resolution are listed below:

1. After all alleles are evaluated, reads from failing alleles are distributed to the remaining (passing) alleles in proportion to their frequencies among the reads attributed to passing alleles.
2. If `Number=A` for a given key supplied in `--INFO_rules` or `--sample_rules`, then there is one value per ALT allele, and the rule is evaluated separately for each ALT (REF cannot be evaluated).
3. If `Number=1` for a given key supplied in `--INFO_rules` or `--sample_rules`, then there is one value for the site, and the rule is applied equally to all alleles (REF and ALT).
4. If a site fails (e.g., due to insufficient coverage), or if all alleles fail at that site (including REF), then all reads are given to the **major** (most common) allele, whether REF or ALT.
5. If an allele fails at a site, then the **first** rule it failed will be reported in the `FILTER` column.
6. For a site with one or more failing alleles, **one flag per allele**, including reference, will be reported in the `FILTER` column (semicolon-separated). This will including `PASS` for any passing alleles. Distinct values will be reported only once, in alphabetical order (order is not meaningful). Thus a site that has one ALT that fails will have two flags: one for the REF, one for the ALT.
7. VCFgenie will add the following data to the `INFO` column:

	- `DECISION`=(`fixedREF`|`fixedALT`|`fail`|`failZeroAC`|`failDP`|`failAC`|`failMAF`|`failINFO`|`failsample`|`failFDR`|`pass`)
	- `STATUS`=(`PASS`|`FAIL`)
	- `MULTIALLELIC` (Flag)
	- `FAIL_REF` (Flag)

8. Any flags other than `PASS` that were present in the `FILTER` column of the VCF input will be kept.
9. `NOCALL`: this flag is added to any site at which a MAJOR ALLELE fails. At the same time, any `PASS` flags from previous passing (minor) alleles are removed. In other words, `PASS` is converted to `NOCALL` at sites where a major allele fails. It is possible for a minor allele to pass while a major allele fails, e.g., when a `--sample_rule` specifies a criterion using a key with `Number=A` (only ALT alleles have a value) and an ALT allele is the major allele. In such cases, **no** `AF` **corrections are performed**; `AC` and `AF` are left as they were prior to evaluation.

## <a name="options"></a>Options
Call **VCFgenie** using the following options: 

**REQUIRED:**

* `-i`, `--VCF_files` **[FILE(S)]**: input variant call format (VCF) file(s)
* `-e`, `--error_per_site` **[float]**: sequencing error rate per site (assumes all nucleotides equally probable)
* `-L`, `--seq_len` **[int]**: length of reference sequence (e.g., contig, chromosome, or genome) in nucleotides
* `-n`, `--num_samples` **[int]**: number of samples (VCF files) in full analysis
* `-f`, `--FDR_cutoff` **[float]**: analysis-wide false discovery rate (FDR) cutoff

**OPTIONAL:**

* `-o`, `--out_dir` **[DIR]**: output directory name
* `-a`, `--AC_key` **[str]**: FORMAT key to use for obtaining variant allele count
* `-A`, `--AC_key_new` **[str]**: FORMAT key to use for new (filtered) variant allele count
* `-v`, `--AF_key` **[str]**: FORMAT key to use for variant allele frequency
* `-V`, `--AF_key_new` **[str]**: FORMAT key to use for new (filtered) variant allele frequency
* `-D`, `--DP_key` **[str]**: FORMAT key to use for read depth (coverage)
* `-d`, `--min_DP` **[int]**: read depth (coverage) cutoff (min allowed)
* `-c`, `--min_AC` **[int]**: allele count cutoff (min reads allowed to support minor allele)
* `-m`, `--min_MAF` **[float]**: minor allele frequency cutoff (min allowed)
* `-I`, `--INFO_rules` **[str]**: custom rules for filtering based on the INFO column, formatted as a comma-separated string in the case of multiple rules; must use keys present in the INFO column, and must use the comparison operators `==`/`!=`/`<`/`<=`/`>=`/`>`
	* Example: `--INFO_rules="STB>0.5,STB<0.9"`
* `-s`, `--sample_rules` **[str]**: custom rules for filtering based on the sample column, formatted as a comma-separated string in the case of multiple rules; must use keys present in the FORMAT column, and must use the comparison operators `==`/`!=`/`<`/`<=`/`>=`/`>`
	* Example: `--sample_rules="FSRF>=5,FSRR>=5,FSAF>=5,FSAR>=5"`

*To be continued...*

## <a name="examples"></a>EXAMPLES

### EXAMPLE 1: A Simple Run

If your analysis will use the VCF file's default data keys for allele count (`AC`), allele frequency (`AF`), and read depth (`DP`), you need only provide the required options:

	VCFgenie.py --seq_len=7857 --error_per_site=0.01103 --num_samples=1 --FDR=0.05  --VCF_files example1.vcf > example1.out

### EXAMPLE 2: Custom AC, AF, and DP Keys

If your analysis will **not** use the VCF file's default data keys for allele count (`AC`), allele frequency (`AF`), and read depth (`DP`), you must provide both the required options and the keys to use for obtaining read depth (`--DP_key`), allele count (`--AC_key`), and allele frequency (`--AF`):

	VCFgenie.pl --seq_len=7857 --error_per_site=0.01103 --num_samples=1 --FDR=0.05 --AC_key=FAO --AF_key=AF --DP_key=FDP --VCF_files example.vcf > VCFgenie.out

### EXAMPLE 3: Miniumum AC, MAF, and DP Values

If, in addition to the specified FDR cutoff (`--FDR`), you wish to specify a minimum acceptable allele count (`--min_AC`), minimum acceptable minor allele frequency (`--min_MAF`), or minimum acceptable read depth (`--min_DP`), they may be specified as follows:

	VCFgenie.pl --seq_len=7857 --error_per_site=0.01103 --num_samples=1 --FDR=0.05 --AC_key=FAO --AF_key=AF --DP_key=FDP --min_AC=10 --min_MAF=0.01 --min_DP=100 --VCF_files example.vcf > VCFgenie.out

### EXAMPLE 4: To Be Continued...

*To be continued...*

## <a name="output"></a>Output

**VCFgenie** outputs the following data:

### <a name="standard-output"></a>Standard Output

At the command line (Terminal), VCFgenie will report a log of the input parameters, the file(s) to be processed, and various summary statistics related to passing and failing alleles.

### <a name="filtered-VCF-files"></a>Filtered VCF Output File(s)

VCFgenie will produce one `*_filtered.vcf` file for each `*.vcf` file in the working directory, to be placed in `--out_dir` (DEFAULT: `VCFgenie_output`)


## <a name="troubleshooting"></a>Troubleshooting

If you have questions about **VCFgenie**, please click on the <a target="_blank" href="https://github.com/chasewnelson/VCFgenie/issues">Issues</a> tab at the top of this page and begin a new thread, so that others might benefit from the discussion. Common questions will be addressed in this section.

## <a name="acknowledgments"></a>Acknowledgments
**VCFgenie** was written and is maintained with support from a Research Fellowship from the National Cancer Institute (NCI), National Institutes of Health (NIH) to C.W.N. (2021-present), Lisa Mirabello group. The logo image was designed by Mitch Lin (2019); copyright-free DNA helix obtained from Pixabay. Thanks to Laurie Burdette, Lisa Mirabello, Sambit Mishra, Maisa Pinheiro, and Meredith Yeager for discussion.

## <a name="citation"></a>Citation

When using this software, please refer to and cite this page:

>https://github.com/chasewnelson/VCFgenie

## <a name="contact"></a>Contact
If you have questions about **VCFgenie**, please click on the <a target="_blank" href="https://github.com/chasewnelson/VCFgenie/issues">Issues</a> tab at the top of this page and begin a new thread, so that others might benefit from the discussion.

Other correspondence should be addressed to Chase W. Nelson: 

* chase.nelson <**AT**> nih <**DOT**> gov