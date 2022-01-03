<img src="https://github.com/chasewnelson/VCFgenie/blob/main/VCFgenie_logo.png?raw=true" title="VCFgenie logo by Mitch Lin" alt="VCFgenie logo by Mitch Lin" align="middle">

# VCFgenie
**VCFgenie** is a Python program for dynamically processing within-sample (pooled-sequencing) variants in a VCF file by implementing arbitrarily compplex filtering rules and controlling for a false-discovery rate.

To test the software with [example data](#examples), execute the program at the Unix command line or Mac Terminal as follows:

### FORMAT:

	VCFgenie.pl --genome_len=7857 --error_per_site=0.01103 --num_samples=1 --FDR=1 \
		--AF_key=AF --AC_key=FAO --DP_key=FDP --min_DP=100 --min_AC=10 --min_MAF=0.01 \
		--VCF_files example.vcf > filter_VCF.out

Find more [examples](#examples) below. 

## <a name="contents"></a>Contents

* [Description](#description)
* [How it Works](#how-it-works)
* [Options](#options)
* [Examples](#examples)
* [Output](#output)
* [Troubleshooting](#troubleshooting)
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)
* [Contact](#contact)
* [References](#references)

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

