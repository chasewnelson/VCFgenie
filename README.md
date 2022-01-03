<img src="https://github.com/chasewnelson/VCFgenie/blob/main/VCFgenie_logo.png?raw=true" title="VCFgenie logo by Mitch Lin" alt="VCFgenie logo by Mitch Lin" align="middle">

# VCFgenie
**VCFgenie** is a Python program for dynamically processing within-sample (pooled-sequencing) variants in a VCF file by implementing arbitrarily compplex filtering rules and controlling for a false-discovery rate.

To test the software with the [example data](#examples), execute the program at the Unix command line or Mac Terminal as follows:

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

[**Variant Call Format (VCF)** files](#http://samtools.github.io/hts-specs/VCFv4.2.pdf) are currently the primary means of storing genetic polymorphism data. These data may summarize either **(1) individual-based sequencing data**, wherein one sequencing reaction corresponds to one individual (any ploidy); or **(2) pooled-sequencing data**, wherein one sequencing reaction corresponds to genetic material from multiple individuals from a population under study (any ploidy). Common uses of pooled-sequencing include analysis of somatic tissue samples containing populations of tumor cells; samples containing known quantities of genetic material from multiple indivudals; and within-host populations of pathogens such as viruses.

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

The program examines one or more user-provided VCF file(s) to filter and modify records based on an arbitrarily complex set of user-provided rules. One record corresponds to one line of the VCF file, which reports a miniumum of two alleles: REF and ALT. Arbitrarily larger numbers of ALT alleles (i.e., MULTIALLELIC sites) are allowed. Criteria for accepting or rejecting each allele are evaluated in the order described in [Options](#options).

At present, VCFgenie explicitly e





VCF file conventions differ widely. Thus, VCFgenie allows the user to specify which FORMAT data key should be used for:

1. `allele count`: the number of reads matching to the ALT allele (DEFAULT: "AC")
2. `allele frequency`: the frequency of the ALT allele (DEFAULT: "AF")
3. `read depth (coverage)`: the number of mapping reads overlapping the site (DEFAULT: "DP")




After reading in the user-provided alignment, OLGenie calculates the number of NN, SN, NS, and SS sites and differences, reporting the mean of all pairwise comparisons. This is done separately for each focal reference codon by considering all unique nonamer (9nt) alleles of which the reference codon is the center, and of which 6nt constitute a minimum overlapping unit: one reference gene codon and its two overlapping alternate gene codons. (Note that sas13 is unique in that one reference codon overlaps exactly one alternate codon.) OLGenie is sufficiently fast that these tasks require no parallelism beyond the level of the single gene alignment. Thus, for datasets with many genes, the user can implement their own parallelization by running numerous alignments (genes) simultaneously.

After results are obtained for each focal codon in the alignment, significant deviations from the null expectation of neutrality (*d*<sub>N</sub> - *d*<sub>S</sub> = 0) may be tested using a *Z*-test, where the standard error is estimated using bootstrapping (focal codon unit). Don't worry — we provide scripts to do it all!