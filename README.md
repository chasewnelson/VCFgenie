<img src="https://github.com/chasewnelson/VCFgenie/blob/master/VCFgenie_logo.png?raw=true" title="VCFgenie logo by Mitch Lin" alt="VCFgenie logo by Mitch Lin" align="middle">

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

**Variant Call Format (VCF)** files are currently the primary means of storing genetic polymorphism data. These data may summarize either **(1) individual-based sequencing data**, wherein one sequencing reaction corresponds to one individual (any ploidy); or **(2) pooled-sequencing data**, wherein one sequencing reaction corresponds to genetic material from multiple individuals from a population under study (any ploidy). Common uses of pooled-sequencing include analysis of somatic tissue samples containing populations of tumor cells; samples containing known quantities of genetic material from multiple indivudals; and within-host populations of pathogens such as viruses.

VCF files are the result of single nucleotide polymorphism (SNP) calling, for which numerous tools exist, including Freebayes, Genome Analysis Tool Kit (GATK), Ion Proton Variant Caller (TVC), Intrahost Variant Analysis of Replicates (iVar), LowFreq, Samtools mpileup, and VarScan 2.

After SNP calling, various tools and libraries also exist for downstream analysis and filtering of VCF files. However, current options are designed for filtering but not modifying VCF records; treat reference (REF) alleles inconsistently; or do not provide options for handling sites that are multiallelic, i.e., more than one alternate (ALT) allele. Thus, there is a need for a tool that implements arbitrary user-provided filtering rules to both filter and modify records, handle failing REF and multiple ALT alleles consistently, and control for a false-discovery rate among multiple samples, specifically in a way that is reproducible (i.e., no custom code needed). VCFgenie was written to meet this need, and allows the aforementioned steps to be performed before downstream (e.g., evolutionary and population genetic) analyses that rely on the called polymorphisms.