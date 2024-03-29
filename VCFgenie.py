#! /usr/bin/env python

"""
Purpose: Dynamically filter within-sample (pooled sequencing) variants to control false-positive count
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/SARS-CoV-2-ORF3d
Date   : 2021-12-04

Details: The script applies a dynamic (site-by-site) filter using the binomial
    distribution and the user-provided parameters. Must be called in a working
    directory with the target VCF files, automatically detected using their
    '.vcf' extension. All VCF files in the directory will be used.

    The user-provided parameter num_samples may correspond to the number of VCF
    files in the directory if none others have been analyzed in the same
    analysis; alternatively, the parameter may exceed the number of VCF files
    if other samples have been analyzed as part of the same analysis.

    The user may specify which FORMAT keys to use for read depth (coverage)
    (DP) and allele frequency (AF). ALT allele counts (AC) will be inferred.

    Multiallelic records (more than one ALT allele) will be printed across
    multiple lines, one per ALT.

    Each VCF file must contain the called SNPs for exactly 1 pooled sequencing
    sample, and therefore must contain exactly 8-10 columns:
        (1) CHROM
        (2) POS
        (3) ID
        (4) REF
        (5) ALT
        (6) QUAL
        (7) FILTER
        (8) INFO
        (9) FORMAT [OPTIONAL]
        (10) <sample_name> [OPTIONAL]

    The binomial filter uses the scipy function:

        binom.cdf(x, n, p)

    where:
        n = total reads
        p = error rate / 3
        x = number of reads of allele
        FP cutoff = (1 - binom.cdf(x - 1, n, p)) * seq_len * num_samples

    Note that x - 1 is used because we subtract the cumulative probability
        BEFORE the given number of variants.

    The equivalent function in base R is: pbinom(x, n, p)

    Definitions:
        -REF = reference sequence against which variants were called
        -ALT = alternative allele
        -AF = allele frequency of alternative allele
        -MAF = minor allele frequency
        -DP = depth (coverage)
        -fixedREF = site fixed for the REF allele (100% reference)
        -fixedALT = site fixed for a particular ALT allele (100% non-reference)
        -fail = allele fails criteria (all-inclusive)
            -failDP = read depth (coverage) fails min_DP
            -failZeroAC = allele fails when 0 reads support it
            -failAC = allele count (whether REF or ALT) fails min_AC
            -failMinAF = allele frequency (whether REF or ALT) fails min_AF
            -failMaxAF = allele frequency (whether REF or ALT) fails max_AF
            -failINFO: one or more of the user-provided INFO_rules fails
            -failsample: one or more of the user-provided sample_rules fails
            -failp = allele fails p_cutoff
        -pass = allele passes all criteria

TODO:
    -whether INFO or INDIVIDUAL SAMPLES examined

EXAMPLES:
    $ VCFgenie.py --error_rate=0.0001 --VCF_files=example.vcf
    $ VCFgenie.py --error_rate=0.0001 --p_cutoff=0.000001 --AC_key=FAO --AC_key_new=FAO --AF_key=AF --AF_key_new=AF --DP_key=FDP --overwrite_INFO --VCF_files=example.vcf

Returns:
    1) summary statistics about the variants that were processed (STDOUT)
    2) one *_filtered.vcf file for each *.vcf file in the working directory (to --out_dir)
"""

import argparse
import numpy as np
import operator as op
import os
import re
import sys
import vcf
from collections import defaultdict
from scipy.stats import binom
from typing import Any, Dict, List, NamedTuple, Optional, TextIO

usage = """# -----------------------------------------------------------------------------
VCFgenie.py - Dynamically filter within-sample (pooled sequencing) variants to control false positive calls
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ VCFgenie.py --help
    $ pydoc ./VCFgenie.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ VCFgenie.py --error_rate=0.0001 --VCF_files=example.vcf
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    # REQUIRED
    VCF_files: List[str]
    error_rate: float
    # seq_len: int
    # num_samples: int
    # FP_cutoff: float

    # OPTIONAL
    out_dir: str
    p_cutoff: float
    AC_key: str
    AC_key_new: str
    AF_key: str
    AF_key_new: str
    DP_key: str
    PVR_key_new: str
    PVA_key_new: str
    min_AC: float
    min_AF: float
    max_AF: float
    min_DP: float
    INFO_rules: str
    sample_rules: str
    overwrite_INFO: bool


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Dynamically filter variants to control false-positive count | EXAMPLE: VCFgenie.py --error_rate=0.0001 --VCF_files=example.vcf',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-i',
                        '--VCF_files',
                        metavar='FILE',
                        help='input, one or more variant call format (VCF) files [REQUIRED]',
                        required=True,  # 'required' won't work for positionals
                        nargs='+',
                        # type=argparse.FileType('r'))  # was 'rt' but want .gz as option?
                        type=str)

    parser.add_argument('-e',
                        '--error_rate',
                        metavar='float',
                        help='sequencing error rate for single nucleotide substitutions per base sequenced [REQUIRED]',
                        required=True,
                        type=float)

    # help='error rate per site expected as the result of sample preparation, library preparation, '
    # 'and sequencing (currently assumes all nucleotide errors equally probable) [REQUIRED]',

    # parser.add_argument('-L',
    #                     '--seq_len',
    #                     metavar='int',
    #                     help='length of reference sequence (e.g., contig, chromosome, or genome) in nucleotides '
    #                          '[REQUIRED]',
    #                     required=True,
    #                     type=int)

    # parser.add_argument('-n',
    #                     '--num_samples',
    #                     metavar='int',
    #                     help='number of samples (VCF files) in full analysis [REQUIRED]',
    #                     required=True,
    #                     type=int)

    # parser.add_argument('-f',
    #                     '--FP_cutoff',
    #                     metavar='float',
    #                     help='analysis-wide false-positive (FP) count cutoff [REQUIRED]',
    #                     required=True,
    #                     type=float)

    # -------------------------------------------------------------------------
    # OPTIONAL - WITH DEFAULT
    parser.add_argument('-o',
                        '--out_dir',
                        metavar='DIR',
                        help='output directory name [OPTIONAL]',
                        required=False,
                        type=str,
                        default='VCFgenie_output')

    parser.add_argument('-p',
                        '--p_cutoff',
                        metavar='float',
                        help='p-value cutoff [OPTIONAL]',
                        required=False,
                        type=float,
                        default=1.)

    parser.add_argument('-c',
                        '--AC_key',
                        metavar='str',
                        help='FORMAT key to use for obtaining ALT allele count [OPTIONAL]',
                        required=False,
                        type=str,
                        default='AC')  # 'FAO' for Ion Torrent

    parser.add_argument('-C',
                        '--AC_key_new',
                        metavar='str',
                        help='FORMAT key to use for new (filtered) ALT allele count [OPTIONAL]',
                        required=False,
                        type=str,
                        default='NAC')

    parser.add_argument('-a',
                        '--AF_key',
                        metavar='str',
                        help='FORMAT key to use for ALT allele frequency [OPTIONAL]',
                        required=False,
                        type=str,
                        default='AF')  # 'AF' remains for Ion Torrent

    parser.add_argument('-A',
                        '--AF_key_new',
                        metavar='str',
                        help='FORMAT key to use for new (filtered) ALT allele frequency [OPTIONAL]',
                        required=False,
                        type=str,
                        default='NAF')

    parser.add_argument('-d',
                        '--DP_key',
                        metavar='str',
                        help='FORMAT key to use for read depth (coverage) [OPTIONAL]',
                        required=False,
                        type=str,
                        default='DP')  # 'FDP' for Ion Torrent

    parser.add_argument('-v',
                        '--PVR_key_new',
                        metavar='str',
                        help='INFO and FORMAT key to use for p-value for REF allele [OPTIONAL]',
                        required=False,
                        type=str,
                        default='PVR')

    parser.add_argument('-V',
                        '--PVA_key_new',
                        metavar='str',
                        help='INFO and FORMAT key to use for p-values for ALT allele(s) [OPTIONAL]',
                        required=False,
                        type=str,
                        default='PVA')

    parser.add_argument('-t',
                        '--min_AC',
                        metavar='int',
                        help='allele count cutoff (min reads allowed to support allele, REF or ALT) [OPTIONAL]',
                        required=False,
                        type=float,  # OK to be a float
                        default=0)

    parser.add_argument('-q',
                        '--min_AF',
                        metavar='float',
                        help='minor allele frequency cutoff (min allowed, REF or ALT) [OPTIONAL]',
                        required=False,
                        type=float,
                        default=0)

    parser.add_argument('-Q',
                        '--max_AF',
                        metavar='float',
                        help='major allele frequency cutoff (max allowed, REF or ALT) [OPTIONAL]',
                        required=False,
                        type=float,
                        default=1)

    parser.add_argument('-D',
                        '--min_DP',
                        metavar='int',
                        help='read depth (coverage) cutoff (min allowed, REF or ALT) [OPTIONAL]',
                        required=False,
                        type=float,  # OK to be a float
                        default=1)

    parser.add_argument('-I',
                        '--INFO_rules',
                        metavar='str',
                        help='custom rules for filtering based on the INFO column, formatted as a comma-'
                             'separated string in the case of multiple rules; must use keys present in the '
                             'INFO column, and must use the comparison operators "=="/"!="/"<"/"<="/">="/">" '
                             '(e.g., --INFO_rules="STB>0.5,STB<0.9") [OPTIONAL]',
                        required=False,
                        type=str,
                        default=None)

    parser.add_argument('-s',
                        '--sample_rules',
                        metavar='str',
                        help='custom rules for filtering based on the sample column, formatted as a comma-'
                             'separated string in the case of multiple rules; must use keys present in the '
                             'FORMAT column, and must use the comparison operators "=="/"!="/"<"/"<="/">="/">" '
                             '(e.g., --sample_rules="FSRF>=5,FSRR>=5,FSAF>=5,FSAR>=5") [OPTIONAL]',
                        required=False,
                        type=str,
                        default=None)

    parser.add_argument('-O',
                        '--overwrite_INFO',
                        help='for any data keys matching between INFO and FORMAT, OVERWRITE INFO data with <sample> '
                             'data [OPTIONAL]',
                        action='store_true')

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments

    # REQUIRED

    # VCF files EXIST
    if bad_files := [file for file in args.VCF_files if not os.path.isfile(file)]:
        parser.error(f'\n### ERROR: VCF file(s) do not exist: {",".join(bad_files)}')

    # VCF files have RECOGNIZABLE HEADER
    re_start_pound = re.compile(r'^#')
    re_header_start = re.compile(r'^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
    for this_VCF_file in args.VCF_files:
        missing_header = True

        this_VCF_fh = open(this_VCF_file)
        for line in this_VCF_fh:
            if not re_start_pound.match(line):
                break

            if re_header_start.match(line):
                missing_header = False
                break

        if missing_header:
            parser.error(f'\n### ERROR: VCF_file="{this_VCF_file}" does not contain a recognizable header. '
                         'Wrong format?\n')

    # VCF files RECOGNIZABLE to pyvcf module
    # TODO validate AF_key and DP_key here?
    for this_VCF_file in args.VCF_files:
        try:
            this_vcf_Reader = vcf.Reader(filename=this_VCF_file)  # opens without error if exists, even if not VCF
        except AttributeError:
            parser.error(f'\n### ERROR: VCF_file="{this_VCF_file}" not readable by pyvcf. VCF file wrong format?\n')

        try:
            next(this_vcf_Reader)
        except ValueError as e:
            parser.error(f'### ERROR: ValueError: {e}\n'
                         f'### Records in VCF_file="{this_VCF_file}" not readable by pyvcf.\n'
                         '### VCF file wrong format?\n')
        except StopIteration as e:
            parser.error(f'\n### ERROR: StopIteration: {e}\n' 
                         f'### VCF_file="{this_VCF_file}" not iterable by pyvcf.\n'
                         '### VCF may not contain any records (rows;  FAOvariants).\n')

    # Appropriate numeric values
    if args.error_rate < 0 or args.error_rate >= 1 or not args.error_rate:
        parser.error(f'\n### ERROR: error_rate "{args.error_rate}" must be >= 0 and < 1')

    # if args.seq_len <= 0 or not args.seq_len:
    #     parser.error(f'\n### ERROR: seq_len "{args.seq_len}" must be > 0 nucleotides')
    #
    # if args.num_samples <= 0 or not args.num_samples:
    #     parser.error(f'\n### ERROR: num_samples "{args.num_samples}" must be > 0')

    # if args.FP_cutoff <= 0 or not args.FP_cutoff:
    #     parser.error(f'\n### ERROR: FP_cutoff "{args.FP_cutoff}" must be > 0')

    # OPTIONAL
    if os.path.isdir(args.out_dir):
        # sys.exit(f'\n### ERROR: out_dir "{args.out_dir}" already exists')
        parser.error(f'\n### ERROR: out_dir="{args.out_dir}" already exists')
    else:
        os.makedirs(args.out_dir)

    if args.p_cutoff < 0 or args.p_cutoff > 1 or not args.p_cutoff:
        parser.error(f'\n### ERROR: p_cutoff "{args.p_cutoff}" must be between 0 and 1 (inclusive)')

    # if args.AC_key == args.AC_key_new:
    #     # parser.error(f'AC_key="{args.AC_key}" may not match AC_key_new="{args.AC_key_new}"')
    #     print(f'\n### WARNING: AC_key="{args.AC_key_new}" already exists and will be OVERWRITTEN')
    #
    # if args.AF_key == args.AF_key_new:
    #     # parser.error(f'AF_key="{args.AF_key}" may not match AF_key_new="{args.AF_key_new}"')
    #     print(f'\n### WARNING: AF_key="{args.AF_key_new}" already exists and will be OVERWRITTEN')

    if args.min_AC < 0:
        parser.error(f'\n### ERROR: min_AC "{args.min_AC}" must be >= 0')

    if args.min_AF < 0:
        parser.error(f'\n### ERROR: min_AF "{args.min_AF}" must be >= 0')

    if args.max_AF > 1:
        parser.error(f'\n### ERROR: min_AF "{args.max_AF}" must be <= 1')

    if args.min_DP < 1:
        parser.error(f'\n### ERROR: min_DP "{args.min_DP}" must be >= 1')

    # INFO_rules and sample_rules will be validated first thing in main() when they're regexed

    return Args(VCF_files=args.VCF_files,
                error_rate=args.error_rate,
                # seq_len=args.seq_len,
                # num_samples=args.num_samples,
                # FP_cutoff=args.FP_cutoff,
                out_dir=args.out_dir,
                p_cutoff=args.p_cutoff,
                AC_key=args.AC_key,
                AC_key_new=args.AC_key_new,
                AF_key=args.AF_key,
                AF_key_new=args.AF_key_new,
                DP_key=args.DP_key,
                PVR_key_new=args.PVR_key_new,
                PVA_key_new=args.PVA_key_new,
                min_AC=args.min_AC,
                min_AF=args.min_AF,
                max_AF=args.max_AF,
                min_DP=args.min_DP,
                INFO_rules=args.INFO_rules,
                sample_rules=args.sample_rules,
                overwrite_INFO=args.overwrite_INFO)


# -----------------------------------------------------------------------------
def main() -> None:
    """ The best criticism of the bad is the practice of the better """

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    VCF_files = args.VCF_files
    error_rate = args.error_rate
    # seq_len = args.seq_len
    # num_samples = args.num_samples
    # FP_cutoff = args.FP_cutoff
    out_dir = args.out_dir  # optional
    p_cutoff = args.p_cutoff  # optional
    AC_key = args.AC_key  # optional
    AC_key_new = args.AC_key_new  # optional
    AF_key = args.AF_key  # optional
    AF_key_new = args.AF_key_new  # optional
    DP_key = args.DP_key  # optional
    PVR_key_new = args.PVR_key_new  # optional
    PVA_key_new = args.PVA_key_new  # optional
    min_AC = args.min_AC  # optional
    min_AF = args.min_AF  # optional
    max_AF = args.max_AF  # optional
    min_DP = args.min_DP  # optional
    INFO_rules = args.INFO_rules  # optional
    sample_rules = args.sample_rules  # optional
    overwrite_INFO = args.overwrite_INFO  # optional

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT AND LOG
    print(usage)

    print('# -----------------------------------------------------------------------------')
    # print('LOG')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')
    print(f'LOG:error_rate="{error_rate}"')
    # print(f'LOG:seq_len="{seq_len}"')
    # print(f'LOG:num_samples="{num_samples}"')
    # print(f'LOG:FP_cutoff="{FP_cutoff}"')
    print(f'LOG:out_dir="{out_dir}"')
    print(f'LOG:p_cutoff="{p_cutoff}"')
    print(f'LOG:AC_key="{AC_key}"')
    print(f'LOG:AC_key_new="{AC_key_new}"')
    print(f'LOG:AF_key="{AF_key}"')
    print(f'LOG:AF_key_new="{AF_key_new}"')
    print(f'LOG:DP_key="{DP_key}"')
    print(f'LOG:PVR_key_new="{PVR_key_new}"')
    print(f'LOG:PVA_key_new="{PVA_key_new}"')
    print(f'LOG:min_AC="{min_AC}"')
    print(f'LOG:min_AF="{min_AF}"')
    print(f'LOG:max_AF="{max_AF}"')
    print(f'LOG:min_DP="{min_DP}"')
    print(f'LOG:INFO_rules="{INFO_rules}"')
    print(f'LOG:sample_rules="{sample_rules}"')
    print(f'LOG:overwrite_INFO="{overwrite_INFO}"', flush=True)
    if overwrite_INFO:
        print('### WARNING: INFO will be overwritten with <sample> data when possible!')

    # WARN about unintended overwriting due to the specified new keys
    if AC_key == AC_key_new:
        print(f'### WARNING: AC_key_new="{AC_key_new}" already exists as AC_key="{AC_key_new}" and will be OVERWRITTEN')

    if AF_key == AF_key_new:
        print(f'### WARNING: AF_key_new="{AF_key_new}" already exists as AF_key="{AF_key_new}" and will be OVERWRITTEN')

    # TODO: check that ANY of the keys being used aren't already present? Perhaps overkill

    # -------------------------------------------------------------------------
    # REGEX & TUPLES

    # match user-supplied rules
    re_rule = re.compile(r'(\w+)([!=<>]+)([.\d]+)')

    # VCF file regex
    # NOTE: DP and AF are sufficient because the allele count, so inconsistently coded across VCFs, can simply be inferred

    # VCF_line_pattern = r"([\w\_\.]+)\t(\d+)\t([^\s]+)\t([acgtuACGTU]+)\t([acgtuACGTU]+)\t(\d+)\t(\w+)\t([\w\=\;\.\,]+)"
    # re_VCF_INFO = re.compile(VCF_INFO_pattern)

    re_VCF_INFO = re.compile(r'(\w+)=([\w\d/,.\-<>]+)')

    # VCF_SB_pattern = r"SB=\d+" # SB=0
    # re_VCF_SB = re.compile(VCF_SB_pattern)
    # VCF_DP4_grppattern = r"DP4=(\d+),(\d+),(\d+),(\d+)" # DP4=27,235,0,4
    # VCF_DP4_grpregex = re.compile(VCF_DP4_pattern)

    # VCF_DP_pattern = r"DP=\d+"  # DP=266
    # re_VCF_DP = re.compile(VCF_DP_pattern)
    # VCF_AF_pattern = r"AF=[\d\.]+"  # AF=0.015038
    # re_VCF_AF = re.compile(r'AF=[\d.,]+')

    # Data regex
    # re_datum = re.compile(r'[\d/.\-]+')
    # re_AF_datum = re.compile(r'([\s\t;])AF=[\d.NA]+')
    # re_DP_datum = re.compile(r'([\s\t;])DP=[\d.NA]+')
    # re_CSV_datum = re.compile(r'[,\d\w.\-]+')

    # VCF metadata block regex
    # re_metadata_block = re.compile(r'^##(\w+)=.+$')
    # ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
    # re_VCF_FORMAT_metadata = re.compile(r'^##(\w+)=<ID=(\w+),Number=(\w+),Type=(\w+),Description="(.+)">$')
    re_VCF_INFO_metadata = re.compile(r'^##(\w+)=<ID=(\w+),Number=(\w+),Type=(\w+),Description="(.+)".*>$')
    # TODO update regex for recording Source and Version:
    # ##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">

    # non-float, non-int regex
    # re_non_int = re.compile(r'\D')  # not a digit or dot or minus sign
    # re_non_float = re.compile(r'^\d^.^-')  # not a digit or dot or minus sign

    operator_choices = ('==', '!=', '<', '<=', '>=', '>')
    decision_choices = ('pass',
                        'fail', 'failZeroAC', 'failDP', 'failAC', 'failMinAF', 'failMaxAF', 'failINFO', 'failsample',
                        'failp',
                        'fixedREF', 'fixedALT',)  # 'failFP',
    decision_choices_FAIL = ('failZeroAC', 'failDP', 'failAC', 'failMinAF', 'failMaxAF', 'failINFO', 'failsample', 'failp')  # 'failFP'

    # -------------------------------------------------------------------------
    # PARSE INFO_rules and sample_rules into LISTS of 3-TUPLES using re_rule
    # rule_FILTER_metadata: List[str] = []
    INFO_rule_FILTER_dict: Dict[str, List[str]] = defaultdict(list)
    sample_rule_FILTER_dict: Dict[str, List[str]] = defaultdict(list)

    # INFO rules
    INFO_rule_lt = []
    if INFO_rules is not None:
        INFO_rule_list = INFO_rules.split(',')

        for i, INFO_rule in enumerate(INFO_rule_list):
            # key, operator, value = re_rule.match(INFO_rule).groups()
            INFO_rule_match = re_rule.match(INFO_rule)

            if INFO_rule_match is not None:
                INFO_rule_match_groups = INFO_rule_match.groups()
                key: str = INFO_rule_match_groups[0]
                operator: str = INFO_rule_match_groups[1]
                value: Any = INFO_rule_match_groups[2]

                if operator in operator_choices:
                    print(f'LOG:INFO_rule_{i}:key="{key}",operator="{operator}",value="{value}"')
                    INFO_rule_lt.append((key, operator, value))
                    INFO_rule_FILTER_dict[key].append(INFO_rule)
                    # rule_FILTER_metadata.append(f'##FILTER=<ID={},Description="User-provided INFO rule: {INFO_rule}">')
                else:
                    sys.exit(f'\n### ERROR: INFO_rule="{INFO_rule}" does not use an acceptable operator {operator_choices}')
            else:
                sys.exit(f'\n### ERROR: INFO_rule="{INFO_rule}" does not use an acceptable format or operator {operator_choices}')

        # Loop cleanup
        del INFO_rule_match_groups
        del key
        del operator
        del value

    # print(f'INFO_rule_lt={INFO_rule_lt}')

    # SAMPLE rules
    sample_rule_lt = []
    if sample_rules is not None:
        sample_rule_list = sample_rules.split(',')

        for i, sample_rule in enumerate(sample_rule_list):
            # key, operator, value = re_rule.match(sample_rule).groups()
            # print(f'LOG:sample_rule_{i}:key="{key}",operator="{operator}",value="{value}"')
            # sample_rule_lt.append((key, operator, value))

            # key, operator, value = re_rule.match(sample_rule).groups()
            sample_rule_match = re_rule.match(sample_rule)

            if sample_rule_match is not None:
                key, operator, value = sample_rule_match.groups()

                if operator in operator_choices:
                    print(f'LOG:sample_rule_{i}:key="{key}",operator="{operator}",value="{value}"')
                    sample_rule_lt.append((key, operator, value))
                    sample_rule_FILTER_dict[key].append(sample_rule)
                    # rule_FILTER_metadata.append(f'##FILTER=<ID=sample_rule_{i},Description="User-provided sample rule: {sample_rule}">')
                else:
                    sys.exit(
                        f'\n### ERROR: sample_rule="{sample_rule}" does not use an acceptable operator {operator_choices}')
            else:
                sys.exit(f'\n### ERROR: sample_rule="{sample_rule}" does not use an acceptable format or operator {operator_choices}')
    # print(f'sample_rule_lt={sample_rule_lt}')

    # VERIFY there are no duplicate keys; may have to deal with this later
    if not set(INFO_rule_FILTER_dict.keys()).isdisjoint(set(sample_rule_FILTER_dict.keys())):
        sys.exit(f'\n### ERROR: overlapping key(s) utilized for --INFO_rules ({set(INFO_rule_FILTER_dict.keys())}) and '
                 f'--sample_rules ({set(sample_rule_FILTER_dict.keys())}). Contact author for update')

    print()

    # NEW LINES TO PRINT TO METADATA
    new_FILTER_lines: List[str] = []

    for key in INFO_rule_FILTER_dict:
        new_FILTER_lines.append(f'##FILTER=<ID={key},Description="User-provided INFO rule(s): {",".join(INFO_rule_FILTER_dict[key])}">')

    for key in sample_rule_FILTER_dict:
        new_FILTER_lines.append(f'##FILTER=<ID={key},Description="User-provided sample rule(s): {",".join(sample_rule_FILTER_dict[key])}">')

    # -------------------------------------------------------------------------
    # INITIALIZE descriptions of summary statistics types
    summaryKey_to_description = {
        'fixedREF': 'site fixed for the REF allele (100% reference)',
        'fixedALT': 'site fixed for a particular ALT allele (100% non-reference)',
        'fail': 'allele fails criteria (all-inclusive)',
        'failZeroAC': 'allele fails when 0 reads support it (REF or ALT)',
        'failDP': 'read depth (coverage) fails --min_DP',
        'failAC': 'allele count (REF or ALT) fails --min_AC',
        'failMinAF': 'allele frequency (REF or ALT) fails --min_AF',
        'failMaxAF': 'allele frequency (REF or ALT) fails --max_AF',
        'failINFO': 'fails one or more of the user-provided --INFO_rules',
        'failsample': 'fails one or more of the user-provided --sample_rules',
        # 'failFP': 'fails --FP_cutoff',
        'failp': 'fails --p_cutoff',
        # 'failToFixedALT': 'failing minor allele matched REF; converted to fixed ALT allele',
        'pass': 'allele passes all criteria'
    }

    # # -------------------------------------------------------------------------
    # # Gather VCF files in working directory
    # VCF_filenames = sorted([name for name in os.listdir('.') if os.path.isfile(name) and name.endswith('.vcf')])
    # print('VCF files to examine, gathered from working directory: ' + str(VCF_filenames))
    # print('Will write output to files with names of the form: <VCF_name>_filtered.vcf')

    # -------------------------------------------------------------------------
    # Name VCF files to examine
    print('# -----------------------------------------------------------------------------')
    # re_file_ext = re.compile(r'.\w+$')
    print('VCF files to process (output files will have names of the form "<VCF_root_name>_filtered.vcf"): '
          '[IN_FILE_NAME] -> [OUT_FILE_NAME]')

    # Store input -> output file names
    VCF_to_out_file_dict: Dict[str, str] = defaultdict(str)
    VCF_files = sorted(VCF_files)
    # for this_VCF_file in sorted(VCF_files):
    for this_VCF_file in VCF_files:
        # Create out_file name
        this_inFile_root, this_inFile_ext = os.path.splitext(this_VCF_file)  # takes str, includes '.'

        # Input VCF file must end with '.vcf' extension and its implied out_file must not already exist
        if this_inFile_ext == '.vcf':
            this_outFile_name = this_inFile_root + '_filtered.vcf'

            if not os.path.exists(this_outFile_name):
                VCF_to_out_file_dict[this_VCF_file] = this_outFile_name
                print(f'{this_VCF_file} -> {this_outFile_name}')
            else:
                sys.exit(f'\n### ERROR: output VCF file {this_VCF_file} already exists\n')
        # elif this_inFile_ext == '.gz':
        #     this_outFile_name = this_inFile_root + '.gz_filtered.vcf'
        #
        #     if not os.path.exists(this_outFile_name):
        #         VCF_to_out_file_dict[this_VCF_file] = this_outFile_name
        #         print(f'{this_VCF_file} -> {this_outFile_name}')
        #     else:
        #         sys.exit(f'\n### ERROR: output VCF file {this_VCF_file} already exists\n')
        else:
            # sys.exit(f'\n### ERROR: input VCF file {this_VCF_file} must end with the extension ".vcf" or ".gz"\n')
            sys.exit(f'\n### ERROR: input VCF file {this_VCF_file} must end with the extension ".vcf"\n')

    print()

    # -------------------------------------------------------------------------
    # Prepare summary statistics
    sample_n = 0
    total_record_n = 0  # lines/rows

    total_ALLELE_n = 0
    total_ALLELE_pass_n = 0
    total_ALLELE_fail_n = 0

    total_REF_n = 0
    total_REF_pass_n = 0
    total_REF_fail_n = 0

    total_ALT_n = 0
    total_ALT_pass_n = 0
    total_ALT_fail_n = 0

    total_MAJOR_n = 0
    total_MAJOR_pass_n = 0
    total_MAJOR_fail_n = 0

    total_MINOR_n = 0
    total_MINOR_pass_n = 0
    total_MINOR_fail_n = 0

    # count_dd: Dict[str, Dict[str, int]] = defaultdict(dict)
    # AF_ddl: Dict[str, Dict[str, float]] = defaultdict(dict)
    # AC_ddl: Dict[str, Dict[str, int]] = defaultdict(dict)
    # DP_ddl: Dict[str, Dict[str, int]] = defaultdict(dict)

    ALLELE_n_dict: Dict[str, int] = defaultdict(int)
    ALLELE_AF_dl: Dict[str, List[float]] = defaultdict(list)
    ALLELE_AC_dl: Dict[str, List[int]] = defaultdict(list)
    ALLELE_DP_dl: Dict[str, List[int]] = defaultdict(list)

    REF_n_dict: Dict[str, int] = defaultdict(int)
    REF_AF_dl: Dict[str, List[float]] = defaultdict(list)
    REF_AC_dl: Dict[str, List[int]] = defaultdict(list)
    REF_DP_dl: Dict[str, List[int]] = defaultdict(list)

    ALT_n_dict: Dict[str, int] = defaultdict(int)
    ALT_AF_dl: Dict[str, List[float]] = defaultdict(list)
    ALT_AC_dl: Dict[str, List[int]] = defaultdict(list)
    ALT_DP_dl: Dict[str, List[int]] = defaultdict(list)

    MAJOR_n_dict: Dict[str, int] = defaultdict(int)
    MAJOR_AF_dl: Dict[str, List[float]] = defaultdict(list)
    MAJOR_AC_dl: Dict[str, List[int]] = defaultdict(list)
    MAJOR_DP_dl: Dict[str, List[int]] = defaultdict(list)

    MINOR_n_dict: Dict[str, int] = defaultdict(int)
    MINOR_AF_dl: Dict[str, List[float]] = defaultdict(list)
    MINOR_AC_dl: Dict[str, List[int]] = defaultdict(list)
    MINOR_DP_dl: Dict[str, List[int]] = defaultdict(list)

    # -------------------------------------------------------------------------
    # LOOP VCF files to EXTRACT SNP DATA
    for this_VCF_file in VCF_files:
        this_VCF_fh = open(this_VCF_file, 'rt')
        # if this_VCF_file.endswith('_filtered.vcf'):
        #     sys.exit(f'\n### ERROR: VCF file {this_VCF_file} may not already end with "_filtered.vcf", because '
        #              'this script will give it that suffix\n')

        # # Check file exists - already did this in arguments
        # if os.path.isfile(this_VCF_file):
        print('# -----------------------------------------------------------------------------')
        print(f'Analyzing file={this_VCF_file}')

        # ---------------------------------------------------------------------
        # OPEN the VCF file object
        # vcf_file_object = vcf.Reader(fsock=this_VCF_fh)  # CAN be .gz!, CAN be fh
        # this_vcf_reader = vcf.Reader(filename=this_VCF_file)

        # OPEN output file for writing
        # this_outFile_hdl = open(VCF_to_out_file_dict[this_VCF_fh.name], "w")
        this_out_file = os.path.join(out_dir, VCF_to_out_file_dict[this_VCF_file])
        this_out_fh = open(this_out_file, "wt")
        # this_vcf_writer = vcf.Writer(this_outFile_hdl, this_vcf_reader)

        # ---------------------------------------------------------------------
        # VCF DOCUMENTATION: http://samtools.github.io/hts-specs/VCFv4.2.pdf

        # INFO
        # ##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">
        # 'Number': an Integer that describes the number of values that can be included with the INFO field
        # 1=always a single number
        # 2=always a pair of numbers
        # A=the field has one value per alternate allele ALT  # TODO this is what we need for collapsing MULTIALLELICS back into one record
        # R=the field has one value for each possible allele (including the reference REF)
        # G=the field has one value for each possible genotype (more relevant to the FORMAT tags)
        # .=the number of possible values varies, is unknown, or is unbounded
        # 0=for ‘Flag’ type, for which the INFO field does not contain a Value entry
        # 'Description': must be surround by "double quotes"

        # Types for INFO fields are: Integer, Float, Flag, Character, and String

        # ##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">

        # FILTER
        # ## FILTER=<ID=ID,Description="description">

        # ---------------------------------------------------------------------
        # Prepare additional metadata entries (##), to print for each file directly above the header (#)
        # new_FILTER_lines = f'##FILTER=<ID=seq_len,Description="Reference sequence genome length (nucleotides): {seq_len}">\n' \
        #                    f'##FILTER=<ID=error_rate,Description="Sequencing error rate per site (assumes all nucleotides equally probable): {error_rate}">\n' \
        #                    f'##FILTER=<ID=num_samples,Description="Number of samples (VCF files) in full analysis: {num_samples}">\n' \
        #                    f'##FILTER=<ID=FP_cutoff,Description="Analysis-wide false discovery rate (FP) cutoff: {FP_cutoff}">\n' \
        #                    f'##FILTER=<ID=min_DP,Description="Read depth (coverage) cutoff (min allowed): {min_DP}">\n' \
        #                    f'##FILTER=<ID=min_AC,Description="Allele count cutoff (min reads allowed to support minor allele): {min_AC}">\n' \
        #                    f'##FILTER=<ID=min_AF,Description="Minor allele frequency cutoff (min allowed): {min_AF}">'

        # new_INFO_lines = f'##INFO=<ID=MULTIALLELIC,Number=0,Type=Flag,Description="Indicates whether a site is multiallelic (more than one ALT allele)\">'
        # new_metadata_lines = f'##FILTER=<ID=error_rate,Description="Sequencing error rate per site (assumes all nucleotides equally probable): {error_rate}">\n' \
        #                      f'##FILTER=<ID=seq_len,Description="Reference sequence genome length (nucleotides): {seq_len}">\n' \
        #                      f'##FILTER=<ID=num_samples,Description="Number of samples (VCF files) in full analysis: {num_samples}">\n' \
        #                      f'##FILTER=<ID=FP_cutoff,Description="Analysis-wide false-positive (FP) count cutoff: {FP_cutoff}">\n' \
        #                      f'##FILTER=<ID=min_DP,Description="Read depth (coverage) cutoff (min allowed): {min_DP}">\n' \
        #                      f'##FILTER=<ID=min_AC,Description="Allele count cutoff (min reads allowed to support allele): {min_AC}">\n' \
        #                      f'##FILTER=<ID=min_AF,Description="Minor allele frequency cutoff (min allowed): {min_AF}">\n' \
        #                      f'##FILTER=<ID=max_AF,Description="Major allele frequency cutoff (min allowed): {max_AF}">\n'

        new_metadata_lines = f'##FILTER=<ID=error_rate,Description="VCFgenie, sequencing error rate for single nucleotide substitutions per base sequenced: {error_rate}">\n' \
                             f'##FILTER=<ID=p_cutoff,Description="VCFgenie, p-value cutoff: {p_cutoff}">\n' \
                             f'##FILTER=<ID=min_DP,Description="VCFgenie, read depth (coverage) cutoff (min allowed): {min_DP}">\n' \
                             f'##FILTER=<ID=min_AC,Description="VCFgenie, allele count cutoff (min reads allowed to support allele): {min_AC}">\n' \
                             f'##FILTER=<ID=min_AF,Description="VCFgenie, minor allele frequency cutoff (min allowed): {min_AF}">\n' \
                             f'##FILTER=<ID=max_AF,Description="VCFgenie, major allele frequency cutoff (min allowed): {max_AF}">\n'

        for new_FILTER_line in new_FILTER_lines:
            new_metadata_lines += f'{new_FILTER_line}\n'

        new_metadata_lines += '##INFO=<ID=DECISION,Number=.,Type=String,Description="VCFgenie, decision(s) made to yield STATUS of PASS and/or FAIL">\n' \
                              '##INFO=<ID=STATUS,Number=A,Type=String,Description="VCFgenie, whether the ALT allele(s) is PASS or FAIL">\n' \
                              '##INFO=<ID=MULTIALLELIC,Number=0,Type=Flag,Description="VCFgenie, site is multiallelic (more than one ALT allele)">\n' \
                              '##INFO=<ID=REF_FAIL,Number=0,Type=Flag,Description="VCFgenie, REF allele failed to meet the criteria">\n' \
                              f'##INFO=<ID={PVR_key_new},Number=1,Type=Float,Description="VCFgenie, p-value for REF allele">\n'\
                              f'##INFO=<ID={PVA_key_new},Number=A,Type=Float,Description="VCFgenie, p-value for ALT allele(s)">\n' \
                              f'##FORMAT=<ID={PVR_key_new},Number=1,Type=Float,Description="VCFgenie, p-value for REF allele">\n' \
                              f'##FORMAT=<ID={PVA_key_new},Number=A,Type=Float,Description="VCFgenie, p-value for ALT allele(s)">'
        # Others may be written if encountered just below for AC_key_new, AF_key_new, or DP_key_new

        # TODO: here, add p_value
        # Questions:
        # 1) does p-value get written to BOTH format AND INFO?
        # 2) should it be REF,ALTs or separate for REF and then ALTs?

        # Keep track of numbers of records/samples/pass
        this_record_n = 0
        # this_variant_n = 0
        # this_variant_pass_n = 0

        this_ALLELE_n = 0
        this_ALLELE_pass_n = 0
        this_ALLELE_fail_n = 0

        this_REF_n = 0
        this_REF_pass_n = 0
        this_REF_fail_n = 0

        this_ALT_n = 0
        this_ALT_pass_n = 0
        this_ALT_fail_n = 0

        this_MAJOR_n = 0
        this_MAJOR_pass_n = 0
        this_MAJOR_fail_n = 0

        this_MINOR_n = 0
        this_MINOR_pass_n = 0
        this_MINOR_fail_n = 0

        sample_n += 1
        this_sample = ''

        INFO_metadata_dd: Dict[str, dict] = {}
        FORMAT_metadata_dd: Dict[str, dict] = {}

        seen_AC_key_new_FLAG = False
        new_AC_FORMAT_line = f'##FORMAT=<ID={AC_key_new},Number=A,Type=Integer,Description="VCFgenie, ALT allele count after processing">'
        seen_AF_key_new_FLAG = False
        new_AF_FORMAT_line = f'##FORMAT=<ID={AF_key_new},Number=A,Type=Float,Description="VCFgenie, ALT allele frequency after processing">'

        for line in this_VCF_fh:
            # chomp newline
            line = line.rstrip()

            if line.startswith("##"):
                # SAVE metadata in the INFO_metadata_dd or FORMAT_metadata_dd
                # r'^##(\w+)=<ID=(\w+),Number=(\w+),Type=(\w+),Description="(.+)"$'
                # r'^##( 0 )=<ID=( 1 ),Number=( 2 ),Type=( 3 ),Description="( 4)"$'
                if match := re_VCF_INFO_metadata.match(line):
                    if line.startswith("##INFO="):
                        # print(f'INFO_line={line}')
                        grps = match.groups()  # 0-based for groups BUT 1-based if using .sub - hate Python
                        ID = grps[1]
                        INFO_metadata_dd[ID] = {}
                        INFO_metadata_dd[ID]['Number'] = grps[2]
                        INFO_metadata_dd[ID]['Type'] = grps[3]
                        INFO_metadata_dd[ID]['Description'] = grps[4]

                        this_out_fh.write(line + "\n")
                    elif line.startswith("##FORMAT="):
                        # print(f'FORMAT_line={line}')
                        grps = match.groups()  # 0 is whole match; groups are 1-based
                        ID = grps[1]
                        FORMAT_metadata_dd[ID] = {}
                        FORMAT_metadata_dd[ID]['Number'] = grps[2]
                        FORMAT_metadata_dd[ID]['Type'] = grps[3]
                        FORMAT_metadata_dd[ID]['Description'] = grps[4]

                        if ID == AC_key_new:
                            print(f'### WARNING: FORMAT key="{ID}" already exists and will be OVERWRITTEN')
                            # WRITE *new* definition
                            this_out_fh.write(f'{new_AC_FORMAT_line}\n')
                        elif ID == AF_key_new:
                            print(f'### WARNING: FORMAT key="{ID}" already exists and will be OVERWRITTEN')
                            # WRITE *new* definition
                            this_out_fh.write(f'{new_AF_FORMAT_line}\n')
                        else:
                            # WRITE existing metadata headers
                            this_out_fh.write(line + "\n")

                else:
                    # WRITE existing metadata headers
                    this_out_fh.write(line + "\n")

            elif line.startswith("#"):
                # WRITE new AC, AF FORMAT metadata if they were not already seen and replaced
                if not seen_AC_key_new_FLAG:
                    this_out_fh.write(f'{new_AC_FORMAT_line}\n')

                if not seen_AF_key_new_FLAG:
                    this_out_fh.write(f'{new_AF_FORMAT_line}\n')

                # WRITE new metadata
                this_out_fh.write(new_metadata_lines + "\n")

                # PRINT gathered metadata information
                # print('INFO_metadata_dd:')
                # pprint(INFO_metadata_dd)
                # print(f'FORMAT_metadata_dd:')
                # pprint(FORMAT_metadata_dd)

                # TEST that the keys given at top of script are present in FORMAT
                needed_keys_set = {AC_key, AF_key, DP_key}
                FORMAT_keys_set = set(FORMAT_metadata_dd.keys())
                if not needed_keys_set.issubset(FORMAT_keys_set):
                    sys.exit(f'\n### ERROR: one or more keys="{needed_keys_set}" are not present in FORMAT="{FORMAT_keys_set}"')

                # WRITE header line
                this_out_fh.write(line + "\n")

                # RECORD sample name
                this_sample = line.split('\t')[9]
                print(f'Analyzing sample={this_sample}')
            else:
                this_record_n += 1
                total_record_n += 1  # all records examined, no matter how they're categorized

                line_list = line.split("\t")

                if len(line_list) != 10:
                    sys.exit('\n# ERROR: each VCF file must contain the called SNPs for exactly 1\n'
                             '# pooled sequencing sample, and therefore must contain exactly 10 columns:\n'
                             '\t(1) CHROM\n'
                             '\t(2) POS\n'
                             '\t(3) ID\n'
                             '\t(4) REF\n'
                             '\t(5) ALT\n'
                             '\t(6) QUAL\n'
                             '\t(7) FILTER\n'
                             '\t(8) INFO\n'
                             '\t(9) FORMAT [OPTIONAL]\n'
                             '\t(10) <sample_name>\n')

                # -------------------------------------------------------------
                # EXTRACT DATA for this SNP
                this_CHROM = line_list[0]
                this_POS = int(line_list[1])
                this_ID = line_list[2]
                this_REF = str(line_list[3])
                this_ALT_rec = str(line_list[4])
                this_QUAL = line_list[5]
                this_FILTER = line_list[6]  # formatted as a semicolon-separated list, so not possible to CSV by variant
                this_INFO_rec = line_list[7]
                this_FORMAT_rec = line_list[8]
                this_sample_rec = line_list[9]

                # -------------------------------------------------------------
                # PROCESS DATA
                # ==> ALT <==
                # this_ALT_rec_comma_n = this_ALT_rec.count(',')
                this_ALT_list = this_ALT_rec.split(',')
                # this_ALT_n = len(this_ALT_list)  # ALREADY USE this_ALT_n for another purpose

                # ==> FILTER <==
                this_FILTER_new_list: List[str] = []

                # ==> INFO <==
                # re_this_INFOgroups = re_VCF_INFO.match(this_INFO)
                this_INFO_list = this_INFO_rec.split(';')
                this_INFO_list_keys = list(this_INFO_list)
                this_INFO_list_keys = [re_VCF_INFO.sub(r'\1', x) for x in this_INFO_list_keys]

                this_INFO_list_values: List[Optional[str]] = [re_VCF_INFO.sub(r'\2', x) if '=' in x else None for x in list(this_INFO_list)]
                # this_INFO_list_values = list(this_INFO_list)
                # this_INFO_list_values = [re_VCF_INFO.sub(r'\2', x) if '=' in x else None for x in this_INFO_list_values]
                this_INFO_data = dict(zip(this_INFO_list_keys, this_INFO_list_values))
                # this_DP = re_VCF_DP.search(this_INFO_rec)
                # this_DP = this_INFO_rec[this_DP.start():this_DP.end()]  # AF=0.367003
                # this_DP = int(this_DP.replace("DP=", ""))
                # this_INFO_DP = this_INFO_data[DP_key]  # TODO: option later?
                # this_AF_rec = re_VCF_AF.search(this_INFO_rec)
                # this_AF_rec = this_INFO_rec[this_AF_rec.start():this_AF_rec.end()]  # AF=0.367003; AF=0,1
                # this_AF_rec = this_AF_rec.replace('AF=', '')
                # this_INFO_AF = this_INFO_data[AF_key]  # TODO: option later?

                # validate INFO_rules
                if len(INFO_rule_lt) > 0:
                    for i, rule in enumerate(INFO_rule_lt):
                        (key, operator, value) = rule

                        if key not in this_INFO_data.keys():
                            sys.exit(f'\n### ERROR: key="{key}" specified in --INFO_rules is not present in INFO')

                # print(f'this_INFO_rec={this_INFO_rec}')
                # print(f'this_INFO_list_keys={this_INFO_list_keys}')
                # print(f'this_INFO_list_values={this_INFO_list_values}')
                # print(f'this_INFO_data={this_INFO_data}')
                # print(f'this_INFO_DP={this_INFO_DP}')
                # print(f'this_INFO_AF={this_INFO_AF}')

                # ==> FORMAT <==
                this_FORMAT_list = this_FORMAT_rec.split(':')  # original ORDER maintained

                # ==> <sample> <==
                this_sample_list = this_sample_rec.split(':')  # original ORDER maintained
                # print(f'this_sample_rec={this_sample_rec}')

                # this_SB = re_VCF_SB.match(this_INFO) # SEARCH INSTEAD
                # this_DP4_grps = VCF_DP4_grpregex.match(this_INFO) # SEARCH INSTEAD

                if len(this_FORMAT_list) != len(this_sample_list):
                    sys.exit('\n### ERROR: number of data entries in FORMAT and <sample> columns do not match: '
                             f'file="{this_VCF_file}";seq="{this_CHROM}";pos="{this_POS}";REF="{this_REF}";ALT_list="{this_ALT_list}"')

                # SAMPLE DATA: create a dictionary
                this_sample_data = dict(zip(this_FORMAT_list, this_sample_list))

                # validate sample_rules
                if len(sample_rule_lt) > 0:
                    for i, rule in enumerate(sample_rule_lt):
                        (key, operator, value) = rule

                        if key not in this_sample_data.keys():
                            sys.exit(f'\n### ERROR: key="{key}" specified in --sample_rules is not present in FORMAT/sample')

                # print(f'this_sample_data={this_sample_data}')
                # print(f'DP=this_sample_data[{DP_key}]={this_sample_data[DP_key]}')
                # print(f'AF=this_sample_data[{AF_key}]={this_sample_data[AF_key]}')

                # if this_ALT_rec_comma_n > 0:
                if len(this_ALT_list) > 1:  # this_ALT_n > 1:
                    print(f'MULTIALLELIC: file={this_VCF_file}'
                          f';seq={this_CHROM}'
                          f';pos={this_POS}'
                          f';REF={this_REF}'
                          f';ALT={this_ALT_list}')

                # for this_sample_i, this_sample in enumerate(line.samples):  TODO multiple samples
                # this_sample_name = this_sample.sample

                # print()
                # print(f'this_sample_i,this_sample_name ==> {this_sample_i},{this_sample_name} ')

                # -------------------------------------------------------------
                # GATHER DATA using the provided KEYS
                # this_sample_data = this_sample.data._asdict()

                # DP
                this_DP = int(this_sample_data[DP_key])  # just one int

                # IF COVERAGE IS ZERO, EXCLUDE
                if this_DP < 1:
                    continue  # SKIP this iteration; variant effectively filtered out

                # AF
                this_AF_rec = this_sample_data[AF_key]  # list if multiple values, else puts it in a list

                # Ensure AF is in a list, even if just one value: convert to list if not one already
                # N.B.: could also have inferred what to do by the number of ALT alleles, this_ALT_n
                # this_AF_list = [this_AF_rec] if isinstance(this_AF_rec, float) else this_AF_rec

                this_AF_list: List[float] = list(map(float, this_AF_rec.split(',')))
                # this_AF_list = this_AF_rec.split(',')
                # this_AF_list = list(map(float, this_AF_list))

                # AC
                this_AC_rec = this_sample_data[AC_key]  # list if multiple values, else puts it in a list
                this_AC_list: List[float] = list(map(float, this_AC_rec.split(',')))

                # this_AC_list = this_AC_rec.split(',')
                # this_AC_list = list(map(float, this_AC_list))

                # print(f'this_DP={this_DP}')
                # # print(f'type(this_DP)={type(this_DP)}')
                # print(f'this_ALT_list={this_ALT_list}')
                # print(f'this_AF_rec={this_AF_rec}')
                # print(f'this_AF_list={this_AF_list}')
                # print(f'this_AC_list={this_AC_list}')

                if len(this_ALT_list) != len(this_AF_list) or len(this_ALT_list) != len(this_AC_list):
                    sys.exit('\n### ERROR: number of ALT, AF, and AC records do not match: '
                             f'file="{this_VCF_file}";seq="{this_CHROM}";pos="{this_POS}";REF="{this_REF}";ALT_list="{this_ALT_list}"')

                # STORE counts for each ALT allele and CHECK AF is as expected
                allele_AC_dict: Dict[str, int] = defaultdict(int)
                allele_index_dict: Dict[str, Optional[int]] = defaultdict(int)
                for i, this_ALT in enumerate(this_ALT_list):
                    if round(this_DP * this_AF_list[i]) != this_AC_list[i]:
                        # sys.exit(f'\n### ERROR: value of {AC_key} ({round(this_AC_list[i])}) conflicts with allele count implied by {DP_key}*{AF_key} ({round(this_DP * this_AF_list[i])}): '
                        #          f'file="{this_VCF_file}";seq="{this_CHROM}";pos="{this_POS}";REF="{this_REF}";ALT_list="{this_ALT_list}";'
                        #          f'{AC_key}="{this_sample_data[AC_key]}";{DP_key}="{this_DP}";{AF_key}="{this_AF_list[i]}";{AC_key}="{this_AC_list[i]}"')

                        print(f'### WARNING: value of {AC_key} ({round(this_AC_list[i])}) conflicts with allele ' 
                              f'count implied by {DP_key}*{AF_key} ({round(this_DP * this_AF_list[i])}): '
                              f'file="{this_VCF_file}";seq="{this_CHROM}";pos="{this_POS}";REF="{this_REF}";' 
                              f'ALT_list="{this_ALT_list}";{AC_key}="{this_sample_data[AC_key]}";{DP_key}="{this_DP}";' 
                              f'{AF_key}="{this_AF_list[i]}";{AC_key}="{this_AC_list[i]}"; '
                              f'this mismatch may be due to an error inherent to SNP calling; '
                              'or, if may be due to error correction (e.g., flow cell correction in Ion Torrent); '
                              f'VALUE OF {AC_key} = ({round(this_AC_list[i])}) WILL BE USED')

                    allele_AC_dict[str(this_ALT)] = int(this_AC_list[i])

                    if str(this_ALT) in allele_index_dict.keys():  # allele_index_dict[str(this_ALT)] is 0, gets False
                        sys.exit(f'\n### ERROR: multiple ALT alleles for "{this_ALT}": '
                                 f'file="{this_VCF_file}";seq="{this_CHROM}";pos="{this_POS}";REF="{this_REF}";ALT_list="{this_ALT_list}"')

                    allele_index_dict[str(this_ALT)] = i

                # STORE COUNT FOR REF allele, which is the remainder from DP, if none of the ALTs match the REF
                if this_REF in allele_AC_dict.keys():
                    sys.exit(f'\n### ERROR: REF="{AC_key}" also an ALT: '
                             f'file="{this_VCF_file}";seq="{this_CHROM}";pos="{this_POS}";REF="{this_REF}";ALT_list="{this_ALT_list}"')
                else:
                    allele_AC_dict[this_REF] = this_DP - sum(list(allele_AC_dict.values()))
                    allele_index_dict[this_REF] = None  # this could have been anywhere; here for organization
                # print(f'dict(allele_AC_dict)={dict(allele_AC_dict)}')

                # Determine MINOR AND MAJOR alleles
                major_allele = ''
                major_allele_AC = 0  # min possible, to top

                # minor_allele = ''
                # minor_allele_AC = this_DP  # max possible, to bottom

                for allele, AC in allele_AC_dict.items():
                    if AC > major_allele_AC:
                        major_allele = allele
                        major_allele_AC = AC
                    # if AC < minor_allele_AC:
                    #     minor_allele = allele
                    #     minor_allele_AC = AC

                # print(f'major_allele={major_allele};major_allele_AC={major_allele_AC}')
                # print(f'minor_allele={minor_allele};minor_allele_AC={minor_allele_AC}')

                # 'vcf' library' does NOT seem to have a way of treating multiallelics separately, e.g., even in
                # the script /Users/cwnelson88/opt/anaconda3/bin/vcf_filter.py

                # -------------------------------------------------------------
                # LOOP AND TEST ALL ALLELES, includes REF (1) and ALT (1 or more)
                allele_decision_dl: Dict[str, List[str]] = defaultdict(list)
                p_value_dict: Dict[str, float] = defaultdict(float)

                for allele, AC in allele_AC_dict.items():
                    this_ALLELE_n += 1
                    total_ALLELE_n += 1

                    # Retrieve SNV error rate
                    error_rate_SNV = error_rate / 3  # TODO: equal just for now, but later allow 4x4 matrices

                    # AF for reuse in this loop
                    this_AF = AC / this_DP  # NOTE: this is BEFORE correction, which is what we want

                    # Calculate analysis-wide FP for this allele count and coverage
                    # FP_implied = (1 - binom.cdf(AC - 1, this_DP, error_rate_SNV)) * seq_len * num_samples
                    # x - 1 because subtracting everything BEFORE this value (cumulative distribution)

                    # Calculate the variant p-value for this allele count and coverage: P(X>=AC) | ~Binom(DP,error)
                    p_value = 1 - binom.cdf(AC - 1, this_DP, error_rate_SNV)
                    p_value_dict[allele] = p_value

                    # ---------------------------------------------------------
                    # CATEGORIZE
                    decision_list: List[str] = []

                    # failZeroAC: it's ALREADY zero-frequency; could also be fixedALT or fixedREF (below)
                    if AC == 0:
                        decision_list.append('failZeroAC')
                        this_FILTER_new_list.append('AC')

                    # failDP: read depth (coverage) is insufficient at this site
                    if this_DP < min_DP:
                        decision_list.append('failDP')
                        this_FILTER_new_list.append('DP')

                    # failAC: allele count (whether REF or ALT) is insufficient at this site
                    if AC < min_AC:
                        decision_list.append('failAC')
                        this_FILTER_new_list.append('AC')

                    # failMinAF: allele frequency (whether REF or ALT) fails min_AF
                    if AC / this_DP < min_AF:
                        decision_list.append('failMinAF')
                        this_FILTER_new_list.append('AF')

                    # failMaxAF: allele frequency (whether REF or ALT) fails max_AF
                    if AC / this_DP > max_AF:
                        decision_list.append('failMaxAF')
                        this_FILTER_new_list.append('AF')

                    # failINFO: one of the user-provided INFO_rules
                    if len(INFO_rule_lt) > 0:
                        for i, rule in enumerate(INFO_rule_lt):
                            # (key, operator, value) = rule
                            INFO_rule_key: str = rule[0]
                            operator: str = rule[1]
                            value: Any = rule[2]

                            data_Number = str(INFO_metadata_dd[INFO_rule_key]['Number'])
                            data_Type = str(INFO_metadata_dd[INFO_rule_key]['Type'])
                            this_INFO_data_value = this_INFO_data[INFO_rule_key]

                            if data_Number == 'A' and allele == this_REF:  # cannot be applied; move on to next rule
                                continue  # next rule; this Number cannot be applied to REF

                            # Number must be: 1=only one value, same for REF and ALT; A=one per ALT
                            if data_Number in ('1', 'A'):
                                if data_Number == 'A':  # could be a list, in this case
                                    this_INFO_data_list = this_INFO_data_value.split(',')
                                    this_INFO_data_value = this_INFO_data_list[allele_index_dict[allele]]

                                # CAST data type for value: Integer, Float, Flag, Character, or String
                                if data_Type == 'Integer':
                                    value = int(value)
                                    this_INFO_data_value = int(this_INFO_data_value)
                                elif data_Type == 'Float':
                                    value = float(value)
                                    this_INFO_data_value = float(this_INFO_data_value)
                                else:
                                    if operator not in ('==', '!='):
                                        sys.exit(
                                            f'\n### ERROR: operator="{operator}" not compatible with rule={rule} of Type={data_Type}')
                                    else:
                                        value = str(value)
                                        this_INFO_data_value = str(this_INFO_data_value)

                                # if decision_list is None:  # it could have changed during last rule
                                if (operator == '==' and not op.eq(this_INFO_data_value, value)) or \
                                        (operator == '!=' and not op.ne(this_INFO_data_value, value)) or \
                                        (operator == '<' and not op.lt(this_INFO_data_value, value)) or \
                                        (operator == '<=' and not op.le(this_INFO_data_value, value)) or \
                                        (operator == '>=' and not op.ge(this_INFO_data_value, value)) or \
                                        (operator == '>' and not op.gt(this_INFO_data_value, value)):
                                    # print('==> failINFO <==\n')
                                    decision_list.append('failINFO')
                                    this_FILTER_new_list.append(f'{INFO_rule_key}')

                    # failsample: one of the user-provided sample_rules
                    if len(sample_rule_lt) > 0:
                        for i, rule in enumerate(sample_rule_lt):
                            # (sample_rule_key, operator, value) = rule
                            sample_rule_key: str = rule[0]
                            operator: str = rule[1]
                            value: Any = rule[2]

                            data_Number = str(FORMAT_metadata_dd[sample_rule_key]['Number'])
                            data_Type = str(FORMAT_metadata_dd[sample_rule_key]['Type'])
                            this_sample_data_value = this_sample_data[sample_rule_key]

                            if data_Number == 'A' and allele == this_REF:  # cannot be applied; move on to next rule
                                continue  # next rule; this Number cannot be applied to REF

                            # Number must be: 1=only one value, same for REF and ALT; A=one per ALT
                            if data_Number in ('1', 'A'):
                                if data_Number == 'A':  # could be a list, in this case
                                    this_sample_data_list = this_sample_data_value.split(',')
                                    this_sample_data_value = this_sample_data_list[allele_index_dict[allele]]

                                # CAST data type for value: Integer, Float, Flag, Character, or String
                                if data_Type == 'Integer':
                                    value = int(value)
                                    this_sample_data_value = int(this_sample_data_value)
                                elif data_Type == 'Float':
                                    value = float(value)
                                    this_sample_data_value = float(this_sample_data_value)
                                else:
                                    if operator not in ('==', '!='):
                                        sys.exit(f'\n### ERROR: operator="{operator}" not compatible with rule={rule} of Type={data_Type}')
                                    else:
                                        value = str(value)
                                        this_sample_data_value = str(this_sample_data_value)

                                # if decision_list is None:  # it could have changed during last rule
                                if (operator == '==' and not op.eq(this_sample_data_value, value)) or \
                                        (operator == '!=' and not op.ne(this_sample_data_value, value)) or \
                                        (operator == '<' and not op.lt(this_sample_data_value, value)) or \
                                        (operator == '<=' and not op.le(this_sample_data_value, value)) or \
                                        (operator == '>=' and not op.ge(this_sample_data_value, value)) or \
                                        (operator == '>' and not op.gt(this_sample_data_value, value)):
                                    # print('==> failsample <==\n')
                                    decision_list.append('failsample')
                                    this_FILTER_new_list.append(f'{sample_rule_key}')

                                    # # DETERMINE data types for values
                                    # # Ex/ '37' > '5' is *False*, so make sure we're using correct type
                                    # # this_sample_data_value
                                    # if re_non_int.search(this_sample_data_value):
                                    #     if re_non_float.search(this_sample_data_value):
                                    #         this_sample_data_value = str(this_sample_data_value)
                                    #     else:
                                    #         this_sample_data_value = float(this_sample_data_value)
                                    # else:
                                    #     this_sample_data_value = int(this_sample_data_value)
                                    #
                                    # # value
                                    # if re_non_int.search(value):
                                    #     if re_non_float.search(value):
                                    #         value = str(value)
                                    #     else:
                                    #         value = float(value)
                                    # else:
                                    #     value = int(value)

                            else:  # data_Number is not '1' or 'A'
                                print(f'### WARNING: no method for implementing rule={rule}, which has Number={data_Number}; will not be implemented')

                    # # failFP: allele FAILS because insufficient reads support it (FP)
                    # if FP_implied > FP_cutoff:
                    #     decision_list.append('failFP')
                    #     this_FILTER_new_list.append('FP')

                    # failp: allele FAILS the p-value cutoff supplied as --p_cutoff
                    if p_value >= p_cutoff:
                        decision_list.append('failp')
                        this_FILTER_new_list.append('p')

                    # OTHERWISE PASS
                    if len(decision_list) == 0:  # PASS
                        decision_list.append('pass')
                        this_FILTER_new_list.append('PASS')

                    # # DIE if no decision was reached
                    # if decision is None:
                    #     sys.exit(f'\n### ERROR decision="None" at '
                    #              f'file="{this_VCF_file}";seq="{this_CHROM}";pos="{this_POS}";REF="{this_REF}";ALT_list="{this_ALT_list}"')

                    # DIE if nonsensical decision reached
                    if len(decision_list) > 1 and 'pass' in decision_list:
                        sys.exit('\n### ERROR: more than one decision made, but PASS among them')

                    # TALLY pass/fail
                    if decision_list[0] == 'pass':
                        this_ALLELE_pass_n += 1
                        total_ALLELE_pass_n += 1
                        # PASS will be looped in the decision_list, BELOW
                        # ALLELE_n_dict['pass'] += 1
                        # ALLELE_AF_dl['pass'].append(this_AF)  # BEFORE correction
                        # ALLELE_AC_dl['pass'].append(AC)
                        # ALLELE_DP_dl['pass'].append(this_DP)

                    else:
                        this_ALLELE_fail_n += 1
                        total_ALLELE_fail_n += 1
                        ALLELE_n_dict['fail'] += 1
                        ALLELE_AF_dl['fail'].append(this_AF)  # BEFORE correction
                        ALLELE_AC_dl['fail'].append(AC)
                        ALLELE_DP_dl['fail'].append(this_DP)

                    # TALLY pass/fail by MAJOR/MINOR
                    if allele == major_allele:
                        this_MAJOR_n += 1
                        total_MAJOR_n += 1

                        if decision_list[0] == 'pass':
                            this_MAJOR_pass_n += 1
                            total_MAJOR_pass_n += 1

                        else:
                            this_MAJOR_fail_n += 1
                            total_MAJOR_fail_n += 1
                            MAJOR_n_dict['fail'] += 1
                            MAJOR_AF_dl['fail'].append(this_AF)  # BEFORE correction
                            MAJOR_AC_dl['fail'].append(AC)
                            MAJOR_DP_dl['fail'].append(this_DP)
                    else:
                        this_MINOR_n += 1
                        total_MINOR_n += 1

                        if decision_list[0] == 'pass':
                            this_MINOR_pass_n += 1
                            total_MINOR_pass_n += 1

                        else:
                            this_MINOR_fail_n += 1
                            total_MINOR_fail_n += 1
                            MINOR_n_dict['fail'] += 1
                            MINOR_AF_dl['fail'].append(this_AF)  # BEFORE correction
                            MINOR_AC_dl['fail'].append(AC)
                            MINOR_DP_dl['fail'].append(this_DP)

                    # TALLY pass/fail by REF/ALT
                    if allele == this_REF:
                        this_REF_n += 1
                        total_REF_n += 1

                        if decision_list[0] == 'pass':
                            this_REF_pass_n += 1
                            total_REF_pass_n += 1

                        else:
                            this_REF_fail_n += 1
                            total_REF_fail_n += 1
                            REF_n_dict['fail'] += 1
                            REF_AF_dl['fail'].append(this_AF)  # BEFORE correction
                            REF_AC_dl['fail'].append(AC)
                            REF_DP_dl['fail'].append(this_DP)
                            this_INFO_data['REF_FAIL'] = None  # REF only

                    else:
                        this_ALT_n += 1
                        total_ALT_n += 1

                        if decision_list[0] == 'pass':
                            this_ALT_pass_n += 1
                            total_ALT_pass_n += 1

                        else:
                            this_ALT_fail_n += 1
                            total_ALT_fail_n += 1
                            ALT_n_dict['fail'] += 1
                            ALT_AF_dl['fail'].append(this_AF)  # BEFORE correction
                            ALT_AC_dl['fail'].append(AC)
                            ALT_DP_dl['fail'].append(this_DP)

                    # STORE DECISION for use at THIS SITE
                    allele_decision_dl[allele] = decision_list

                    for decision in decision_list:  # this includes 'pass'
                        # STORE DECISION for reporting
                        # ALL
                        ALLELE_n_dict[decision] += 1
                        ALLELE_AF_dl[decision].append(this_AF)  # BEFORE correction
                        ALLELE_AC_dl[decision].append(AC)
                        ALLELE_DP_dl[decision].append(this_DP)

                        # By MAJOR/MINOR status
                        if allele == major_allele:
                            MAJOR_n_dict[decision] += 1
                            MAJOR_AF_dl[decision].append(this_AF)  # BEFORE correction
                            MAJOR_AC_dl[decision].append(AC)
                            MAJOR_DP_dl[decision].append(this_DP)
                        else:
                            MINOR_n_dict[decision] += 1
                            MINOR_AF_dl[decision].append(this_AF)  # BEFORE correction
                            MINOR_AC_dl[decision].append(AC)
                            MINOR_DP_dl[decision].append(this_DP)

                        # By REF/ALT status
                        if allele == this_REF:
                            REF_n_dict[decision] += 1
                            REF_AF_dl[decision].append(this_AF)  # BEFORE correction
                            REF_AC_dl[decision].append(AC)
                            REF_DP_dl[decision].append(this_DP)
                        else:
                            ALT_n_dict[decision] += 1
                            ALT_AF_dl[decision].append(this_AF)  # BEFORE correction
                            ALT_AC_dl[decision].append(AC)
                            ALT_DP_dl[decision].append(this_DP)

                        # REPORT FAILING to STDOUT
                        if decision in decision_choices_FAIL:
                            # print(f'==> {decision} <==\n')

                            # # By REF/ALT status
                            # REForALT = None
                            # if allele == this_REF:
                            #     REForALT = 'REF'
                            # else:
                            #     REForALT = 'ALT'

                            # PRINT output string for log for FAILS
                            out_string = f'FAILED: file="{this_VCF_file}"' \
                                         f';seq="{this_CHROM}"' \
                                         f';pos="{this_POS}"' \
                                         f';allele_FAILED="{allele}"' \
                                         f';REF="{this_REF}"' \
                                         f';ALT="{",".join(this_ALT_list)}"' \
                                         f';AF="{round(this_AF, 5)}"' \
                                         f';AC="{AC}"' \
                                         f';DP="{this_DP}"' \
                                         f';p="{p_value}"' \
                                         f';DECISION="{decision}"' \
                                         ';STATUS="FAIL"'  # f';allele_FAILED="{REForALT}"'  # f';FP="{round(FP_implied, 5)}"' \
                            print(out_string)

                    # if round(this_AF, 3) == 0.414:
                    #     sys.exit(f'\nfile="{this_VCF_file}"'
                    #              f'\nseq="{this_CHROM}"'
                    #              f'\npos="{this_POS}"'
                    #              f'\nallele_AC_dict={dict(allele_AC_dict)}')

                # PERFORM AC, AF CORRECTION given the decisions made
                # fail_n = len([x for x in allele_decision_dl.values() if x in decision_choices_FAIL])

                # alleles_fail = [allele for allele, decision in allele_decision_dl.items() if decision in decision_choices_FAIL]
                alleles_fail = [allele for allele, decision_list in allele_decision_dl.items() if decision_list[0] != 'pass']
                # fail_n = len(alleles_fail)

                # alleles_pass = [allele for allele, decision in allele_decision_dl.items() if not decision in decision_choices_FAIL]
                alleles_pass = [allele for allele, decision_list in allele_decision_dl.items() if decision_list[0] == 'pass']
                pass_n = len(alleles_pass)

                # SAVE a PRE-CORRECTED copy of allele counts
                allele_AC_dict_preCorr = dict(allele_AC_dict)

                # MAJOR fails but ONE OR MORE MINORS do NOT
                if major_allele in alleles_fail and pass_n > 0:
                    # sys.exit(f'\n### ERROR: major_allele={major_allele} failed but 1 or more minor alleles did not; unexpected behavior: '
                    #          f'file="{this_VCF_file}";seq="{this_CHROM}";pos="{this_POS}";REF="{this_REF}";ALT_list="{this_ALT_list}"')

                    # PRINT output string for log for NOCALL
                    out_string = f'NOCALL: major_allele="{major_allele}" failed but minor allele(s) did not: ' \
                                 f'file="{this_VCF_file}"' \
                                 f';seq="{this_CHROM}"' \
                                 f';pos="{this_POS}"' \
                                 f';REF="{this_REF}"' \
                                 f';ALT="{this_ALT_list}"' \
                                 f';DP="{this_DP}"' \
                                 ';FILTER="NOCALL"'
                    print(out_string)

                    # REMOVE any PASS that exists, APPEND the NOCALL
                    this_FILTER_new_list = [x for x in this_FILTER_new_list if x != 'PASS']
                    this_FILTER_new_list.append('NOCALL')

                else:  # otherwise, employ CORRECTION
                    fail_AC_sum = sum([AC for allele, AC in allele_AC_dict.items() if allele in alleles_fail])
                    pass_AC_sum = this_DP - fail_AC_sum  # which is like the 'passing DP'; if ALL fail, assigned to major

                    # print(f'alleles_fail="{alleles_fail}"')
                    # print(f'alleles_pass="{alleles_pass}"')

                    # REDEFINE the allele counts in allele_AC_dict after assigning by frequency in passing reads
                    # allele_AC_dict = {allele: this_DP * (allele_AC_dict[allele] / pass_AC_sum) for allele in alleles_pass}
                    # print(f'BEFORE CORRECTION, allele_AC_dict={allele_AC_dict}')
                    for this_allele, this_AC in allele_AC_dict.items():
                        if this_allele in alleles_fail:
                            allele_AC_dict[this_allele] = 0
                        elif this_allele in alleles_pass:
                            allele_AC_dict[this_allele] = round(this_DP * (allele_AC_dict[this_allele] / pass_AC_sum))
                        else:
                            sys.exit('\n### ERROR: paradox in which allele neither passes nor fails: '
                                     f'file="{this_VCF_file}";seq="{this_CHROM}";pos="{this_POS}";REF="{this_REF}";ALT_list="{this_ALT_list}"')
                    # print(f'AFTER CORRECTION, allele_AC_dict={allele_AC_dict}')

                # IF ALL alleles fail after correction, assign ALL READS to the MAJOR ALLELE (due to ANY REASON, e.g., failDP)
                # if sum(allele_AC_dict.values()) == 0 and allele_decision_dl[this_REF] == 'failDP' and this_REF == major_allele:
                #     allele_AC_dict[this_REF] = this_DP
                new_AC_sum = sum(allele_AC_dict.values())

                # DIE if the above step didn't take care of it
                if this_DP != new_AC_sum:  # or this_POS == 7529:
                    if new_AC_sum == 0:  # and 'failDP' in allele_decision_dl.values():
                        allele_AC_dict[major_allele] = this_DP
                        # print(f'AFTER CORRECTION 2, allele_AC_dict={allele_AC_dict}')

                    elif this_DP - new_AC_sum == 1:  # if only 1 different, rounding error; give to major allele
                        # print(f'\n### WARNING: rounding error, added 1 to major allele AC')
                        allele_AC_dict[major_allele] += 1

                    elif this_DP - new_AC_sum == -1:  # if only 1 different, rounding error; taken from major allele
                        # print(f'\n### WARNING: rounding error, subtracted 1 from major allele AC')
                        allele_AC_dict[major_allele] -= 1

                    else:
                        sys.exit(f'\n### ERROR: depth differs before ({this_DP}) and after ({sum(allele_AC_dict.values())}) correction for failing alleles: '
                                 f'file="{this_VCF_file}";seq="{this_CHROM}";pos="{this_POS}";REF="{this_REF}";ALT_list="{this_ALT_list}"\n'
                                 f'allele_decision_dl={dict(allele_decision_dl)}\n'
                                 f'allele_AC_dict={dict(allele_AC_dict)}')

                # TODO: make sure that if NEW KEY is identical to the ALREADY KEY, we don't overwrite it too soon

                # -------------------------------------------------------------
                # LOOP ALL ALLELES FOR TRACKING FIXED ONES (REF or ALT)
                for this_allele, this_AC in allele_AC_dict.items():
                    this_AF_preCorr = allele_AC_dict_preCorr[this_allele] / this_DP
                    this_AF_corr = this_AC / this_DP  # AFTER correction
                    if round(this_AF_corr, 3) == 1.000:
                        # print(f'this_REF={this_REF};this_allele={this_allele}')
                        if this_allele == this_REF:
                            ALLELE_n_dict['fixedREF'] += 1
                            ALLELE_AF_dl['fixedREF'].append(this_AF_preCorr)
                            ALLELE_AC_dl['fixedREF'].append(this_AC)
                            ALLELE_DP_dl['fixedREF'].append(this_DP)

                            REF_n_dict['fixedREF'] += 1
                            REF_AF_dl['fixedREF'].append(this_AF_preCorr)
                            REF_AC_dl['fixedREF'].append(this_AC)
                            REF_DP_dl['fixedREF'].append(this_DP)

                            if this_allele == major_allele:
                                MAJOR_n_dict['fixedREF'] += 1
                                MAJOR_AF_dl['fixedREF'].append(this_AF_preCorr)
                                MAJOR_AC_dl['fixedREF'].append(this_AC)
                                MAJOR_DP_dl['fixedREF'].append(this_DP)
                            else:
                                MINOR_n_dict['fixedREF'] += 1
                                MINOR_AF_dl['fixedREF'].append(this_AF_preCorr)
                                MINOR_AC_dl['fixedREF'].append(this_AC)
                                MINOR_DP_dl['fixedREF'].append(this_DP)
                        else:
                            ALLELE_n_dict['fixedALT'] += 1
                            ALLELE_AF_dl['fixedALT'].append(this_AF_preCorr)
                            ALLELE_AC_dl['fixedALT'].append(this_AC)
                            ALLELE_DP_dl['fixedALT'].append(this_DP)

                            ALT_n_dict['fixedALT'] += 1
                            ALT_AF_dl['fixedALT'].append(this_AF_preCorr)
                            ALT_AC_dl['fixedALT'].append(this_AC)
                            ALT_DP_dl['fixedALT'].append(this_DP)

                            if this_allele == major_allele:
                                MAJOR_n_dict['fixedALT'] += 1
                                MAJOR_AF_dl['fixedALT'].append(this_AF_preCorr)
                                MAJOR_AC_dl['fixedALT'].append(this_AC)
                                MAJOR_DP_dl['fixedALT'].append(this_DP)
                            else:
                                MINOR_n_dict['fixedALT'] += 1
                                MINOR_AF_dl['fixedALT'].append(this_AF_preCorr)
                                MINOR_AC_dl['fixedALT'].append(this_AC)
                                MINOR_DP_dl['fixedALT'].append(this_DP)

                # -------------------------------------------------------------
                # LOOP the ALT alleles (usually just one)
                # REF_fails = allele_AC_dict[this_REF] == 0  # True or False
                decision_list_ALL: List[str] = []  # 'pass' included
                status_list: List[str] = []  # 'pass' included
                new_AC_list: List[int] = []
                new_AF_list: List[float] = []
                PVA_list: List[float] = []

                for i, this_ALT in enumerate(this_ALT_list):  # ordered TODO <== do I mean the ORIGINAL order? I think so
                    # this_variant_n += 1
                    # total_ALT_n += 1

                    # Determine if the ALT is MAJOR
                    # this_AF = float(this_AF_list[i])
                    this_AF_corr = allele_AC_dict[this_ALT] / this_DP
                    ALT_IS_MAJOR = this_AF_corr >= 0.5  # True or False

                    # if ALT_IS_MAJOR and round(this_AF, 3) == 1.000:
                    #     ALLELE_n_dict['fixedALT'] += 1
                    #     ALT_n_dict['fixedALT'] += 1

                    # Check it matches the previously inferred major allele (technically possible no allele AF >= 0.5)
                    if ALT_IS_MAJOR and this_ALT != major_allele:
                        sys.exit(f'### ERROR: major allele dubious: this_ALT="{this_ALT}"; major_allele="{major_allele}" at: '
                                 f'file="{this_VCF_file}";seq="{this_CHROM}";pos="{this_POS}";REF="{this_REF}";ALT_list="{this_ALT_list}"')

                    # Construct new records
                    decision_list_ALT = allele_decision_dl[this_ALT]
                    # decision_list_ALL.append(decision)
                    decision_list_ALL.extend(decision_list_ALT)
                    # status = 'FAIL' if decision in decision_choices_FAIL else 'PASS'  # TODO FIX THIS
                    status = 'FAIL' if decision_list_ALT[0] != 'pass' else 'PASS'  # fixed bug from '_ALL'
                    status_list.append(status)
                    new_AC_list.append(allele_AC_dict[this_ALT])
                    new_AF_list.append(round(this_AF_corr, 7))
                    PVA_list.append(p_value_dict[this_ALT])

                # -------------------------------------------------------------
                # REPLACE or ADD DATA for SAMPLE

                # ADD new keys if not present
                if AC_key_new not in this_FORMAT_list:
                    this_FORMAT_list.append(AC_key_new)
                if AF_key_new not in this_FORMAT_list:
                    this_FORMAT_list.append(AF_key_new)
                if PVR_key_new not in this_FORMAT_list:
                    this_FORMAT_list.append(PVR_key_new)
                if PVA_key_new not in this_FORMAT_list:
                    this_FORMAT_list.append(PVA_key_new)

                # REPLACE data
                this_sample_data[AC_key_new] = ','.join(map(str, new_AC_list))  # will usually just be one (no comma)
                this_sample_data[AF_key_new] = ','.join(map(str, new_AF_list))  # will usually just be one (no comma)
                this_sample_data[PVR_key_new] = p_value_dict[this_REF]
                this_sample_data[PVA_key_new] = ','.join(map(str, PVA_list))  # will usually just be one (no comma)

                # -------------------------------------------------------------
                # WRITE THE MODIFIED RECORD

                # INITIALIZE new line list
                new_line_list = [this_CHROM, this_POS, this_ID, this_REF, this_ALT_rec, this_QUAL]  # , this_FILTER]

                # ASSEMBLE FILTER record
                # Carry over any FILTER flags that were NOT 'PASS'
                this_FILTER_set = set(this_FILTER.split(';'))
                this_FILTER_set.difference_update({'PASS'})

                # set because, e.g., if it's all PASS, only one PASS necessary
                this_FILTER_new_set = set(this_FILTER_new_list)
                # if len(this_FILTER_new_set) == 1 and list(this_FILTER_new_set)[0] == 'PASS':
                #     this_FILTER_new_list = ['PASS']

                # Combine previous non-PASS FILTER flags with the new ones
                this_FILTER_new_set.update(this_FILTER_set)

                # convert to sorted (auto-generates list)
                this_FILTER_new_list = sorted(this_FILTER_new_set)

                # combine and append
                new_FILTER_record = ';'.join(this_FILTER_new_list)
                new_line_list.append(new_FILTER_record)

                # ASSEMBLE INFO record
                # recall: this_INFO_data = dict(zip(this_INFO_list_keys, this_INFO_list_values))

                # STRIPPED VERSION of keys for debugging
                # this_INFO_list_keys = ['AF', 'AO', 'DP', 'FAO', 'FDP', 'FRO', 'FSAF', 'FSAR', 'FSRF', 'FSRR', 'REF_FAIL']  # 'FR',

                # -------------------------------------------------------------
                # OVERWRITE INFO DATA WITH SAMPLE DATA if requested by user
                if overwrite_INFO:
                    for INFO_key in this_INFO_data.keys():
                        # print(f'INFO_key="{INFO_key}" ', end='')
                        if INFO_key in this_sample_data:
                            # print(f'INFO={this_INFO_data[INFO_key]} sample={this_sample_data[INFO_key]} | ', end='')
                            this_INFO_data[INFO_key] = this_sample_data[INFO_key]
                            # print(f'NEW: INFO={this_INFO_data[INFO_key]} sample={this_sample_data[INFO_key]}')
                        # else:
                        #     print()

                # UPDATE key list in case any were added, preserving order of the first ones
                this_INFO_list_keys.extend(set(this_INFO_data.keys()).difference(this_INFO_list_keys))

                new_INFO_list: List[str] = []
                for key in this_INFO_list_keys:
                    new_INFO_list.append(key if this_INFO_data[key] is None else f'{key}={this_INFO_data[key]}')

                # Add entries for DECISION, STATUS, and MULTIALLELIC
                new_INFO_list.append(f'DECISION={",".join(map(str, decision_list_ALL))}')
                new_INFO_list.append(f'STATUS={",".join(map(str, status_list))}')

                if this_ALT_n > 1:
                    new_INFO_list.append('MULTIALLELIC')

                # Add p-values
                new_INFO_list.append(f'{PVR_key_new}={p_value_dict[this_REF]}')
                new_INFO_list.append(f'{PVA_key_new}={",".join(map(str, PVA_list))}')

                # append INFO to the line
                new_INFO_record = ';'.join(new_INFO_list)
                new_line_list.append(new_INFO_record)

                # ASSEMBLE FORMAT record
                new_FORMAT_record = ':'.join(this_FORMAT_list)
                new_line_list.append(new_FORMAT_record)

                # ASSEMBLE SAMPLE record
                # recall: this_sample_data = dict(zip(this_FORMAT_list, this_sample_list))
                new_sample_list: List[str] = []
                for key in this_FORMAT_list:
                    new_sample_list.append(f'{this_sample_data[key]}')  # NO key for the sample bit
                new_sample_record = ':'.join(new_sample_list)
                new_line_list.append(new_sample_record)

                # PRINT NEW LINE
                this_out_fh.write('\t'.join(map(str, new_line_list)) + '\n')

            # FINISH PROCESSING and WRITING line

        # FINISH PROCESSING VCF file
        this_out_fh.close()

        # VERIFY NEW VCF file is RECOGNIZABLE to pyvcf module
        try:
            this_out_vcf_Reader = vcf.Reader(filename=this_out_file)  # opens without error if exists, even if not VCF
        except AttributeError:
            sys.exit(f'\n### ERROR: new VCF_file="{this_out_file}" not readable by pyvcf. VCF file wrong format?\n')

        try:
            next(this_out_vcf_Reader)
            # first_rec = next(this_out_vcf_Reader)
            # print(f'first_rec.INFO={first_rec.INFO}')
        except ValueError:
            sys.exit(f'\n### ERROR: records in new VCF_file="{this_out_file}" not readable by pyvcf. '
                     'VCF file wrong format?\n')

        # Print passing variants for this VCF
        print()
        print(f'FILE SUMMARY for file={this_VCF_file}:')
        print(f'{this_VCF_file}:record_n={this_record_n}')
        print(f'{this_VCF_file}:ALLELE_n={this_ALLELE_n}')
        print(f'{this_VCF_file}:REF_n={this_REF_n}')
        print(f'{this_VCF_file}:ALT_n={this_ALT_n}')
        # print(f'Found {this_ALT_n} ALT alleles in {this_record_n} records (rows)')

        print(f'{this_VCF_file}:ALLELES_pass={this_ALLELE_pass_n}/{this_ALLELE_n}={round(100 * this_ALLELE_pass_n / this_ALLELE_n, 1)}%')
        print(f'{this_VCF_file}:REF_pass={this_REF_pass_n}/{this_REF_n}={round(100 * this_REF_pass_n / this_REF_n, 1)}%')
        print(f'{this_VCF_file}:ALT_pass={this_ALT_pass_n}/{this_ALT_n}={round(100 * this_ALT_pass_n / this_ALT_n, 1)}%')
        print()

    # -------------------------------------------------------------------------
    # Debugging log
    # PAP288913_F1.tvc_no_pad.vcf

    # VCFgenie.out - each line is an ALLELE
    # allele_FAILED=REF ==> 26 (at most one per line)
    # allele_FAILED=ALT ==> 9 (at most one per line)
    # STATUS=FAIL ==> 35
    # 'FAILED: ' ==> 35

    # PAP288913_F1.tvc_no_pad_filtered.vcf - each line is a SITE
    # REF_FAIL ==> 27 (at most once per line) XX WHOLE FILE INCLUDING METADATA
    # [=,]FAIL ==> 10 (could be more than one per line) XX WHOLE FILE INCLUDING METADATA

    # but 1 or more minor alleles did not ==> 26 (at most 1 per line)
    # [=,]FAIL ==> 9 (could be more than 1 per line)
    # SUM is 35  - Q.E.D.

    # num records (lines) is 297

    # -------------------------------------------------------------------------
    # FINISH AND PRINT SUMMARY STATISTICS
    print()
    print('# -----------------------------------------------------------------------------')
    print('============================> SUMMARY STATISTICS <=============================')
    print('# -----------------------------------------------------------------------------')
    print(f'Total samples (files) examined: {sample_n}')
    print(f'Total records (lines) examined: {total_record_n}')
    print(f'Total ALLELES examined: {total_ALLELE_n} ({total_REF_n} REF + {total_ALT_n} ALT; {total_MAJOR_n} MAJOR + {total_MINOR_n} MINOR)')
    print(f'Total ALLELES pass: {ALLELE_n_dict["pass"]}/{total_ALLELE_n}={round(100 * total_ALLELE_pass_n / total_ALLELE_n, 1)}%')
    print(f'Total REF pass: {REF_n_dict["pass"]}/{total_REF_n}={round(100 * total_REF_pass_n / total_REF_n, 1)}%')
    print(f'Total ALT pass: {ALT_n_dict["pass"]}/{total_ALT_n}={round(100 * total_ALT_pass_n / total_ALT_n, 1)}%')
    print(f'Total MAJOR pass: {MAJOR_n_dict["pass"]}/{total_MAJOR_n}={round(100 * total_MAJOR_pass_n / total_MAJOR_n, 1)}%')
    print(f'Total MINOR pass: {MINOR_n_dict["pass"]}/{total_MINOR_n}={round(100 * total_MINOR_pass_n / total_MINOR_n, 1)}%')
    print()
    # print(f'Total variants fail: {ALT_n_dict["fail"]}')
    # print(f'Fraction variants fail: {round(100 * ALT_n_dict["fail"] / total_ALT_n, 1)}%')

    print('=============> ALL ALLELES (at least 2 per record, REF and ALT) <==============')
    print(f'ALLELE_n={total_ALLELE_n}')
    # print(f'ALLELE_fail_n={total_ALLELE_fail_n}={round(100 * total_ALLELE_fail_n / total_ALLELE_n, 1)}%')
    # print(f'ALLELE_pass_n={total_ALLELE_pass_n}={round(100 * total_ALLELE_pass_n / total_ALLELE_n, 1)}%')
    for decision in decision_choices:
        print(f'ALLELE_{decision}={ALLELE_n_dict[decision]}={round(100 * ALLELE_n_dict[decision] / total_ALLELE_n, 1)}%')

    for decision in decision_choices:  # ALLELE_n_dict.keys():
        print('')
        # print('# -----------------------------------------------------------------------------')
        print(f'==> "{decision}": {summaryKey_to_description[decision]}')
        print(f'ALLELE_{decision}_n={ALLELE_n_dict[decision]}')

        if ALLELE_n_dict[decision] > 0:
            print(f'AF_{decision}: {summary_string(ALLELE_AF_dl[decision])}')
            print(f'AC_{decision}: {summary_string(ALLELE_AC_dl[decision])}')
            print(f'DP_{decision}: {summary_string(ALLELE_DP_dl[decision])}')
    print()

    print('=========================> REF ALLELES (1 per record) <========================')
    print(f'REF_n={total_REF_n}')
    # print(f'REF_fail_n={total_REF_fail_n}={round(100 * total_REF_fail_n / total_REF_n, 1)}%')
    # print(f'REF_pass_n={total_REF_pass_n}={round(100 * total_REF_pass_n / total_REF_n, 1)}%')
    for decision in decision_choices:
        if decision not in 'fixedALT':
            print(f'REF_{decision}={REF_n_dict[decision]}={round(100 * REF_n_dict[decision] / total_REF_n, 1)}%')

    for decision in decision_choices:  # REF_n_dict.keys():
        if decision not in 'fixedALT':
            print('')
            # print('# -----------------------------------------------------------------------------')
            print(f'==> "{decision}": {summaryKey_to_description[decision]}')
            print(f'REF_{decision}_n={REF_n_dict[decision]}')

            if REF_n_dict[decision] > 0:
                print(f'AF_{decision}: {summary_string(REF_AF_dl[decision])}')
                print(f'AC_{decision}: {summary_string(REF_AC_dl[decision])}')
                print(f'DP_{decision}: {summary_string(REF_DP_dl[decision])}')

    print()

    print('=====================> ALT ALLELES (1 or more per record) <====================')
    print(f'ALT_n={total_ALT_n}')
    # print(f'ALT_fail_n={total_ALT_fail_n}={round(100 * total_ALT_fail_n / total_ALT_n, 1)}%')
    # print(f'ALT_pass_n={total_ALT_pass_n}={round(100 * total_ALT_pass_n / total_ALT_n, 1)}%')
    for decision in decision_choices:
        if decision not in 'fixedREF':
            print(f'ALT_{decision}={ALT_n_dict[decision]}={round(100 * ALT_n_dict[decision] / total_ALT_n, 1)}%')

    for decision in decision_choices:  # ALT_n_dict.keys():
        if decision not in 'fixedREF':
            print('')
            # print('# -----------------------------------------------------------------------------')
            print(f'==> "{decision}": {summaryKey_to_description[decision]}')
            print(f'ALT_{decision}_n={ALT_n_dict[decision]}')

            if ALT_n_dict[decision] > 0:
                print(f'AF_{decision}: {summary_string(ALT_AF_dl[decision])}')
                print(f'AC_{decision}: {summary_string(ALT_AC_dl[decision])}')
                print(f'DP_{decision}: {summary_string(ALT_DP_dl[decision])}')

    print()

    print('=====================> MAJOR ALLELES (1 or more per record) <====================')
    print(f'MAJOR_n={total_MAJOR_n}')
    # print(f'MAJOR_fail_n={total_MAJOR_fail_n}={round(100 * total_MAJOR_fail_n / total_MAJOR_n, 1)}%')
    # print(f'MAJOR_pass_n={total_MAJOR_pass_n}={round(100 * total_MAJOR_pass_n / total_MAJOR_n, 1)}%')
    for decision in decision_choices:
        print(f'MAJOR_{decision}={MAJOR_n_dict[decision]}={round(100 * MAJOR_n_dict[decision] / total_MAJOR_n, 1)}%')

    for decision in decision_choices:  # MAJOR_n_dict.keys():
        print('')
        # print('# -----------------------------------------------------------------------------')
        print(f'==> "{decision}": {summaryKey_to_description[decision]}')
        print(f'MAJOR_{decision}_n={MAJOR_n_dict[decision]}')

        if MAJOR_n_dict[decision] > 0:
            print(f'AF_{decision}: {summary_string(MAJOR_AF_dl[decision])}')
            print(f'AC_{decision}: {summary_string(MAJOR_AC_dl[decision])}')
            print(f'DP_{decision}: {summary_string(MAJOR_DP_dl[decision])}')

    print()

    print('=====================> MINOR ALLELES (1 or more per record) <====================')
    print(f'MINOR_n={total_MINOR_n}')
    # print(f'MINOR_fail_n={total_MINOR_fail_n}={round(100 * total_MINOR_fail_n / total_MINOR_n, 1)}%')
    # print(f'MINOR_pass_n={total_MINOR_pass_n}={round(100 * total_MINOR_pass_n / total_MINOR_n, 1)}%')
    for decision in decision_choices:
        print(f'MINOR_{decision}={MINOR_n_dict[decision]}={round(100 * MINOR_n_dict[decision] / total_MINOR_n, 1)}%')

    for decision in decision_choices:  # MINOR_n_dict.keys():
        print('')
        # print('# -----------------------------------------------------------------------------')
        print(f'==> "{decision}": {summaryKey_to_description[decision]}')
        print(f'MINOR_{decision}_n={MINOR_n_dict[decision]}')

        if MINOR_n_dict[decision] > 0:
            print(f'AF_{decision}: {summary_string(MINOR_AF_dl[decision])}')
            print(f'AC_{decision}: {summary_string(MINOR_AC_dl[decision])}')
            print(f'DP_{decision}: {summary_string(MINOR_DP_dl[decision])}')
    print()
    print('# -----------------------------------------------------------------------------')

    # -------------------------------------------------------------------------
    # DONE message
    print('\n# -----------------------------------------------------------------------------')
    print('DONE')


# -----------------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
def summary_string(num_list: list, decimals: int = 3, suffix: str = "") -> str:
    """ Calculate summary statistics on an input list of numbers """
    n_num = len(num_list)

    # All numbers?
    all_numeric = all(str(x).replace('-', '').replace('.', '').isnumeric() for x in num_list)

    if n_num == 0 or not all_numeric:
        min_num = 'NA'
        Q1_num = 'NA'
        mean_num = 'NA'
        std_num = 'NA'
        median_num = 'NA'
        Q3_num = 'NA'
        max_num = 'NA'
    else:  # all work if only 1 element
        min_num = np.round(np.min(num_list), decimals)
        Q1_num = np.round(np.quantile(num_list, 0.25), decimals)
        mean_num = np.round(np.mean(num_list), decimals)
        std_num = np.round(np.std(num_list), decimals)
        median_num = np.round(np.median(num_list), decimals)
        Q3_num = np.round(np.quantile(num_list, 0.75), decimals)
        max_num = np.round(np.max(num_list), decimals)

    return f'n={n_num};min={min_num}{suffix};Q1={Q1_num}{suffix};mean={mean_num}{suffix};std={std_num}{suffix}' \
           f';median={median_num}{suffix};Q3={Q3_num}{suffix};max={max_num}{suffix}'


# -----------------------------------------------------------------------------
def test_summary_string() -> None:
    """ Test summary_string() """
    assert summary_string([]) == 'n=0;min=NA;Q1=NA;mean=NA;std=NA;median=NA;Q3=NA;max=NA'
    assert summary_string(['']) == 'n=1;min=NA;Q1=NA;mean=NA;std=NA;median=NA;Q3=NA;max=NA'
    assert summary_string(['Hello']) == 'n=1;min=NA;Q1=NA;mean=NA;std=NA;median=NA;Q3=NA;max=NA'
    assert summary_string([1]) == 'n=1;min=1;Q1=1.0;mean=1.0;std=0.0;median=1.0;Q3=1.0;max=1'
    assert summary_string([1, 5.6, 298.3, 55, 9]) == 'n=5;min=1.0;Q1=5.6;mean=73.78;std=113.933;median=9.0;Q3=55.0;max=298.3'


# -----------------------------------------------------------------------------
# CALL MAIN
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()
