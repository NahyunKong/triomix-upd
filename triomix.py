#################################################################
# Triomix: Trio contamination / UPD analysis pipeline
#
# Original author: Chris Yoon (https://triomix.readthedocs.io/en/latest/)
# UPD mode written by: Nahyun Kong
# Last update Mar 2nd 2026
# contact: nahyun@wustl.edu
#################################################################

#!/usr/bin/env python3
import sys
import os
import subprocess
import shlex
import re
import argparse
import json
import multiprocessing as mp
import pysam
import gzip
import random
import pandas as pd
import numpy as np

VERSION = '0.0.1'


# ============================================================
# 1) ARGUMENT PARSING
# ============================================================
def argument_parser():
    parser = argparse.ArgumentParser(prog='triomix')
    parser.add_argument('--version', action='version', version='%(prog)s v' + VERSION)

    # ---- Inputs ----
    # BAM/CRAM mode
    parser.add_argument('-f', '--father', required=False, help="Father BAM/CRAM")
    parser.add_argument('-m', '--mother', required=False, help="Mother BAM/CRAM")
    parser.add_argument('-c', '--child',  required=False, help="Child BAM/CRAM")

    # Joint VCF mode
    parser.add_argument('--joint_vcf', required=False, default=None,
                        help="Joint multi-sample VCF(.vcf/.vcf.gz/.bcf) with trio in one file.")
    parser.add_argument('--fatherSample', required=False, default=None,
                        help="Father sample name inside --joint_vcf")
    parser.add_argument('--motherSample', required=False, default=None,
                        help="Mother sample name inside --joint_vcf")
    parser.add_argument('--childSample', required=False, default=None,
                        help="Child sample name inside --joint_vcf")

    parser.add_argument('-r', '--reference', required=True, help="Reference FASTA file")
    parser.add_argument('-s', '--snp', required=False, default=None,
                        help="Optional SNP sites BED (or BED.gz). Highly recommended for VCF mode.")
    parser.add_argument('-t', '--thread', required=False, default=1, type=int)
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd())
    parser.add_argument('-p', '--prefix', required=False, default=None)
    parser.add_argument('--runmode', required=False, default='all', choices=['single', 'joint', 'all'])
    parser.add_argument("-u", '--upd', default=1, choices=[0, 1], type=int)
    parser.add_argument('--parent', required=False, action='store_true')
    parser.add_argument('--updMode', required=False, action='store_true')
    parser.add_argument('--dmrMode', required=False, action='store_true')
    parser.add_argument('-d', '--downsample', required=False, default=0.1, type=float)

    args = vars(parser.parse_args())

    # ---- Validate input mode ----
    if args['joint_vcf'] is not None:
        # joint VCF mode
        if not is_vcf(args['joint_vcf']):
            raise SystemExit(f"--joint_vcf must be .vcf/.vcf.gz/.bcf, got: {args['joint_vcf']}")

        missing = [k for k in ('fatherSample', 'motherSample', 'childSample') if not args.get(k)]
        if missing:
            raise SystemExit("--joint_vcf requires trio sample names: missing " +
                             ", ".join("--" + m for m in missing))

        if args['prefix'] is None:
            args['prefix'] = args['childSample']

    else:
        # BAM/CRAM mode
        if not (args.get('father') and args.get('mother') and args.get('child')):
            raise SystemExit("Provide either (--joint_vcf + --fatherSample/--motherSample/--childSample) "
                             "OR all of (-f/-m/-c).")

        for p in (args['father'], args['mother'], args['child']):
            if not is_bam_or_cram(p):
                raise SystemExit(f"BAM/CRAM mode expects .bam/.cram only: {p}")

        if args['prefix'] is None:
            args['prefix'] = sampleNameBam(args['child'])

    return (args['father'], args['mother'], args['child'], args['reference'],
            args['snp'], args['thread'], args['output_dir'], args['prefix'],
            args['runmode'], args['upd'], args['parent'], args['updMode'],
            args['dmrMode'], args['downsample'],
            args['joint_vcf'], args['fatherSample'], args['motherSample'], args['childSample'])

def sampleNameBam(bamFile):
    """get @RG SM: information as sample name from BAM header"""
    bam = pysam.AlignmentFile(bamFile)
    name = bam.header['RG'][0]['SM']
    return name



def is_bam_or_cram(path: str) -> bool:
    p = path.lower()
    return p.endswith(".bam") or p.endswith(".cram")

def is_vcf(path: str) -> bool:
    p = path.lower()
    return p.endswith(".vcf") or p.endswith(".vcf.gz") or p.endswith(".bcf")

def identify_chromosomes(fasta_file):
    """given fasta file, identify the chromosomes that are autosomal"""
    fai_file = fasta_file + '.fai'
    if (os.path.exists(fai_file)):
        with open(fai_file, 'r') as f:
            for line in f:
                chrom, length, *args = line.split('\t')
                if re.search(r'^chr[0-9XY]+$|^[0-9XY]+$', chrom):
                    yield (chrom, float(length))
    else:
        print(f'There is no index file for {fasta_file}. Exiting...')
        sys.exit(1)


def parse_par_bed_file(par_bed):
    """get start and end position of pseudoautosomal regions from a given bed file"""
    par_exclusion_list = []
    with open(par_bed, 'r') as f:
        for line in f:
            chrom, start, end = line.strip().split('\t')
            par_exclusion_list.append((float(start), float(end)))
    return par_exclusion_list


def split_regions(fasta_file, segment_length):
    """splits chromosome into segment lengths"""
    chr_regions = []
    for chrom, chr_length in identify_chromosomes(fasta_file):
        for i in range(0, max(1, int(chr_length / segment_length))):  # fixed here
            start = i * segment_length
            end = (i + 1) * segment_length
            if (chr_length - end <= segment_length):
                end = chr_length
            chr_regions.append(f'{chrom}:{start:.0f}-{end:.0f}')
    return chr_regions


def check_gzip_file(file_path):
    """checks if the file is binary or not"""
    return file_path.endswith(".gz")


def check_region_and_snp_bed(region, snp_bed):
    """if a region does not contain anything on the bed varscan output has no variant, it will cause error downstream when merging"""
    region_chromosome, start_end = region.split(':')
    region_start, region_end = start_end.split('-')
    region_start = float(region_start)
    region_end = float(region_end)
    variant_count = 0

    is_snp_bed_gzip = check_gzip_file(snp_bed)
    f = gzip.open(snp_bed, 'rt') if is_snp_bed_gzip else open(snp_bed, 'r')

    for line in f:
        chrom, start, end, *args = line.strip().split('\t')
        if chrom == region_chromosome:
            if float(start) > region_start and float(start) < region_end:
                variant_count += 1
            if variant_count > 0:
                f.close()
                return True
    f.close()
    return False


def filter_regions_with_snv(region_list, snp_bed):
    """using check_region_and_snp_bed remove regions with no overlaping snp in the region"""
    keep_list = []
    for region in region_list:
        if check_region_and_snp_bed(region, snp_bed):
            keep_list.append(region)
        else:
            print(f'{region} filtered out since it does not have any of the SNPs in the BED file')
    return keep_list


# ============================================================
# 2) BAM/CRAM MODE: mpileup -> counts 
# ============================================================
def mpileup(father_bam, mother_bam, child_bam, region, output_dir, snp_bed, prefix):
    """Run mpileup in a defined region of interest"""
    region_string = re.sub(r':|-', '_', region)

    if snp_bed is not None:
        snp_bed_string = f' -l {snp_bed} '
    else:
        snp_bed_string = ''

    output_file_compressed = os.path.join(output_dir, f'{prefix}_{region_string}.mpileup.gz')

    if not os.path.isfile(output_file_compressed):
        cmd = f'{SAMTOOLS} mpileup -B -Q 20 -q 20 {snp_bed_string}-r {region} -f {REFERENCE} {father_bam} {mother_bam} {child_bam} | gzip -f > {output_file_compressed}'
        os.system(cmd)

    return output_file_compressed


def vaf(alt_count, depth):
    if depth > 0:
        return float(alt_count / depth)
    else:
        return 0


def count_int(count):
    """handle '.' in alt/ref counts"""
    if count is None or count < 0:
        return 0
    else:
        return count


def random_sample_selection(sampling_rate):
    """random sampling function, avoid sampling if rate==1"""
    if sampling_rate == 1:
        return True
    else:
        prob = random.uniform(0, 1)
        if prob < sampling_rate:
            return True
        else:
            return False


def parse_mpileup(mpileup_line, homoref_sampling_rate):
    """Parse mpileup result into counts string"""
    chrom, pos, ref, n_alt_base_father, n_alt_base_mother, n_alt_base_child, \
        father_alt_base, mother_alt_base, child_alt_base, \
        father_depth, mother_depth, child_depth, trio_alt_counts = parse_mpileup_line(mpileup_line)

    if (n_alt_base_father == 1 and n_alt_base_mother == 0):
        alt_base = ''.join(father_alt_base)
        alt_parent = 'F'
        alt_father = trio_alt_counts['father'][alt_base]
        alt_mother = trio_alt_counts['mother'][alt_base]
        child_alt = trio_alt_counts['child'][alt_base]
        father_vaf = vaf(alt_father, father_depth)
        mother_vaf = vaf(alt_mother, mother_depth)

    elif (n_alt_base_father == 0 and n_alt_base_mother == 1):
        alt_base = ''.join(mother_alt_base)
        alt_parent = 'M'
        alt_father = trio_alt_counts['father'][alt_base]
        alt_mother = trio_alt_counts['mother'][alt_base]
        child_alt = trio_alt_counts['child'][alt_base]
        father_vaf = vaf(alt_father, father_depth)
        mother_vaf = vaf(alt_mother, mother_depth)

    elif (n_alt_base_father == 0 and n_alt_base_mother == 0):
        alt_base = 'N'  # parents homoref
        alt_parent = 'NA'
        alt_father = 0
        alt_mother = 0
        # if parents are homoref/homoref, any non ref bases are errors
        child_alt = sum(trio_alt_counts['child'].values())
        father_vaf = vaf(alt_father, father_depth)
        mother_vaf = vaf(alt_mother, mother_depth)

    else:
        father_vaf = -1  # arbitrary to filter out
        mother_vaf = -1  # arbitrary to filter out
        # ambiguous site
        pass

    snvcount = ''
    if (father_vaf > 0.4 and father_vaf < 0.6 and mother_vaf == 0) or (mother_vaf > 0.4 and mother_vaf < 0.6 and father_vaf == 0):
        if father_depth > 10 and mother_depth > 10:
            child_vaf = vaf(child_alt, child_depth)
            snvcount = f'{chrom}\t{pos:.0f}\t{ref}\t{alt_base}\t{child_alt}\t{child_depth}\t{child_vaf}\t{alt_parent}\tNA\t{alt_father}\t{father_depth}\t{father_vaf}\t{alt_mother}\t{mother_depth}\t{mother_vaf}\n'
    elif (father_vaf == 0 and mother_vaf == 1) or (father_vaf == 1 and mother_vaf == 0):
        if father_depth > 10 and mother_depth > 10:
            child_vaf = vaf(child_alt, child_depth)
            snvcount = f'{chrom}\t{pos:.0f}\t{ref}\t{alt_base}\t{child_alt}\t{child_depth}\t{child_vaf}\tNA\t{alt_parent}\t{alt_father}\t{father_depth}\t{father_vaf}\t{alt_mother}\t{mother_depth}\t{mother_vaf}\n'
    elif (father_vaf == 0 and mother_vaf == 0) and random_sample_selection(homoref_sampling_rate):
        if father_depth > 10 and mother_depth > 10:
            child_vaf = vaf(child_alt, child_depth)
            snvcount = f'{chrom}\t{pos:.0f}\t{ref}\t{alt_base}\t{child_alt}\t{child_depth}\t{child_vaf}\tNA\tNA\t{alt_father}\t{father_depth}\t{father_vaf}\t{alt_mother}\t{mother_depth}\t{mother_vaf}\n'
    return snvcount


def parse_mpileup_line(mpileup_line):
    """due to repeated use, decided to make this into a function"""
    split_mpileup = mpileup_line.split('\t')
    chrom = split_mpileup[0]
    pos = float(split_mpileup[1])
    ref = split_mpileup[2]

    trio_alt_counts = dict({'father': None, 'mother': None, 'child': None})
    trio_depth_counts = dict({'father': None, 'mother': None, 'child': None})

    for individual, i in zip(['father', 'mother', 'child'], range(1, 4)):
        mismatch_dict = dict({'A': 0, 'C': 0, 'G': 0, 'T': 0, 'ins': 0, 'del': 0, 'depth': 0})

        bases_index = 3 * i + 1
        depths_index = 3 * i
        depths = int(split_mpileup[depths_index])
        mpiledup = split_mpileup[bases_index].upper()
        insertions = re.findall(r'\+[0-9]+[ACGTNacgtn]+', mpiledup)
        deletions = re.findall(r'-[0-9]+[ACGTNacgtn]+', mpiledup)

        mismatch_dict['ins'] = len(insertions)
        mismatch_dict['del'] = len(deletions)

        mpileupsnv = re.sub(r'\+[0-9]+[ACGTNacgtn]+|-[0-9]+[ACGTNacgtn]+', '', mpiledup)

        mismatch_dict['A'] = mpileupsnv.count('A')
        mismatch_dict['T'] = mpileupsnv.count('T')
        mismatch_dict['G'] = mpileupsnv.count('G')
        mismatch_dict['C'] = mpileupsnv.count('C')
        trio_depth_counts.update({individual: depths})
        trio_alt_counts.update({individual: mismatch_dict})

    father_depth = trio_depth_counts['father']
    mother_depth = trio_depth_counts['mother']
    child_depth = trio_depth_counts['child']

    father_alt_base = [k for k, v in trio_alt_counts['father'].items() if v > 0]
    mother_alt_base = [k for k, v in trio_alt_counts['mother'].items() if v > 0]
    child_alt_base = [k for k, v in trio_alt_counts['child'].items() if v > 0]

    n_alt_base_father = len(father_alt_base)
    n_alt_base_mother = len(mother_alt_base)
    n_alt_base_child = len(child_alt_base)

    return chrom, pos, ref, n_alt_base_father, n_alt_base_mother, n_alt_base_child, father_alt_base, mother_alt_base, child_alt_base, father_depth, mother_depth, child_depth, trio_alt_counts


def parse_mpileup_child_homoalt(mpileup_line):
    """Parse mpileup result into counts string
    Only find places where child is homoalt"""
    chrom, pos, ref, n_alt_base_father, n_alt_base_mother, n_alt_base_child, \
        father_alt_base, mother_alt_base, child_alt_base, \
        father_depth, mother_depth, child_depth, trio_alt_counts = parse_mpileup_line(mpileup_line)

    if (n_alt_base_child == 1):
        alt_base = ''.join(child_alt_base)
        alt_father = trio_alt_counts['father'][alt_base]
        alt_mother = trio_alt_counts['mother'][alt_base]
        child_alt = trio_alt_counts['child'][alt_base]
        father_vaf = vaf(alt_father, father_depth)
        mother_vaf = vaf(alt_mother, mother_depth)
        child_vaf = vaf(child_alt, child_depth)
    else:
        father_vaf = -1
        mother_vaf = -1
        child_vaf = -1
        pass

    snvcount = ''
    if child_vaf == 1:
        if father_depth > 10 and mother_depth > 10:
            snvcount = f'{chrom}\t{pos:.0f}\t{ref}\t{alt_base}\t{child_alt}\t{child_depth}\t{child_vaf}\t{alt_father}\t{father_depth}\t{father_vaf}\t{alt_mother}\t{mother_depth}\t{mother_vaf}\n'

    return snvcount


def get_counts_childhomoalt(mpileup_file):
    """parse mpileup results into a table"""
    output_counts_region = mpileup_file + '.childhomoalt.counts'
    with open(output_counts_region, 'w') as g:
        g.write('chrom\tpos\trefbase\taltbase\talt\tdepth\tvaf\tfather_alt\tfather_depth\tfather_vaf\tmother_alt\tmother_depth\tmother_vaf\n')
        with gzip.open(mpileup_file, 'rt') as f:
            for line in f:
                g.write(parse_mpileup_child_homoalt(line))
        f.close()
    g.close()
    return output_counts_region


def get_child_count(mpileup_file, homoref_sampling_rate):
    """parse mpileup results into a table"""
    output_counts_region = mpileup_file + '.counts'
    with open(output_counts_region, 'w') as g:
        g.write('chrom\tpos\trefbase\taltbase\talt\tdepth\tvaf\thetero_parent\thomoalt_parent\tfather_alt\tfather_depth\tfather_vaf\tmother_alt\tmother_depth\tmother_vaf\n')
        with gzip.open(mpileup_file, 'rt') as f:
            for line in f:
                g.write(parse_mpileup(line, homoref_sampling_rate))
        f.close()
    g.close()
    return output_counts_region


# ============================================================
# 3) VCF MODE HELPERS (joint VCF function)
# ============================================================
def _get_ad_dp(sample):
    """
    Returns (ref_count, alt_count, dp) for a bi-allelic SNP.
    Requires AD and (DP or AD-sum).
    """
    ad = sample.get("AD", None)
    dp = sample.get("DP", None)

    if ad is None:
        return None, None, None
    if len(ad) < 2:
        return None, None, None

    refc = ad[0] if ad[0] is not None else 0
    altc = ad[1] if ad[1] is not None else 0

    if dp is None:
        dp = refc + altc
    return int(refc), int(altc), int(dp)


def _is_simple_snp(rec):
    # only handle bi-allelic SNVs
    if len(rec.alts or []) != 1:
        return False
    if len(rec.ref) != 1:
        return False
    if len(rec.alts[0]) != 1:
        return False
    if rec.ref not in "ACGT":
        return False
    if rec.alts[0] not in "ACGT":
        return False
    return True


def _iter_records_with_optional_bed(vcf, snp_bed):
    """Yield VCF records either by full scan or restricted to BED intervals."""
    if snp_bed is None:
        for rec in vcf.fetch():
            yield rec
    else:
        opener = gzip.open if check_gzip_file(snp_bed) else open
        with opener(snp_bed, "rt") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                chrom, start, end, *rest = line.rstrip("\n").split("\t")
                start = int(start)
                end = int(end)
                for rec in vcf.fetch(chrom, start, end):
                    yield rec

def build_counts_from_joint_vcf(joint_vcf_path, father_sample, mother_sample, child_sample,
                                out_counts_path, snp_bed=None, homoref_sampling_rate=1.0):
    """
    joint multi-sample VCF mode.

    - one VCF/BCF file
    - user provides trio sample names
    - same output columns as mpileup-derived counts
    """
    vcf = pysam.VariantFile(joint_vcf_path)
    header_samples = set(list(vcf.header.samples))

    for s in (father_sample, mother_sample, child_sample):
        if s not in header_samples:
            raise SystemExit(
                f"Sample '{s}' not found in joint VCF header. "
                f"Available samples: {', '.join(sorted(header_samples))}"
            )

    header = (
        "chrom\tpos\trefbase\taltbase\talt\tdepth\tvaf\thetero_parent\thomoalt_parent\t"
        "father_alt\tfather_depth\tfather_vaf\tmother_alt\tmother_depth\tmother_vaf\n"
    )

    with open(out_counts_path, "w") as out:
        out.write(header)

        for rec in _iter_records_with_optional_bed(vcf, snp_bed):
            if not _is_simple_snp(rec):
                continue

            fs = rec.samples[father_sample]
            ms = rec.samples[mother_sample]
            cs = rec.samples[child_sample]

            f_ref, f_alt, f_dp = _get_ad_dp(fs)
            m_ref, m_alt, m_dp = _get_ad_dp(ms)
            c_ref, c_alt, c_dp = _get_ad_dp(cs)

            if None in (f_alt, f_dp, m_alt, m_dp, c_alt, c_dp):
                continue
            if f_dp <= 10 or m_dp <= 10:
                continue

            father_vaf = vaf(f_alt, f_dp)
            mother_vaf = vaf(m_alt, m_dp)
            child_vaf = vaf(c_alt, c_dp)

            hetero_parent = "NA"
            homoalt_parent = "NA"

            if (0.4 < father_vaf < 0.6 and mother_vaf == 0):
                hetero_parent = "F"
            elif (0.4 < mother_vaf < 0.6 and father_vaf == 0):
                hetero_parent = "M"
            elif (father_vaf == 0 and mother_vaf == 1):
                homoalt_parent = "M"
            elif (father_vaf == 1 and mother_vaf == 0):
                homoalt_parent = "F"
            elif (father_vaf == 0 and mother_vaf == 0):
                # match BAM behavior: sample homoref/homoref sites unless snp_bed provided
                keep = True if snp_bed is not None else random_sample_selection(homoref_sampling_rate)
                if not keep:
                    continue
            else:
                continue

            out.write(
                f"{rec.chrom}\t{rec.pos}\t{rec.ref}\t{rec.alts[0]}\t{c_alt}\t{c_dp}\t{child_vaf}\t"
                f"{hetero_parent}\t{homoalt_parent}\t"
                f"{f_alt}\t{f_dp}\t{father_vaf}\t{m_alt}\t{m_dp}\t{mother_vaf}\n"
            )

    print("[JOINT VCF->counts] Head:\n" + "".join(open(out_counts_path).readlines()[:6]), flush=True)
    return out_counts_path


# ============================================================
# 4) RSCRIPTS / SEGMENTATION / PLOTTING
# ============================================================
def run_mle_rscript(count_table, output_dir, runmode, upd):
    """run mle script"""
    cmd = f'{RSCRIPT} {MLE_RSCRIPT} -i {count_table} -o {output_dir} -r {runmode} -u {upd}'
    print(cmd)
    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()
    return 0


def run_mle_parent_rscript(x_combined_counts, output_dir, sexchrom):
    """Run mle script for parent contam by child"""
    cmd = f'{RSCRIPT} {MLE_PARENT_RSCRIPT} -i {x_combined_counts} -o {output_dir} -x {sexchrom}'
    print(cmd)
    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()
    return 0


def run_mle_parent_child_homoalt_rscript(child_homoalt_counts, output_dir):
    """Run mle script for parent contam by child"""
    cmd = f'{RSCRIPT} {MLE_PARENT_CHILD_HOMOALT_RSCRIPT} -i {child_homoalt_counts} -o {output_dir}'
    print(cmd)
    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()
    return 0


def run_plot_rscript(count_table, output_dir, reference, downsample=0.1):
    """run plot script"""
    cmd = f'{RSCRIPT} {PLOT_RSCRIPT} -i {count_table} -o {output_dir} -r {reference} -d {downsample}'
    print(cmd)
    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()
    return 0


def run_plot_parent_rscript(count_table, output_dir, reference, downsample=0.1):
    """run plot script"""
    cmd = f'{RSCRIPT} {PLOT_RSCRIPT_PARENT} -i {count_table} -o {output_dir} -r {reference} -d {downsample}'
    print(cmd)
    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()
    return 0


def run_segmentation(SEGMENTATION_RSCRIPT, count_table, output_dir, segment_length):
    """run segmentation on homo-ref + homo-alt sites for each parental SNPs"""
    cmd = f'{RSCRIPT} {SEGMENTATION_RSCRIPT} -i {count_table} -o {output_dir} -s {segment_length}'
    print(cmd)
    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()
    return 0


def get_paths(path_config):
    """configures the paths to SAMTOOLS AND VARSCAN"""
    with open(path_config) as f:
        path = json.load(f)
    return path['SAMTOOLS'], path['RSCRIPT'], path['GZIP']


# ============================================================
# 5) chrX utilities 
#    NOTE: parse_mpileup_chrX was referenced; adding a working implementation.
# ============================================================
def position_pseudoautosomal(pos, par_list):
    """True if position is pseudoautosomal, False if position is not pseudoautosomal"""
    status = False
    for start, end in par_list:
        if (pos - start) * (pos - end) < 0:
            status = True
    return status


def parse_mpileup_chrX(mpileup_line, sexchrom, par_regions):
    """
    Minimal implementation so get_child_count_chrX works
    Groups:
      - PAR vs non-PAR
      - sexchrom can be "XX", "XY", etc. 
    Output columns expected by get_child_count_chrX:
      chrom pos refbase altbase alt depth vaf chrX_group father_alt father_depth father_vaf mother_alt mother_depth mother_vaf
    """
    chrom, pos, ref, n_alt_base_father, n_alt_base_mother, n_alt_base_child, \
        father_alt_base, mother_alt_base, child_alt_base, \
        father_depth, mother_depth, child_depth, trio_alt_counts = parse_mpileup_line(mpileup_line)

    # chrX_group label
    is_par = position_pseudoautosomal(pos, par_regions) if par_regions is not None else False
    chrX_group = "PAR" if is_par else "nonPAR"

    # choose alt base if exactly one alt in child; else skip
    if n_alt_base_child != 1:
        return ""

    alt_base = ''.join(child_alt_base)
    alt_father = trio_alt_counts['father'].get(alt_base, 0)
    alt_mother = trio_alt_counts['mother'].get(alt_base, 0)
    child_alt = trio_alt_counts['child'].get(alt_base, 0)

    father_vaf = vaf(alt_father, father_depth)
    mother_vaf = vaf(alt_mother, mother_depth)
    child_vaf = vaf(child_alt, child_depth)

    return f'{chrom}\t{pos:.0f}\t{ref}\t{alt_base}\t{child_alt}\t{child_depth}\t{child_vaf}\t{chrX_group}\t{alt_father}\t{father_depth}\t{father_vaf}\t{alt_mother}\t{mother_depth}\t{mother_vaf}\n'


def get_child_count_chrX(mpileup_file, sexchrom, par_regions):
    """parse mpileup results into a table"""
    output_counts_region = mpileup_file + '.chrX.counts'
    with open(output_counts_region, 'w') as g:
        g.write('chrom\tpos\trefbase\taltbase\talt\tdepth\tvaf\tchrX_group\tfather_alt\tfather_depth\tfather_vaf\tmother_alt\tmother_depth\tmother_vaf\n')
        with gzip.open(mpileup_file, 'rt') as f:
            for line in f:
                g.write(parse_mpileup_chrX(line, sexchrom, par_regions))
        f.close()
    g.close()
    return output_counts_region


def x_to_autosome_ratio(counts_file, individual):
    """measures the depth of autosome to chrX to determine the sex of each individual
    individual can be 'father', 'mother', or '' (child)
    """
    df = pd.read_csv(counts_file, sep='\t')
    if individual != 'child':
        individual = individual + '_'
    else:
        individual = ''

    x_depth = np.nanmean(df[df['chrom'].str.contains('^X$|^chrX$') == True][[f'{individual}depth']])
    autosome_depth = np.nanmean(df[df['chrom'].str.contains('^X$|^chrX$') == False][[f'{individual}depth']])
    x_to_auto = float(x_depth / autosome_depth)
    return x_to_auto


def sexchrom_ratio(counts_file):
    child_x_ratio = x_to_autosome_ratio(counts_file, 'child')
    father_x_ratio = x_to_autosome_ratio(counts_file, 'father')
    mother_x_ratio = x_to_autosome_ratio(counts_file, 'mother')
    return child_x_ratio, father_x_ratio, mother_x_ratio


def combine_count_files(file_list, output_file):
    """combines multiple split up count file into one single file"""
    count = 0
    with open(output_file, 'w') as f:
        for afile in file_list:
            with open(afile, 'r') as g:
                if count != 0:
                    next(g)  # header only once
                count += 1
                for line in g:
                    f.write(line)
    return output_file


# ============================================================
# 6) METHYLATION / DMR (UNCHANGED)
# ============================================================
def run_DeepVar(bam, chrom, OUTDIR, name, REFERENCE):
    cmd = f"/opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref={REFERENCE}  --reads={bam} --regions {chrom} --output_vcf={OUTDIR}/{name}/deepvar_{name}.{chrom}.vcf.gz --num_shards=16"
    print(cmd)
    os.system(cmd)
    return (f"{OUTDIR}/{name}/deepvar_{name}.{chrom}.vcf.gz")


def run_Hiphase(bam, vcf, outdir, name, REFERENCE):
    cmd = f"hiphase --bam {bam} --output-bam {outdir}/{name}/hiphase_{name}.bam --vcf {vcf} --output-vcf {outdir}/{name}/hiphase_{name}.vcf.gz --reference {REFERENCE} --threads 32 --ignore-read-groups"
    print(cmd)
    os.system(cmd)
    return (f"{outdir}/{name}/hiphase_{name}.bam")


def run_aligned_bam_to_cpg_scores(bam, outdir, name):
    cmd = f"aligned_bam_to_cpg_scores --bam {bam} --output-prefix {outdir}/{name}/methylation_{name} --model /opt/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite --threads 16"
    print(cmd)
    os.system(cmd)


def run_Methylation(output_dir, bam, UPDbed, name, thread, REFERENCE):
    os.system(f'mkdir -p {output_dir}/{name}')

    # Step 0: Filter bad reads
    print(f"{name}: filter Bam")
    filtered = f"{output_dir}/{name}/{os.path.basename(bam).replace('.bam', '.filtered.bam')}"
    cmd = f"samtools view -b -q 6 {bam} | samtools sort -o {filtered}"
    cmd2 = f"samtools index {filtered}"
    os.system(cmd)
    os.system(cmd2)

    # Step 1: get chromosomes that have UPD
    print(f"{name}: setting up run")
    chromosomes = set()
    with open(UPDbed, 'r') as bed_file:
        next(bed_file)
        for line in bed_file:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) > 0:
                chromosomes.add(fields[0])

    # Step 2: Run deepvar for each chromosome
    print(f"{name}: running deepvar")
    deepVarArgs = []
    for chrom in sorted(chromosomes):
        deepVarArgs.append([filtered, chrom, output_dir, name, REFERENCE])

    with mp.Pool(thread) as pool:
        deepvar_vcfs = pool.starmap(run_DeepVar, deepVarArgs)

    # Step 3: merge deepvar vcfs
    mergedVCF = f"{output_dir}/{name}/deepvar_{name}.merged.vcf"
    with pysam.VariantFile(deepvar_vcfs[0]) as first_vcf:
        with pysam.VariantFile(mergedVCF, 'w', header=first_vcf.header) as out_vcf:
            for record in first_vcf:
                out_vcf.write(record)
            for vcf_path in deepvar_vcfs[1:]:
                with pysam.VariantFile(vcf_path) as vcf:
                    for record in vcf:
                        out_vcf.write(record)

    cmdGzip = f"bgzip {mergedVCF}"
    cmdTabix = f"tabix {mergedVCF}.gz"
    os.system(cmdGzip)
    os.system(cmdTabix)

    # Step 4: hiphase
    print(f"{name}: running hiphase")
    phasebam = run_Hiphase(filtered, f"{mergedVCF}.gz", output_dir, name, REFERENCE)

    os.system(f'cd {output_dir}/{name}')

    # Step 5: methylation
    print(f"{name}: running methylation")
    run_aligned_bam_to_cpg_scores(phasebam, output_dir, name)


def DMR(outdir, fatherBam, motherBam, childBam, updRegions, threads, script_dir, REFERENCE):
    os.system(f'mkdir -p {outdir}')
    MethylationArgs = []
    MethylationArgs.append((f"{outdir}", fatherBam, updRegions, "father", threads, REFERENCE))
    MethylationArgs.append((f"{outdir}", motherBam, updRegions, "mother", threads, REFERENCE))
    MethylationArgs.append((f"{outdir}", childBam, updRegions, "child", threads, REFERENCE))

    processes = []
    for run in MethylationArgs:
        p = mp.Process(target=run_Methylation, args=run)
        p.start()
        processes.append(p)

    for p in processes:
        p.join()

    cmd = f'bash {DMR_SCRIPT} -u {updRegions} -f {outdir}/father/methylation_father.hap1.bed -g {outdir}/father/methylation_father.hap2.bed -m {outdir}/mother/methylation_mother.hap1.bed -n {outdir}/mother/methylation_mother.hap2.bed -c {outdir}/child/methylation_child.combined.bed -o {outdir}/dmr -t 0.3 -z {script_dir}'
    os.system(cmd)


# ============================================================
# 7) MAIN (UPDATED: joint VCF mode wiring)
# ============================================================
def main():
    global SAMTOOLS, REFERENCE, RSCRIPT, MLE_RSCRIPT, GZIP, PLOT_RSCRIPT, SEGMENTATION_RSCRIPT, MLE_PARENT_RSCRIPT, MLE_PARENT_CHILD_HOMOALT_RSCRIPT, DMR_SCRIPT, PLOT_RSCRIPT_PARENT, VBID

    father_in, mother_in, child_in, REFERENCE, snp_bed, thread, output_dir, prefix, runmode, upd, parent, UPDMode, DMRmode, downsample, joint_vcf, fatherSample, motherSample, childSample = argument_parser()

    # Determine input mode
    mode = "vcf_joint" if joint_vcf is not None else "bam"

    if mode == "bam":
        father_bam, mother_bam, child_bam = father_in, mother_in, child_in
    else:
        # vcf_joint uses joint_vcf + sample names
        pass

    # ---- mode guards ----
    if mode.startswith("vcf") and DMRmode:
        raise SystemExit("DMR mode requires BAM/CRAM inputs (methylation workflow).")

    if mode.startswith("vcf") and parent:
        raise SystemExit("--parent requires BAM/CRAM inputs (mpileup workflow).")

    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    if DMRmode and not UPDMode:
        print("can only run DMR mode if in UPD mode")
        return 19

    # configure paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    path_config = os.path.join(script_dir, 'path_config.json')

    SAMTOOLS, RSCRIPT, GZIP = get_paths(path_config)
    print(f'RSCRIPT: {RSCRIPT}')

    # R scripts
    MLE_RSCRIPT = os.path.join(script_dir, 'mle.R')
    MLE_PARENT_RSCRIPT = os.path.join(script_dir, 'mle_parent.R')
    MLE_PARENT_CHILD_HOMOALT_RSCRIPT = os.path.join(script_dir, 'mle_parent.R')

    PLOT_RSCRIPT = os.path.join(script_dir, 'plot_variant.R')
    PLOT_RSCRIPT_PARENT = os.path.join(script_dir, 'plot_variant_parent.R')

    SEGMENTATION_RSCRIPT1 = os.path.join(script_dir, 'upd_segmentation.R')
    SEGMENTATION_RSCRIPT2 = os.path.join(script_dir, 'upd_segmentation_2_R.r')

    DMR_SCRIPT = os.path.join(script_dir, 'DMRs.sh')

    # split regions
    segment_length = 50000000
    region_splits = split_regions(REFERENCE, segment_length)

    # optionally filter regions using SNP bed
    if snp_bed is not None:
        region_splits = filter_regions_with_snv(region_splits, snp_bed)

    combined_counts = os.path.join(output_dir, f'{prefix}.child.counts')

    if mode == "bam":
        # ----------------------------
        # BAM/CRAM mode: mpileup -> counts
        # ----------------------------
        tmp_region_dir = os.path.join(output_dir, prefix + '_tmp')
        os.system(f'mkdir -p {tmp_region_dir}')

        arg_list = []
        for region in region_splits:
            arg_list.append((father_bam, mother_bam, child_bam, region, tmp_region_dir, snp_bed, prefix))

        print('variant calling (mpileup)')
        with mp.Pool(thread) as pool:
            mpileup_files = pool.starmap(mpileup, arg_list)

        print('parsing mpileup')
        homoref_sampling_rate = 1 if snp_bed is not None else 0.01

        arg_list = [(mpf, homoref_sampling_rate) for mpf in mpileup_files]
        with mp.Pool(thread) as pool:
            counts_split_files = pool.starmap(get_child_count, arg_list)

        combine_count_files(counts_split_files, combined_counts)


    else:
        # ----------------------------
        # VCF mode (NEW): joint VCF -> counts
        # ----------------------------
        print('VCF mode (NEW): building counts table from joint multi-sample VCF')
        homoref_sampling_rate = 1.0 if snp_bed is not None else 0.01
        build_counts_from_joint_vcf(
            joint_vcf_path=joint_vcf,
            father_sample=fatherSample,
            mother_sample=motherSample,
            child_sample=childSample,
            out_counts_path=combined_counts,
            snp_bed=snp_bed,
            homoref_sampling_rate=homoref_sampling_rate
        )

    # plot variants
    print('plotting variants')
    run_plot_rscript(combined_counts, output_dir, REFERENCE, downsample=downsample)

    print(f'UPD mode {UPDMode}')
    if UPDMode:
        print('upd segmentation')
        upd_segment_length = 1000000  # 1mb
        run_segmentation(SEGMENTATION_RSCRIPT2, combined_counts, output_dir, segment_length=upd_segment_length)

        # Only run DMR if BAM mode + user requested dmrMode
        if mode == "bam" and DMRmode:
            UPDRegions = os.path.join(output_dir, os.path.basename(combined_counts) + ".upd.classification.bed")
            DMR(f"{output_dir}/DMR", father_bam, mother_bam, child_bam, UPDRegions, 32, script_dir, REFERENCE)
        else:
            print("Skipping DMR (requires BAM/CRAM + --dmrMode).")

    else:
        print('running MLE')
        run_mle_rscript(combined_counts, output_dir, runmode, upd)

        # upd segmentation
        print('upd segmentation')
        upd_segment_length = 1000000  # 1mb
        run_segmentation(SEGMENTATION_RSCRIPT1, combined_counts, output_dir, upd_segment_length)

        # get the chrX to autosome depth ratio
        child_x_ratio, father_x_ratio, mother_x_ratio = sexchrom_ratio(combined_counts)
        with open(os.path.join(output_dir, prefix + '.x2a.depth.tsv'), 'w') as f:
            f.write(f'child\t{child_x_ratio:.3f}\nfather\t{father_x_ratio:.3f}\nmother\t{mother_x_ratio:.3f}')

    # run parent DNA contamination if --parent is set (BAM mode only)
    print(f'Parent mode {parent}')
    if parent is True:
        if mode != "bam":
            raise SystemExit("--parent requires BAM/CRAM inputs (mpileup workflow).")
        with mp.Pool(thread) as pool:
            parent_counts_split_files = pool.map(get_counts_childhomoalt, mpileup_files)

        parent_counts = os.path.join(output_dir, f'{prefix}.parent.counts')
        combine_count_files(parent_counts_split_files, parent_counts)
        print('running MLE on parent contamination')
        run_mle_parent_child_homoalt_rscript(parent_counts, output_dir)
        print('plotting parent variants')
        run_plot_parent_rscript(parent_counts, output_dir, REFERENCE, downsample=downsample)

    print('done')


if __name__ == '__main__':
    python_version = sys.version_info[0] + 0.1 * sys.version_info[1]
    if python_version < 3.5:
        print('Triomix requires python version 3.5 or later.')
        print('You can specify the python interpreter such as python triomix -f father.bam -m mother.bam -c child.bam')
        print('Exiting...')
    main()