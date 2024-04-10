import sys
import pysam
import numpy as np
import os
import re

def is_perfect_match(cigar_string):
    return re.fullmatch(r"\d+M", cigar_string) is not None

def is_autosome_or_chrX(chromosome_name):
    return re.match(r"^chr(\d+|X)$", chromosome_name) is not None

def gc_content(seq):
    gc_count = (seq.count('G') + seq.count('C')) / len(seq) if seq else 0
    return round(gc_count, 5)

def write_read_info(file, read, chromosome, frag_start, frag_end, read_name, gc_content_value):
    file.write(f"{chromosome}\t{frag_start}\t{frag_end}\t{read_name}\t{gc_content_value}\n")

def find_all_fragments(bam_file, out_dir):
    basename = os.path.splitext(os.path.basename(bam_file))[0]
    output_pos_pass = os.path.join(out_dir, f"{basename}_pos_pass.txt")
    output_neg_pass = os.path.join(out_dir, f"{basename}_neg_pass.txt")
    output_pos_fail = os.path.join(out_dir, f"{basename}_pos_fail.txt")
    output_neg_fail = os.path.join(out_dir, f"{basename}_neg_fail.txt")


    with pysam.AlignmentFile(bam_file, "rb") as bam, \
         open(output_pos_pass, "w") as pos_pass_file, \
         open(output_neg_pass, "w") as neg_pass_file, \
         open(output_pos_fail, "w") as pos_fail_file, \
         open(output_neg_fail, "w") as neg_fail_file:

        for read in bam:

            if read.is_unmapped:
                continue

            chromosome = read.reference_name
            tlen = read.template_length
            cigar_string = read.cigarstring
            base_qualities = read.query_qualities
            mean_quality = np.mean(base_qualities) if base_qualities else 0
            mapq = read.mapping_quality
            gc_content_value = gc_content(read.query_sequence)

            read_start = read.reference_start
            read_end = read.reference_end
            strand = '-' if read.is_reverse else '+'

            if strand == '+':
                frag_start = read_start # 0-based in pysam
                frag_end = frag_start + abs(tlen) # end position, BED-like (exclusive)
            elif strand == '-':
                frag_end = read_end # read_end is essentially 5' end of negative strand -> later used for motif calculation
                frag_start = frag_end - abs(tlen) # but because of bedtools convention, where start < end position, I am doing this for fragment

            flag = read.flag
            read_name = f"{read.query_name}_{flag}_{chromosome}_{frag_start+1}_{frag_end}"

            passes_filters = not read.is_unmapped and is_autosome_or_chrX(chromosome) and \
                             70 <= abs(tlen) <= 700 and is_perfect_match(cigar_string) and \
                             mean_quality >= 20 and mapq >= 40
            
            if passes_filters:
                output_file = pos_pass_file if strand == '+' else neg_pass_file
            else:
                output_file = pos_fail_file if strand == '+' else neg_fail_file

            write_read_info(output_file, read, chromosome, frag_start, frag_end, read_name, gc_content_value)

if __name__ == "__main__":
    bam_file = sys.argv[1]
    out_dir = sys.argv[2]  
    find_all_fragments(bam_file, out_dir)
