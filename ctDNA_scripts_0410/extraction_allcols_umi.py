import sys
import pysam
import numpy as np
import os
import re

def gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq) if seq else 0

def is_perfect_match(cigar_string):
    return re.fullmatch(r"\d+M", cigar_string) is not None

def is_autosome_or_chrX(chromosome_name):
    return re.match(r"^chr(\d+|X)$", chromosome_name) is not None

def find_all_fragments(bam_file):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            # Skip unmapped reads
            if read.is_unmapped:
                continue

            chromosome = read.reference_name
            read_start = read.reference_start
            read_end = read.reference_end
            tlen = read.template_length
            strand = '-' if read.is_reverse else '+'

            if strand == '+':
                frag_start = read_start  # 0-based in pysam
                frag_end = frag_start + abs(tlen)  # end position, BED-like (exclusive)
            elif strand == '-':
                frag_end = read_end  # read_end is essentially 5' end of negative strand -> later used for motif calculation
                frag_start = frag_end - abs(tlen)  # but because of bedtools convention, where start < end position, I am doing this for fragment

            flag = read.flag
            read_name = f"{read.query_name}_{flag}_{chromosome}_{frag_start+1}"  # 1-based for downstream use
            umi = read.get_tag('RX') if read.has_tag('RX') else 'NA'  # Extract UMI

            cigar_string = read.cigarstring
            clipped = 0 if is_perfect_match(cigar_string) else 1
            gc_content_value = gc_content(read.query_sequence)
            mapq = read.mapping_quality
            base_qualities = read.query_qualities
            mean_quality = np.mean(base_qualities) if base_qualities else 0
            consensus_min_depth = read.get_tag('cM') if read.has_tag('cM') else 0
            consensus_error_rate = read.get_tag('cE') if read.has_tag('cE') else 0

            print(f"{chromosome}\t{frag_start}\t{frag_end}\t{read_name}\t{strand}\t{read_start}\t{read_end}\t{tlen}\t{clipped}\t{cigar_string}\t{gc_content_value}\t{mapq}\t{mean_quality}\t{consensus_min_depth}\t{consensus_error_rate}\t{umi}")

if __name__ == "__main__":
    bam_file = sys.argv[1]
    find_all_fragments(bam_file)
