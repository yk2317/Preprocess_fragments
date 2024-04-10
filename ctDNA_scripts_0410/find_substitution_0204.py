import sys
import os
import pysam
import numpy as np
import re

def is_perfect_match(cigar_string):
    return re.fullmatch(r"\d+M", cigar_string) is not None

def gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq) if seq else 0

def find_substitution(bam_file, chromosome, position, ref_base, alt_base):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        pysam_position = position - 1  # Adjust for pysam's 0-based indexing
        for read in bam.fetch(chromosome, pysam_position, pysam_position + 1):
            if read.is_unmapped:
                continue
                
            read_sequence = read.query_sequence
            read_position = read.get_reference_positions(full_length=True)
            try:
                read_base_idx = read_position.index(pysam_position)
                read_base = read_sequence[read_base_idx]

                if read_base.upper() == alt_base.upper():
                    
                    read_start = read.reference_start
                    read_end = read.reference_end
                    tlen = read.template_length
                    strand = '-' if read.is_reverse else '+'

                    if strand == '+':
                        frag_start = read_start  # 0-based in pysam
                        frag_end = frag_start + abs(tlen)  # end position, BED-like (exclusive)
                    elif strand == '-':
                        frag_end = read_end  # read_end is 5' end of negative strand
                        frag_start = frag_end - abs(tlen)
                        
                    flag = read.flag
                    read_name = f"{read.query_name}_{flag}_{chromosome}_{frag_start + 1}"  # 1-based for downstream use

                    distance_to_start = read_base_idx
                    distance_to_end = len(read_sequence) - read_base_idx - 1
                    
                    cigar_string = read.cigarstring
                    clipped = 0 if is_perfect_match(cigar_string) else 1
                    gc_content_value = gc_content(read.query_sequence)
                    mapq = read.mapping_quality
                    base_qualities = read.query_qualities
                    mean_quality = np.mean(base_qualities) if base_qualities is not None else 0
                    consensus_min_depth = read.get_tag('cM') if read.has_tag('cM') else 0
                    consensus_error_rate = read.get_tag('cE') if read.has_tag('cE') else 0

                    print(f"{read.reference_name}\t{frag_start+1}\t{frag_end}\t{read_name}\t{strand}\t{tlen}\t{position}\t{ref_base}\t{alt_base}\t{distance_to_start}\t{distance_to_end}\t{clipped}\t{gc_content_value}\t{mapq}\t{mean_quality}\t{consensus_min_depth}\t{consensus_error_rate}")
            except ValueError:
                pass

def process_vcf(vcf_file, bam_file_path):
    sample_id = os.path.basename(bam_file_path).split("raw.bam")[0]

    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            parts = line.strip().split('\t')
            
            # Assuming the sample ID is now correctly matched to parts[7] or another column if needed
            if parts[7] == sample_id:  # Adjust this index based on where the sample ID is found
                chromosome = parts[0]
                position = int(parts[1])
                ref_base = parts[3]
                alt_base = parts[4]
                find_substitution(bam_file_path, chromosome, position, ref_base, alt_base)

if __name__ == "__main__":
    vcf_file = sys.argv[1]
    bam_file_path = sys.argv[2]

    process_vcf(vcf_file, bam_file_path)
