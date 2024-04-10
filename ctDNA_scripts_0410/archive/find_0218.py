import os
import sys
import numpy as np
import re
import pysam

def is_perfect_match(cigar_string):
    return re.fullmatch(r"\d+M", cigar_string) is not None

def is_autosome_or_chrX(chromosome_name):
    return re.match(r"^chr(\d+|X)$", chromosome_name) is not None

def gc_content(seq):
    gc_count = (seq.count('G') + seq.count('C')) / len(seq) if seq else 0
    return round(gc_count, 5)

def write_read_info(file, read, chromosome, frag_start, frag_end, read_name, gc_content_value, distance_to_start, distance_to_end):
    file.write(f"{chromosome}\t{frag_start}\t{frag_end}\t{read_name}\t{gc_content_value}\t{distance_to_start}\t{distance_to_end}\n")
    file.flush()  # Explicitly flush the file buffer
    print("Data written to file.")

def find_substitution(bam_file, out_dir, chromosome, position, ref_base, alt_base, output_files):
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
                    distance_to_start = read_base_idx
                    distance_to_end = len(read_sequence) - read_base_idx - 1

                    chromosome = read.reference_name
                    tlen = read.template_length
                    cigar_string = read.cigarstring
                    base_qualities = read.query_qualities
                    mean_quality = np.mean(base_qualities) if base_qualities else 0
                    mapq = read.mapping_quality
                    gc_content_value = gc_content(read.query_sequence)
                    consensus_min_depth = read.get_tag('cM') if read.has_tag('cM') else 0
                    consensus_error_rate = read.get_tag('cE') if read.has_tag('cE') else 0

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
                        mean_quality >= 20 and mapq >= 40 and \
                        consensus_min_depth >= 5 and consensus_error_rate < 0.03

                    if passes_filters:
                        output_file = output_files['pos_pass_file'] if strand == '+' else output_files['neg_pass_file']
                    else:
                        output_file = output_files['pos_fail_file'] if strand == '+' else output_files['neg_fail_file']

                    write_read_info(output_file, read, chromosome, frag_start, frag_end, read_name, gc_content_value, distance_to_start, distance_to_end)
            except ValueError:
                pass

def process_vcf(vcf_file, bam_file_path, out_dir):
    sample_id = os.path.basename(bam_file_path).split("raw.bam")[0]
    output_files = {
        'pos_pass_file': None,
        'neg_pass_file': None,
        'pos_fail_file': None,
        'neg_fail_file': None
    }

    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            parts = line.strip().split('\t')

            if parts[7] == sample_id:  # Adjust this index based on where the sample ID is found
                chromosome = parts[0]
                position = int(parts[1])
                ref_base = parts[3]
                alt_base = parts[4]

                basename = os.path.splitext(os.path.basename(bam_file_path))[0]
                output_files['pos_pass_file'] = open(os.path.join(out_dir, f"{basename}_pos_pass.txt"), "w")
                output_files['neg_pass_file'] = open(os.path.join(out_dir, f"{basename}_neg_pass.txt"), "w")
                output_files['pos_fail_file'] = open(os.path.join(out_dir, f"{basename}_pos_fail.txt"), "w")
                output_files['neg_fail_file'] = open(os.path.join(out_dir, f"{basename}_neg_fail.txt"), "w")

                find_substitution(bam_file_path, out_dir, chromosome, position, ref_base, alt_base, output_files)

                for file in output_files.values():
                    file.close()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py vcf_file bam_file_path out_dir")
        sys.exit(1)

    vcf_file = sys.argv[1]
    bam_file_path = sys.argv[2]
    out_dir = sys.argv[3]

    process_vcf(vcf_file, bam_file_path, out_dir)
