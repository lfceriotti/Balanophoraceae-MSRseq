#!/usr/bin/env python3

import sys
import csv
import re
from collections import defaultdict

def load_anticodon_positions(annotation_file):
    anticodon_positions = defaultdict(set)
    with open(annotation_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["Note"].strip().lower() == "anticodon":
                gene = row["GeneName"].replace(".trna1", "")
                pos = int(row["Position"])
                anticodon_positions[gene].add(pos)
    return anticodon_positions

def load_expected_mismatches(mismatch_file):
    expected = {}
    with open(mismatch_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row["GeneName"].replace(".trna1", "")
            ref_seq = row["Anticodon"].strip().upper()
            alt_seq = row["Anticodon_modified"].strip().upper()
            expected[gene] = (ref_seq, alt_seq)
    return expected

def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

def parse_md_and_cigar(cigar, md_tag, pos):
    ref_to_read = {}
    mismatches = []
    deletions = []

    read_index = 0
    ref_pos = pos

    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    md_parts = re.findall(r'(\d+)|(\^[A-Z]+)|([A-Z])', md_tag)

    md_stream = []
    for m in md_parts:
        if m[0]:
            md_stream.extend(['='] * int(m[0]))
        elif m[1]:
            md_stream.extend(['^' + b for b in m[1][1:]])
        elif m[2]:
            md_stream.append(m[2])

    md_ptr = 0

    for length_str, op in cigar_ops:
        length = int(length_str)
        if op == 'M':
            for _ in range(length):
                if md_ptr >= len(md_stream):
                    break
                state = md_stream[md_ptr]
                if state.startswith('^'):
                    deletions.append((ref_pos, state[1]))
                    ref_pos += 1
                elif state == '=':
                    ref_to_read[ref_pos] = read_index
                    ref_pos += 1
                    read_index += 1
                else:
                    ref_to_read[ref_pos] = read_index
                    mismatches.append((ref_pos, state))
                    ref_pos += 1
                    read_index += 1
                md_ptr += 1
        elif op == 'I':
            for _ in range(length):
                read_index += 1
        elif op == 'D':
            pass
        elif op == 'S':
            read_index += length

    return ref_to_read, mismatches, deletions

def filter_sam_by_anticodon_mismatch(sam_file, output_file, anticodon_positions, expected_mismatches):
    with open(sam_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith("@"):
                outfile.write(line)
                continue

            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue

            qname, flag, rname, pos, _, cigar, _, _, _, seq, qual, *opt = fields
            pos = int(pos)
            flag = int(flag)
            gene = rname.replace(".trna1", "")
            is_reverse = (flag & 0x10) != 0
            seq = seq.upper()
            if is_reverse:
                seq = reverse_complement(seq)

            if gene not in anticodon_positions or gene not in expected_mismatches:
                continue

            anti_pos = sorted(anticodon_positions[gene])
            ref_acodon, alt_acodon = expected_mismatches[gene]
            if len(anti_pos) != 3 or len(ref_acodon) != 3 or len(alt_acodon) != 3:
                continue
            ref2anticodon = dict(zip(anti_pos, zip(ref_acodon, alt_acodon)))

            md_tag = next((x[5:] for x in opt if x.startswith("MD:Z:")), None)
            if not md_tag:
                continue

            ref_to_read, mismatches, deletions = parse_md_and_cigar(cigar, md_tag, pos)

            # Skip read if any deletion in anticodon
            if any(ref in ref2anticodon for ref, _ in deletions):
                continue

            total_anticodon_mismatches = 0
            matching_expected = 0

            for ref_pos, ref_base in mismatches:
                if ref_pos not in ref2anticodon:
                    continue

                total_anticodon_mismatches += 1
                read_index = ref_to_read.get(ref_pos, -1)
                if read_index == -1 or read_index >= len(seq):
                    continue

                read_base = seq[read_index]
                _, expected_base = ref2anticodon[ref_pos]
                if read_base == expected_base:
                    matching_expected += 1

            if total_anticodon_mismatches == 1 and matching_expected == 1:
                outfile.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: script.py annotation.tsv mismatches.tsv input.sam output.sam")
        sys.exit(1)

    annotation_file = sys.argv[1]
    mismatch_file = sys.argv[2]
    input_sam = sys.argv[3]
    output_sam = sys.argv[4]

    anticodon_positions = load_anticodon_positions(annotation_file)
    expected_mismatches = load_expected_mismatches(mismatch_file)

    filter_sam_by_anticodon_mismatch(input_sam, output_sam, anticodon_positions, expected_mismatches)
