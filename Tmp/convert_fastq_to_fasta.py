#!/usr/bin/env python3
import argparse
import os
import gzip
import re
import sys

def process_fastq(input_stream, output_stream, detect_direction, width):
    line_count = 0
    record_count = 0
    last_header = ""
    dot_count = 0  # Track dots printed in the current row

    for line in input_stream:
        line = line.strip()
        line_count += 1

        if line_count % 4 == 1:  # Header line
            record_count += 1
            
            if record_count % 1_000_000 == 0:
                sys.stderr.write(".")
                dot_count += 1
                if dot_count == 100:
                    sys.stderr.write("\n")
                    dot_count = 0
                sys.stderr.flush()
            
            if detect_direction:
                m = re.search(r'^@(\S+?)[ .](1|2)(\s+(.*))?$', line)
                if m:
                    read_id = m.group(1)  # Base read ID
                    pair_suffix = m.group(2) if m.group(2) else ""  # Paired-end indicator
                    metadata = m.group(4) if m.group(4) else ""  # Remaining metadata
                    last_header = f">{read_id}{'/' + pair_suffix if pair_suffix else ''} {metadata}".strip()
                else:
                    raise Exception(f"ERROR: FASTQ header line didn't match expected format: {line}")
            else:
                last_header = f">{line[1:]}"  # Remove leading '@' properly
        
        elif line_count % 4 == 2:  # Sequence line
            if width:
                output_stream.write(f"{last_header}\n{wrapped(line, width)}\n")
            else:
                output_stream.write(f"{last_header}\n{line}\n")
    
    if dot_count > 0:
        sys.stderr.write("\n")
    
    return record_count

def wrapped(string, every=60):
    string = string.rstrip()
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))

def main():
    parser = argparse.ArgumentParser(description='Convert FASTQ to FASTA format')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to an input file or directory')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to an output file or directory')
    parser.add_argument('-d', '--detect_direction', action='store_true', help='Auto-detect mate pair direction')
    parser.add_argument('-w', '--width', type=int, required=False, help='Defines the width of FASTA sequence lines')
    args = parser.parse_args()

    if os.path.isdir(args.input) and os.path.isdir(args.output):
        input_files = [f for f in os.listdir(args.input) if f.endswith(('.fq', '.fastq', '.fq.gz', '.fastq.gz'))]
        
        for filename in input_files:
            input_path = os.path.join(args.input, filename)
            basename = re.sub(r'(.fq|.fastq)(.gz)?$', '', filename)
            output_path = os.path.join(args.output, f"{basename}.fna")
            
            sys.stderr.write(f"Processing file '{filename}'\n")
            sys.stderr.flush()
            
            with (gzip.open(input_path, 'rt') if filename.endswith('.gz') else open(input_path, 'r')) as input_stream, open(output_path, 'w') as output_stream:
                record_count = process_fastq(input_stream, output_stream, args.detect_direction, args.width)
                print(f"{record_count} records written to {output_path}", file=sys.stderr)
    else:
        with (gzip.open(args.input, 'rt') if args.input.endswith('.gz') else open(args.input, 'r')) as input_stream, open(args.output, 'w') as output_stream:
            record_count = process_fastq(input_stream, output_stream, args.detect_direction, args.width)
            print(f"{record_count} records written to {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()
