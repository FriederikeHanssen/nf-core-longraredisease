#!/usr/bin/env python3
import argparse
import re

def extract_intervals_from_co(co_string):
    """
    Extract genomic intervals from CO field string.
    Handles comma-separated regions like: chr22_13704-chr22_13741,chr22_13730-chr22_13851
    """
    intervals = []
    for region in co_string.split(","):
        region = region.strip()  # Remove any whitespace
        # Match pattern: chr[name]_start-chr[name]_end or chr[name]_start-end
        m = re.match(r"^(chr[^\s_]+)_(\d+)-(chr[^\s_]+_)?(\d+)$", region)
        if m:
            chrom = m.group(1)
            start = int(m.group(2)) - 1  # Convert to 0-based BED format
            end = int(m.group(4))
            intervals.append((chrom, start, end))
        else:
            print(f"Warning: Could not parse CO region: {region}")
    return intervals

def main(input_vcf, output_bed):
    """
    Extract all genomic intervals from SURVIVOR VCF file.
    Includes both main END intervals and all CO field intervals.
    """
    bed_set = set()
    
    with open(input_vcf, 'r') as vcf:
        for line_num, line in enumerate(vcf, 1):
            if line.startswith("#"):
                continue
                
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                print(f"Warning: Line {line_num} has insufficient fields, skipping")
                continue
                
            chrom = fields[0]
            pos = int(fields[1])
            info = fields[7]
            
            # Extract main END from INFO field
            end_match = re.search(r"END=(\d+)", info)
            if end_match:
                end = int(end_match.group(1))
                bed_set.add((chrom, pos-1, end))  # Convert to 0-based
            
            # Process sample-specific CO fields if present
            if len(fields) > 9:
                format_keys = fields[8].split(":")
                
                for sample_idx, sample_field in enumerate(fields[9:]):
                    fmt_fields = sample_field.split(":")
                    
                    if len(fmt_fields) != len(format_keys):
                        print(f"Warning: Line {chrom}:{pos} Sample {sample_idx} has mismatched fields: {sample_field}")
                        continue
                    
                    fmt_dict = dict(zip(format_keys, fmt_fields))
                    
                    # Extract CO intervals if present and valid
                    if "CO" in fmt_dict and fmt_dict["CO"] not in ("NAN", ".", "", "0"):
                        print(f"Extracting CO from {chrom}:{pos} Sample {sample_idx} -> {fmt_dict['CO']}")
                        
                        try:
                            intervals = extract_intervals_from_co(fmt_dict["CO"])
                            for interval in intervals:
                                bed_set.add(interval)
                        except Exception as e:
                            print(f"Error processing CO field '{fmt_dict['CO']}' at {chrom}:{pos}: {e}")
    
    # Write sorted BED output
    with open(output_bed, 'w') as bed:
        for chrom, start, end in sorted(bed_set, key=lambda x: (x[0], x[1], x[2])):
            bed.write(f"{chrom}\t{start}\t{end}\n")
    
    print(f"Extracted {len(bed_set)} unique intervals to {output_bed}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract all main and CO intervals from SURVIVOR VCF to BED format",
        epilog="This script processes SURVIVOR VCF files and extracts genomic intervals from both "
               "the main END field and sample-specific CO fields, including comma-separated regions."
    )
    parser.add_argument("-i", "--input", required=True, 
                       help="Input SURVIVOR VCF file")
    parser.add_argument("-o", "--output", required=True, 
                       help="Output BED file")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Enable verbose output")
    
    args = parser.parse_args()
    
    if args.verbose:
        print(f"Processing VCF file: {args.input}")
        print(f"Output BED file: {args.output}")
    
    main(args.input, args.output)