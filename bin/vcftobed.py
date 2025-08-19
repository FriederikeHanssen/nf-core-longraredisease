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

def merge_overlapping_intervals(intervals):
    """
    Merge overlapping and adjacent intervals.
    Input: list of (chrom, start, end) tuples
    Output: list of merged (chrom, start, end) tuples
    """
    if not intervals:
        return []
    
    # Sort by chromosome and start position
    sorted_intervals = sorted(intervals, key=lambda x: (x[0], x[1], x[2]))
    merged = [sorted_intervals[0]]
    
    for current in sorted_intervals[1:]:
        last = merged[-1]
        
        # If same chromosome and overlapping/adjacent (allowing 1bp gap)
        if current[0] == last[0] and current[1] <= last[2] + 1:
            # Merge by extending the end position
            merged[-1] = (last[0], last[1], max(last[2], current[2]))
        else:
            merged.append(current)
    
    return merged

def main(input_vcf, output_bed, verbose=False):
    """
    Extract all genomic intervals from SURVIVOR VCF file.
    Includes both main END intervals and all CO field intervals.
    Merges overlapping intervals automatically.
    """
    bed_set = set()
    
    with open(input_vcf, 'r') as vcf:
        for line_num, line in enumerate(vcf, 1):
            if line.startswith("#"):
                continue
                
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                if verbose:
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
                        if verbose:
                            print(f"Warning: Line {chrom}:{pos} Sample {sample_idx} has mismatched fields")
                        continue
                    
                    fmt_dict = dict(zip(format_keys, fmt_fields))
                    
                    # Extract CO intervals if present and valid
                    if "CO" in fmt_dict and fmt_dict["CO"] not in ("NAN", ".", "", "0"):
                        if verbose:
                            print(f"Extracting CO from {chrom}:{pos} Sample {sample_idx} -> {fmt_dict['CO']}")
                        
                        try:
                            intervals = extract_intervals_from_co(fmt_dict["CO"])
                            for interval in intervals:
                                bed_set.add(interval)
                        except Exception as e:
                            if verbose:
                                print(f"Error processing CO field '{fmt_dict['CO']}' at {chrom}:{pos}: {e}")
    
    # Convert set to list and merge overlapping intervals
    intervals_list = list(bed_set)
    merged_intervals = merge_overlapping_intervals(intervals_list)
    
    # Write merged BED output
    with open(output_bed, 'w') as bed:
        for chrom, start, end in merged_intervals:
            bed.write(f"{chrom}\t{start}\t{end}\n")
    
    print(f"Extracted {len(bed_set)} intervals, merged to {len(merged_intervals)} non-overlapping intervals")
    if verbose:
        print(f"Reduction: {len(bed_set) - len(merged_intervals)} overlapping intervals merged")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract and merge intervals from SURVIVOR VCF to BED format",
        epilog="This script processes SURVIVOR VCF files and extracts genomic intervals from both "
               "the main END field and sample-specific CO fields, automatically merging overlapping intervals."
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
    
    main(args.input, args.output, args.verbose)