#!/usr/bin/env python3

import sys

def round_depth_to_integer(input_lines):
    output_lines = []
    for line in input_lines:
        if line.startswith('#'):
            output_lines.append(line.rstrip('\n') + '\n')
        else:
            parts = line.strip().split('\t')
            format_col = parts[8].split(':')

            # Look for DP or CD in the FORMAT column. Otherwise gave me an error
            if 'DP' in format_col:
                depth_field = 'DP'
            elif 'CD' in format_col:
                depth_field = 'CD'
            else:
                output_lines.append(line + '\n')
                continue  # No depth field to process

            depth_index = format_col.index(depth_field)

            # Get the values, split if comma-separated (should handle most VCFs)
            sample_values = parts[9].split(':')
            depth_values = sample_values[depth_index].split(',')

            # Round, preserving 'nan' - was there in the original code
            rounded_values = []
            for val in depth_values:
                if val.lower() == 'nan':
                    rounded_values.append('nan')
                else:
                    try:
                        rounded_values.append(str(round(float(val))))
                    except ValueError:
                        rounded_values.append(val)  # If value can't be converted, leave as is

            sample_values[depth_index] = ','.join(rounded_values)
            parts[9] = ':'.join(sample_values)
            output_lines.append('\t'.join(parts) + '\n')
    return output_lines

def main():
    input_lines = sys.stdin.readlines()
    rounded_vcf_lines = round_depth_to_integer(input_lines)

    for line in rounded_vcf_lines:
        print(line, end='')

if __name__ == "__main__":
    main()
