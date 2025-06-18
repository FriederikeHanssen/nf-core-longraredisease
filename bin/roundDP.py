import sys

def round_DP_to_integer(input_lines):
    output_lines = []
    for line in input_lines:
        if line.startswith('#'):
            # If it's a header line, remove leading and trailing spaces and append it to the output
            output_lines.append(line.strip() + '\n')
        else:
            # If it's a data line, split it by tabs
            parts = line.strip().split('\t')
            # Find the index of the DP field in the FORMAT column
            format_col = parts[8].split(':')
            dp_index = format_col.index('DP')
            # Split the DP values
            dp_values = parts[9].split(':')[dp_index].split(',')
            
            # Round DP values or preserve 'nan'
            rounded_dp_values = []
            for dp in dp_values:
                if dp.lower() == 'nan':
                    rounded_dp_values.append('nan')  # Preserve 'nan'
                else:
                    rounded_dp_values.append(str(round(float(dp))))  # Round the value
            
            # Join the rounded or preserved DP values back
            dp_values_str = ','.join(rounded_dp_values)
            
            # Replace the DP field with the updated values
            parts[9] = parts[9].replace(parts[9].split(':')[dp_index], dp_values_str)
            # Append the modified line to the output
            output_lines.append('\t'.join(parts) + '\n')
    return output_lines

def main():
    input_lines = sys.stdin.readlines()
    rounded_vcf_lines = round_DP_to_integer(input_lines)

    for line in rounded_vcf_lines:
        print(line, end='')

if __name__ == "__main__":
    main()

