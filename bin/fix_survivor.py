#!/usr/bin/env python
"""
Script to process SURVIVOR VCF output:
- Use SUPP_VEC to identify contributing samples and extract their IDs
- Create combined ID field from contributing sample IDs
- Only consider samples marked in SUPP_VEC for selection
- Keep the best sample from the contributing ones (least NAs, most numeric data)
- Rename sample columns by removing suffixes
"""

import sys
import re
from collections import defaultdict

def parse_info_field(info_str):
    """Parse INFO field and return dictionary of key-value pairs"""
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict

def extract_sample_ids_from_supp_vec(sample_data, sample_names, supp_vec):
    """Extract sample IDs based on SUPP_VEC positions"""
    sample_ids = []
    tools = set()
    
    if not supp_vec or len(supp_vec) != len(sample_data):
        # Fallback: extract from all samples if SUPP_VEC is invalid
        supp_vec = '1' * len(sample_data)
    
    for i, (sample_name, sample_info, is_supported) in enumerate(zip(sample_names, sample_data, supp_vec)):
        # Only process samples that are marked as supported in SUPP_VEC
        if is_supported == '1':
            fields = sample_info.split(':')
            if len(fields) > 7:  # ID field is typically at index 7
                variant_id = fields[7]
                if variant_id not in ['NaN', 'NA', 'NAN', '.', '--', '']:
                    sample_ids.append(variant_id)
                    
                    # Extract tool name from ID
                    variant_id_lower = variant_id.lower()
                    if 'svim' in variant_id_lower:
                        tools.add('svim')
                    if 'sniffles' in variant_id_lower:
                        tools.add('sniffles')
                    if 'cutesv' in variant_id_lower:
                        tools.add('cutesv')
                    # if 'manta' in variant_id_lower:
                    #     tools.add('manta')
                    # if 'delly' in variant_id_lower:
                    #     tools.add('delly')
    
    return sample_ids, tools

def create_combined_id(sample_ids, tools, svtype):
    """Create combined ID from all contributing sample IDs"""
    if not sample_ids:
        # Fallback to tool-based naming if no sample IDs found
        if tools:
            sorted_tools = sorted(list(tools))
            tools_str = '_'.join(sorted_tools)
            return f"{tools_str}_{svtype}"
        else:
            return f"SURVIVOR_{svtype}"
    
    # Remove duplicates while preserving order
    unique_ids = []
    seen = set()
    for sample_id in sample_ids:
        if sample_id not in seen:
            unique_ids.append(sample_id)
            seen.add(sample_id)
    
    return ';'.join(unique_ids)

def count_missing_values(sample_info):
    """Count number of missing/NA values in sample info"""
    fields = sample_info.split(':')
    missing_count = 0
    
    for field in fields:
        if field in ['NaN', 'NA', 'NAN', '.', '--', '']:
            missing_count += 1
    
    return missing_count

def count_numeric_values(sample_info):
    """Count number of numeric values in sample info"""
    fields = sample_info.split(':')
    numeric_count = 0
    
    for field in fields:
        if field not in ['NaN', 'NA', 'NAN', '.', '--', '']:
            # Try to parse as number
            try:
                float(field)
                numeric_count += 1
            except ValueError:
                # Check if it contains numbers (like coordinates "chr1_95077-chr1_95077")
                if re.search(r'\d+', field):
                    numeric_count += 1
    
    return numeric_count

def count_informative_fields(sample_info):
    """Count number of fields with actual information"""
    fields = sample_info.split(':')
    info_count = 0
    
    for field in fields:
        if field not in ['NaN', 'NA', 'NAN', '.', '--', '']:
            info_count += 1
    
    return info_count

def get_quality_score(sample_info):
    """Extract quality score from sample info"""
    fields = sample_info.split(':')
    
    # Quality is typically at index 5 (QV field)
    if len(fields) > 5:
        try:
            qual = float(fields[5])
            return qual
        except (ValueError, TypeError):
            pass
    
    return 0

def score_sample_info(sample_info, debug=False):
    """Score sample based on least NAs and most numeric annotations"""
    # Count different types of information
    missing_count = count_missing_values(sample_info)
    numeric_count = count_numeric_values(sample_info)
    info_count = count_informative_fields(sample_info)
    quality = get_quality_score(sample_info)
    
    # Check if genotype is called
    has_genotype = not sample_info.startswith('./.')
    
    # Scoring priority:
    # 1. Has genotype call (major bonus)
    # 2. Fewer missing values (heavy penalty for NAs)
    # 3. More numeric annotations
    # 4. More informative fields overall
    # 5. Higher quality
    
    score = 0
    
    if has_genotype:
        score += 10000  # Major bonus for having a genotype call
    
    score -= missing_count * 1000  # Heavy penalty for missing data
    score += numeric_count * 100   # Bonus for numeric data
    score += info_count * 10       # Bonus for any informative field
    score += min(quality, 100)     # Quality bonus (capped)
    
    if debug:
        print(f"Sample: {sample_info[:50]}...")
        print(f"  Has genotype: {has_genotype}")
        print(f"  Missing count: {missing_count}")
        print(f"  Numeric count: {numeric_count}")
        print(f"  Info count: {info_count}")
        print(f"  Quality: {quality}")
        print(f"  Final score: {score}")
        print()
    
    return score, has_genotype, missing_count, numeric_count, info_count, quality

def find_best_sample_from_supp_vec(sample_data, sample_names, supp_vec, debug=False):
    """Find the best sample from only those marked in SUPP_VEC"""
    if not supp_vec or len(supp_vec) != len(sample_data):
        # Fallback: consider all samples if SUPP_VEC is invalid
        supp_vec = '1' * len(sample_data)
    
    # Filter to only consider samples marked in SUPP_VEC
    contributing_samples = []
    for i, (sample_name, sample_info, is_supported) in enumerate(zip(sample_names, sample_data, supp_vec)):
        if is_supported == '1':
            contributing_samples.append((i, sample_name, sample_info))
    
    if not contributing_samples:
        # Fallback: if no samples are marked, take the first one
        return 0, sample_names[0] if sample_names else "", sample_data[0] if sample_data else ""
    
    if debug:
        print(f"SUPP_VEC: {supp_vec}")
        print(f"Evaluating {len(contributing_samples)} contributing samples:")
    
    best_score = -float('inf')
    best_idx = contributing_samples[0][0]
    best_sample_name = contributing_samples[0][1]
    best_sample_data = contributing_samples[0][2]
    
    candidates = []
    
    # Score only the contributing samples
    for orig_idx, sample_name, sample_info in contributing_samples:
        score, has_genotype, missing_count, numeric_count, info_count, quality = score_sample_info(sample_info, debug)
        candidates.append((score, has_genotype, missing_count, numeric_count, info_count, quality, orig_idx, sample_name, sample_info))
        
        if debug:
            print(f"Contributing sample {orig_idx} ({sample_name}): score={score}")
    
    # Sort by: score (desc), has_genotype (desc), missing_count (asc), numeric_count (desc), info_count (desc), quality (desc)
    candidates.sort(key=lambda x: (-x[0], -x[1], x[2], -x[3], -x[4], -x[5]))
    
    if debug:
        print("Sorted contributing candidates:")
        for i, cand in enumerate(candidates):
            print(f"  {i+1}. Sample {cand[6]} ({cand[7]}): score={cand[0]}, genotype={cand[1]}, missing={cand[2]}, numeric={cand[3]}")
        print(f"Selected: Sample {candidates[0][6]} ({candidates[0][7]})")
        print()
    
    # Return the best candidate
    best = candidates[0]
    return best[6], best[7], best[8]  # index, name, data

def get_base_sample_name(sample_name):
    """Extract base sample name by keeping only the first two parts"""
    return '_'.join(sample_name.split('_')[:2])

def process_vcf(input_file, output_file, debug=False):
    """Process VCF file according to requirements"""
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sample_names = []
        
        for line_num, line in enumerate(infile, 1):
            line = line.strip()
            
            # Handle header lines
            if line.startswith('##'):
                outfile.write(line + '\n')
                continue
            
            # Handle column header line
            if line.startswith('#CHROM'):
                fields = line.split('\t')
                # Store original sample names (columns 9 onwards)
                sample_names = fields[9:]
                
                if debug:
                    print(f"Found sample names: {sample_names}")
                
                # Create new header with single sample column
                if sample_names:
                    base_name = get_base_sample_name(sample_names[0])
                    new_header = fields[:9] + [base_name]
                    outfile.write('\t'.join(new_header) + '\n')
                else:
                    outfile.write(line + '\n')
                continue
            
            # Process variant lines
            if not line.startswith('#'):
                fields = line.split('\t')
                
                if len(fields) < 9:
                    outfile.write(line + '\n')
                    continue
                
                chrom, pos, variant_id, ref, alt, qual, filter_field, info = fields[:8]
                format_field = fields[8] if len(fields) > 8 else ""
                sample_data = fields[9:] if len(fields) > 9 else []
                
                if debug and line_num <= 2:  # Debug first few variants
                    print(f"\nProcessing variant at line {line_num}: {chrom}:{pos}")
                    print(f"All sample data: {sample_data}")
                
                # Parse INFO field to get SVTYPE and SUPP_VEC
                info_dict = parse_info_field(info)
                svtype = info_dict.get('SVTYPE', 'UNK')
                supp_vec = info_dict.get('SUPP_VEC', '')
                
                # Extract sample IDs based on SUPP_VEC
                sample_ids, tools = extract_sample_ids_from_supp_vec(sample_data, sample_names, supp_vec)
                
                if debug and line_num <= 2:
                    print(f"Extracted sample IDs: {sample_ids}")
                
                # Create combined ID from contributing sample IDs
                new_id = create_combined_id(sample_ids, tools, svtype)
                
                if debug and line_num <= 2:
                    print(f"New ID: {new_id}")
                
                # Find best sample from only those marked in SUPP_VEC    
                if sample_data:
                    best_idx, best_sample_name, best_sample_data = find_best_sample_from_supp_vec(
                        sample_data, sample_names, supp_vec, debug and line_num <= 2
                    )
                    base_sample_name = get_base_sample_name(best_sample_name)

                    # Remove unwanted FORMAT/sample fields by name (e.g. ID, RAL, AAL, CO)
                    unwanted_fields = {'ID', 'RAL', 'AAL', 'CO'}

                    format_keys = format_field.split(':')
                    best_sample_fields = best_sample_data.split(':')

                    # Build new format and sample strings excluding unwanted fields
                    if len(format_keys) == len(best_sample_fields):
                        filtered = [
                            (f, v) for f, v in zip(format_keys, best_sample_fields)
                            if f not in unwanted_fields
                        ]
                        format_keys, best_sample_fields = zip(*filtered) if filtered else ([], [])

                    format_field_trimmed = ':'.join(format_keys)
                    best_sample_data_trimmed = ':'.join(best_sample_fields)

                    # Final output line
                    new_fields = [chrom, pos, new_id, ref, alt, qual, filter_field, info,
                                  format_field_trimmed, best_sample_data_trimmed]                    
                else:
                    new_fields = fields[:9]
                
                outfile.write('\t'.join(new_fields) + '\n')

def main():
    if len(sys.argv) < 3:
        print("Usage: python process_survivor_vcf.py <input.vcf> <output.vcf> [--debug]")
        print("\nThis script processes SURVIVOR VCF output to:")
        print("- Use SUPP_VEC to identify contributing samples")
        print("- Extract IDs from contributing samples and create combined ID")
        print("- Only consider samples marked in SUPP_VEC for selection")
        print("- Keep the best contributing sample (most info, least NAs)")
        print("- Rename sample columns by removing suffixes")
        print("\nSample selection criteria (in order of priority):")
        print("1. Has genotype call (not ./.)")
        print("2. Fewer missing/NA values")
        print("3. More numeric annotations")
        print("4. More informative fields overall")
        print("5. Higher quality scores")
        print("\nSupported tools: svim, sniffles, cutesv, manta, delly")
        print("\nUse --debug flag to see detailed scoring information")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    debug = len(sys.argv) > 3 and sys.argv[3] == '--debug'
    
    try:
        process_vcf(input_file, output_file, debug)
        print(f"Successfully processed {input_file} -> {output_file}")
    except Exception as e:
        print(f"Error processing VCF file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()