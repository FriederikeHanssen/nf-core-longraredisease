#!/usr/bin/env python
"""
Enhanced SURVIVOR VCF processor:
- Uses SUPP_VEC to identify contributing samples
- Extracts variant IDs from contributing samples to make a combined ID
- Chooses the best contributing sample (most info, least NAs)
- Expands output to one line per 'CO' interval from contributing samples
- Removes unwanted FORMAT fields (ID, RAL, AAL, CO)
- Renames sample columns by removing suffixes
"""

import sys
import re

def parse_info_field(info_str):
    """Parse INFO field into dict."""
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict

def extract_sample_ids_from_supp_vec(sample_data, sample_names, supp_vec):
    """Extract variant IDs and tool names from SUPP_VEC contributing samples."""
    sample_ids = []
    tools = set()
    if not supp_vec or len(supp_vec) != len(sample_data):
        supp_vec = '1' * len(sample_data)  # Fallback: all samples

    for sample_name, sample_info, is_supported in zip(sample_names, sample_data, supp_vec):
        if is_supported != '1':
            continue
        fields = sample_info.split(':')
        if len(fields) > 7:
            variant_id = fields[7]
            if variant_id not in ['NaN', 'NA', 'NAN', '.', '--', '']:
                sample_ids.append(variant_id)
                lowid = variant_id.lower()
                if 'svim' in lowid:
                    tools.add('svim')
                if 'sniffles' in lowid:
                    tools.add('sniffles')
                if 'cutesv' in lowid:
                    tools.add('cutesv')
    return sample_ids, tools

def create_combined_id(sample_ids, tools, svtype):
    """Make combined ID from sample IDs or tool names."""
    if not sample_ids:
        return f"{'_'.join(sorted(tools))}_{svtype}" if tools else f"SURVIVOR_{svtype}"
    unique_ids = []
    seen = set()
    for sid in sample_ids:
        if sid not in seen:
            unique_ids.append(sid)
            seen.add(sid)
    return ';'.join(unique_ids)

def score_sample_info(sample_info):
    """Score a sample for best selection."""
    fields = sample_info.split(':')
    missing = sum(f in ['NaN', 'NA', 'NAN', '.', '--', ''] for f in fields)
    numeric = sum(_is_numeric(f) for f in fields)
    info_count = sum(f not in ['NaN', 'NA', 'NAN', '.', '--', ''] for f in fields)
    qual = _safe_float(fields[5]) if len(fields) > 5 else 0
    has_gt = not sample_info.startswith('./.')
    score = 0
    if has_gt:
        score += 10000
    score -= missing * 1000
    score += numeric * 100
    score += info_count * 10
    score += min(qual, 100)
    return score

def _is_numeric(val):
    try:
        float(val)
        return True
    except ValueError:
        return bool(re.search(r'\d', val))

def _safe_float(val):
    try:
        return float(val)
    except ValueError:
        return 0

def find_best_sample(sample_data, sample_names, supp_vec):
    """Find best contributing sample from SUPP_VEC."""
    if not supp_vec or len(supp_vec) != len(sample_data):
        supp_vec = '1' * len(sample_data)
    contrib = [(i, n, d) for i, (n, d, s) in enumerate(zip(sample_names, sample_data, supp_vec)) if s == '1']
    if not contrib:
        return 0, sample_names[0], sample_data[0]
    scored = [(score_sample_info(d), i, n, d) for i, n, d in contrib]
    scored.sort(key=lambda x: (-x[0], x[1]))
    _, idx, name, data = scored[0]
    return idx, name, data

def get_base_sample_name(sample_name):
    """Shorten sample name to first two underscore parts."""
    return '_'.join(sample_name.split('_')[:2])

def extract_co_intervals(sample_data, sample_names, supp_vec):
    """Extract CO intervals from contributing SUPP_VEC samples."""
    intervals = []
    if not supp_vec or len(supp_vec) != len(sample_data):
        supp_vec = '1' * len(sample_data)
    for sample_name, sample_info, is_supported in zip(sample_names, sample_data, supp_vec):
        if is_supported != '1':
            continue
        fields = sample_info.split(':')
        if len(fields) > 8:  # CO usually last
            co_val = fields[-1]
            if co_val not in ['NaN', 'NA', 'NAN', '.', '--', '']:
                for part in co_val.split(','):
                    m = re.match(r"([^_]+)_(\d+)-([^_]+)_(\d+)", part)
                    if m:
                        intervals.append((m.group(1), m.group(2)))
    return intervals

def process_vcf(input_file, output_file, debug=False):
    """Main processing loop."""
    with open(input_file) as infile, open(output_file, 'w') as outfile:
        sample_names = []
        for line in infile:
            line = line.strip()
            if line.startswith('##'):
                outfile.write(line + '\n')
                continue
            if line.startswith('#CHROM'):
                fields = line.split('\t')
                sample_names = fields[9:]
                base_name = get_base_sample_name(sample_names[0]) if sample_names else ""
                outfile.write('\t'.join(fields[:9] + ([base_name] if base_name else [])) + '\n')
                continue

            fields = line.split('\t')
            if len(fields) < 9:
                outfile.write(line + '\n')
                continue

            chrom, pos, _, ref, alt, qual, filt, info, fmt = fields[:9]
            sample_data = fields[9:]

            info_dict = parse_info_field(info)
            svtype = info_dict.get('SVTYPE', 'UNK')
            supp_vec = info_dict.get('SUPP_VEC', '')

            sample_ids, tools = extract_sample_ids_from_supp_vec(sample_data, sample_names, supp_vec)
            new_id = create_combined_id(sample_ids, tools, svtype)

            best_idx, best_name, best_data = find_best_sample(sample_data, sample_names, supp_vec)
            base_name = get_base_sample_name(best_name)

            unwanted = {'ID', 'RAL', 'AAL', 'CO'}
            fmt_keys = fmt.split(':')
            best_vals = best_data.split(':')
            if len(fmt_keys) == len(best_vals):
                fmt_keys, best_vals = zip(*[(f, v) for f, v in zip(fmt_keys, best_vals) if f not in unwanted])
            fmt_trim = ':'.join(fmt_keys)
            vals_trim = ':'.join(best_vals)

            intervals = extract_co_intervals(sample_data, sample_names, supp_vec)
            if intervals:
                for ch, ps in intervals:
                    outfile.write('\t'.join([ch, ps, new_id, ref, alt, qual, filt, info, fmt_trim, vals_trim]) + '\n')
            else:
                outfile.write('\t'.join([chrom, pos, new_id, ref, alt, qual, filt, info, fmt_trim, vals_trim]) + '\n')

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <input.vcf> <output.vcf> [--debug]")
        sys.exit(1)
    process_vcf(sys.argv[1], sys.argv[2], '--debug' in sys.argv)

if __name__ == "__main__":
    main()
