#!/usr/bin/env python
import sys
import logging
from collections import defaultdict
import re

"""
Keep only VCF records whose (CHR, POS) appear in the BED AND share ≥1 ID,
replace the VCF ID field with the BED IDs, then deduplicate records 
with the same replaced ID by keeping the best quality one.

BED format: CHR  POS  ID1;ID2;...
VCF IDs may be semicolon-delimited or '.'.

Notes:
- Chromosome names normalized (chr22 == 22).
- Position matching is exact (no window).
- Deduplication happens AFTER ID replacement.
"""

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s',
    stream=sys.stderr
)
logger = logging.getLogger(__name__)

def norm_chr(c):
    c = c.strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return c

def load_bed(bed_path):
    bed_map = defaultdict(set)  # (chr, pos) -> set(ids)
    with open(bed_path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            parts = ln.strip().split()
            if len(parts) < 2:
                continue
            chrom, pos = parts[0], parts[1]
            ids_str = parts[2] if len(parts) >= 3 else ""
            ids = {i for i in ids_str.split(";") if i and i != "."}
            if ids:
                bed_map[(norm_chr(chrom), int(pos))].update(ids)
    return bed_map

def vcf_id_set(id_field):
    if id_field == "." or id_field == "":
        return set()
    return {i for i in id_field.split(";") if i}

def extract_support_value(info_field):
    """Extract support value from INFO field with fallback options"""
    support = None
    
    # Try various support field formats
    support_patterns = [
        r'RE=(\d+)',           # cuteSV
        r'SUPPORT=(\d+)',      # Sniffles
        r'DR=(\d+)',           # Some callers use DR (supporting reads)
        r'DV=(\d+)',           # Some use DV (variant reads)
        r'AD=\d+,(\d+)',       # Allelic depth (take alt allele count)
    ]
    
    for pattern in support_patterns:
        match = re.search(pattern, info_field)
        if match:
            try:
                value = int(match.group(1))
                support = max(support or 0, value)
            except (ValueError, IndexError):
                continue
    
    return support 

def get_quality_score(cols):
    """Extract and return quality metrics for comparison"""
    qual_str, filter_field, info_field = cols[5], cols[6], cols[7]
    
    # Parse quality score
    try:
        qual = float(qual_str) if qual_str != '.' else 0.0
    except ValueError:
        qual = 0.0
    
    # Filter status (PASS is best)
    filter_priority = 1 if filter_field == "PASS" else 0
    
    # Support value
    support = extract_support_value(info_field)
    
    # Position (lower is better for tiebreaking)
    pos = int(cols[1])
    
    return (filter_priority, support, qual, -pos)  # negative pos for ascending sort

def main(vcf_in, bed_in, vcf_out):
    bed_map = load_bed(bed_in)

    # First pass: collect all valid records that pass the BED filter and replace IDs
    valid_records_with_replaced_ids = []
    total = 0
    
    with open(vcf_in) as fin:
        headers = []
        for ln in fin:
            if ln.startswith("#"):
                headers.append(ln)
                continue

            total += 1
            cols = ln.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue  # malformed, drop

            chrom, pos_str, vcf_id = cols[0], cols[1], cols[2]
            key = (norm_chr(chrom), int(pos_str))

            # require (chr,pos) in BED
            if key not in bed_map:
                continue

            bed_ids = bed_map[key]
            v_ids = vcf_id_set(vcf_id)

            # require at least one shared ID
            if not (v_ids & bed_ids):
                continue

            # Replace ID with the BED IDs
            replaced_id = ";".join(sorted(bed_ids))
            cols_copy = cols.copy()  # Make a copy to avoid modifying original
            cols_copy[2] = replaced_id
            
            quality_metrics = get_quality_score(cols_copy)
            valid_records_with_replaced_ids.append((cols_copy, quality_metrics, replaced_id))

    logger.info(f"Found {len(valid_records_with_replaced_ids)} valid records after BED filtering and ID replacement")

    # Log detailed record information at DEBUG level
    for i, (cols, quality_metrics, replaced_id) in enumerate(valid_records_with_replaced_ids):
        logger.debug(f"Record {i}: pos={cols[1]}, replaced_id='{replaced_id}', quality={quality_metrics}")

    # Second pass: deduplicate by replaced ID
    id_to_records = defaultdict(list)
    
    for i, (cols, quality_metrics, replaced_id) in enumerate(valid_records_with_replaced_ids):
        id_to_records[replaced_id].append((i, cols, quality_metrics))

    # Find and report duplicates
    duplicates_found = 0
    final_records = []
    
    for replaced_id, records in id_to_records.items():
        if len(records) > 1:
            duplicates_found += len(records) - 1
            logger.debug(f"Replaced ID '{replaced_id}' appears in {len(records)} records:")
            
            # Sort by quality (best first)
            records.sort(key=lambda x: x[2], reverse=True)
            
            for i, (record_idx, cols, quality_metrics) in enumerate(records):
                status = "KEEP" if i == 0 else "SKIP"
                logger.debug(f"  Record {record_idx}: pos={cols[1]}, quality={quality_metrics} -> {status}")
            
            # Keep only the best record
            final_records.append(records[0][1])  # records[0][1] is the cols of the best record
        else:
            # Only one record with this replaced ID
            final_records.append(records[0][1])

    # Write output
    kept = replaced = 0
    with open(vcf_out, "w") as fout:
        # Write headers
        for header in headers:
            fout.write(header)
        
        # Write deduplicated records (IDs already replaced)
        for cols in final_records:
            fout.write("\t".join(cols) + "\n")
            kept += 1
            replaced += 1

    logger.info(f"Examined {total} records; kept {kept} after BED filtering and deduplication; replaced IDs in {replaced}")
    logger.info(f"Deduplication: {len(valid_records_with_replaced_ids)} records after BED filtering -> {len(final_records)} unique records ({duplicates_found} duplicates removed)")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        logger.error(f"Usage: {sys.argv[0]} <input.vcf> <input.bed> <output.vcf>")
        sys.exit(1)
    
    # For debugging 
    # logging.getLogger().setLevel(logging.DEBUG)
    
    main(sys.argv[1], sys.argv[2], sys.argv[3])