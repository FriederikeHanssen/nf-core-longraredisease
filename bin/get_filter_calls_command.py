#!/usr/bin/env python
"""Construct SV call filtering command."""

import argparse
import sys
import subprocess
import logging

threshold_lookup = ['0'] + ['2'] * 10 + ['3'] * 9 + ['5'] * 20 + ['8'] * 100


def get_vcf_info_fields(vcf_path):
    """Get INFO field names from VCF header."""
    try:
        result = subprocess.run(['bcftools', 'view', '-h', vcf_path], 
                              capture_output=True, text=True, check=True)
        info_fields = []
        for line in result.stdout.split('\n'):
            if line.startswith('##INFO=<ID='):
                field_id = line.split('ID=')[1].split(',')[0]
                info_fields.append(field_id)
        return info_fields
    except subprocess.CalledProcessError:
        return []


def detect_support_field(vcf_path):
    """Detect the appropriate support field for different SV callers."""
    info_fields = get_vcf_info_fields(vcf_path)
    
    # Common support field names used by different SV callers
    support_field_candidates = [
        'SUPPORT',     # Sniffles
        'RE',          # CuteSV (supporting reads)
        'SUPPORT_READS', # Some versions
        'DR',          # SVIM (supporting reads)
        'DV',          # SVIM (variant reads)
        'READS',       # Generic
        'SR',          # Some callers use SR
    ]
    
    for candidate in support_field_candidates:
        if candidate in info_fields:
            return candidate
    
    # If no standard support field found, return None
    return None


def import_total_depth(path):
    """Get the average read depth."""
    with open(path, "r") as fh:
        for line in fh:
            if "total" in line:
                return float(line.strip().split("\t")[3])


def parse_arguments():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--vcf",
        required=True
    )

    parser.add_argument(
        "--target_bedfile",
        help=(
            "Provide the path to a bedfile containing regions"
            " in which to retain SV calls."
        ),
        required=False
    )

    parser.add_argument(
        "--depth_summary",
        help=(
            "Provide the path to a bedfile (e.g. from mosdepth)"
            " containing depth of coverage by region."
        ),
        required=True
    )

    parser.add_argument(
        "--bcftools_threads",
        help=(
            "Number of threads to use for bcftools view (default: 1)."
        ),
        default=1
    )

    parser.add_argument(
        "--min_read_support",
        help=(
            "Set the lower cutoff for read support, "
            "or 'auto' to enable autodetection."
        ),
        required=True,
    )

    parser.add_argument(
        "--min_read_support_limit",
        help=(
            "Set the absolute lower cutoff for read support to"
            " fall back to when using autodetection."
        ),
        required=True,
        type=int
    )

    parser.add_argument(
        "--contigs",
        help=(
            "Comma-separated list of contigs to keep."
        ),
        required=False,
        default=None
    )

    parser.add_argument(
        "--support_field",
        help=(
            "Specify the INFO field name for read support "
            "(auto-detected if not provided)."
        ),
        required=False,
        default=None
    )

    # NEW ARGUMENT: Make PASS filter optional
    parser.add_argument(
        "--filter_pass",
        help=(
            "Apply PASS filter to keep only variants that passed all filters. "
            "Use --filter_pass to enable, --no-filter_pass to disable."
        ),
        action='store_true',
        default=True  # Default to True (apply PASS filter)
    )

    parser.add_argument(
        "--no-filter_pass",
        help=(
            "Disable PASS filter (keep all variants regardless of FILTER field)."
        ),
        dest='filter_pass',
        action='store_false'
    )

    return parser.parse_args()


def main():
    """Run the entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format='[%(levelname)s] %(message)s'
    )

    args = parse_arguments()

    logging.info(f"Input VCF: {args.vcf}")
    logging.info(f"Read support setting: {args.min_read_support}")
    logging.info(f"Read support limit: {args.min_read_support_limit}")
    
    # Detect or use provided support field
    support_field = args.support_field
    if not support_field:
        support_field = detect_support_field(args.vcf)
        logging.info(f"Auto-detected support field: {support_field}")
    else:
        logging.info(f"Using user-specified support field: {support_field}")

    if not support_field:
        logging.error("Could not detect support field in VCF. Please specify --support_field manually.")
        sys.exit(1)

    # Determine read support threshold
    min_read_support = args.min_read_support_limit
    if args.min_read_support == 'auto':
        avg_depth = import_total_depth(args.depth_summary)
        if avg_depth is None:
            logging.warning("Could not extract average depth. Defaulting to upper threshold.")
            avg_depth = len(threshold_lookup) - 1
        else:
            logging.info(f"Extracted average depth: {avg_depth:.2f}")
        
        detected_read_support = int(threshold_lookup[round(avg_depth)])
        logging.info(f"Detected read support from lookup: {detected_read_support}")

        if detected_read_support > args.min_read_support_limit:
            min_read_support = detected_read_support
        else:
            logging.info(f"Detected threshold below limit. Using min_read_support_limit: {args.min_read_support_limit}")

    logging.info(f"Final minimum read support threshold: {min_read_support}")

    filter_min_read_support = f'INFO/{support_field} >= {min_read_support}'
    filter_string = f"-i '{filter_min_read_support}'"

    # Add optional filters
    if args.target_bedfile:
        filter_string = f"-T {args.target_bedfile} --targets-overlap 1 {filter_string}"
        logging.info(f"Filtering within target regions from BED file: {args.target_bedfile}")

    if args.contigs:
        filter_string = f"{filter_string} -r {args.contigs} "
        logging.info(f"Filtering by contigs: {args.contigs}")

    pass_filter = "-f PASS" if args.filter_pass else ""
    logging.info(f"Applying PASS filter: {args.filter_pass}")

    command = (
        f"bcftools view {pass_filter} --threads "
        f"{args.bcftools_threads} {filter_string} {args.vcf}"
    ).strip()

    logging.info(f"Final bcftools command: {command}")
    
    sys.stdout.write(command)


if __name__ == '__main__':
    main()