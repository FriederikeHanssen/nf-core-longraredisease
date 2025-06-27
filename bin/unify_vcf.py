#!/usr/bin/env python
"""Unify VCFs for tertiary analyses."""

import gzip
import os
import argparse

# Assuming these are available in your environment
# from .util import get_named_logger, wf_parser  # noqa: ABS101

# Simple logger for standalone use
import logging
logger = logging.getLogger("unify_vcf")
logging.basicConfig(level=logging.INFO)


def get_lines(file_path):
    """
    Return lines in a file.

    This is done in order to handle both zipped and unzipped files.
    """
    try:
        with gzip.open(file_path, "rt") as file_h:
            lines = file_h.readlines()
    except gzip.BadGzipFile:
        with open(file_path) as file_h:
            lines = file_h.readlines()
    return lines


def modify_repeat_line(repeat_line):
    """
    Add SVTYPE=REP to INFO in the vcf.

    This is necessary for repeat vcf for providers that are not ONT, to enable
    Geneyx application reading these lines as REP effect.
    ONT repeat files already contain SVTYPE=STR in their INFO section, which is
    also recognized by Geneyx as REP effect
    """
    line_parts = repeat_line.strip().split("\t")
    if len(line_parts) >= 8:
        chrom, pos, vid, ref, alt, qual, vfilter, info = line_parts[:8]
        remaining = line_parts[8:] if len(line_parts) > 8 else []
        
        info = info + ";SVTYPE=REP"
        repeat_line = "\t".join([chrom, pos, vid, ref, alt, qual, vfilter, info] + remaining) + "\n"
    
    return repeat_line


def is_valid_line(line):
    """
    Check that the line contains SVTYPE=.

    this is Geneyx way to recognizing structural variants
    also checks that the gt (genotype) is not "./." (which means a no-call)
    a line is valid if it contains SVTYPE= and it had a valid genotype
    """
    # check if line is header or null / empty
    if line.startswith("#") or not line.strip():
        return True

    cells = line.split("\t")
    if len(cells) < 10:
        return True  # Skip malformed lines
        
    info_cell = cells[7]
    format_key = cells[8]
    format_value = cells[9]

    # if the "svtype" line exists the line is valid
    if "SVTYPE=" in info_cell:
        return True

    try:
        gt_index = format_key.split(":").index("GT")
        gt_value = format_value.split(":")[gt_index]

        # if "svtype" not present and GT is 0 or "./." the line is invalid
        if gt_value == "0/0" or gt_value == "./.":
            return False
    except (ValueError, IndexError):
        # If GT field is not found or malformed, consider line valid
        pass

    return True


def write_vcf_content(ftype, lines_dict, output_h, skip_svtype):
    """
    Write VCF.

    Writes the content of a file into the output file skips the header
    of the file modifies lines when necessary
    """
    # Handle multiple files for the same type (e.g., multiple SV files)
    if isinstance(lines_dict[ftype], list):
        files_to_process = lines_dict[ftype]
    else:
        files_to_process = [lines_dict[ftype]] if lines_dict[ftype] is not None else []
    
    for file_lines in files_to_process:
        if file_lines is None:
            continue
            
        start_read_content_flag = False
        for vcf_line in file_lines:
            if start_read_content_flag:
                if is_valid_line(vcf_line):
                    if ftype == "repeat" and not skip_svtype:
                        vcf_line = modify_repeat_line(vcf_line)
                    output_h.write(vcf_line)
            if vcf_line.startswith("#CHROM"):
                start_read_content_flag = True

        if not start_read_content_flag:
            logger.warning(f"No vcf header for given {ftype} file")


def get_vcf_header_from_lines(vcf_lines):
    """Get VCF lines and returns the header."""
    vcf_header = ""
    for line in vcf_lines:
        vcf_header = vcf_header + line
        if line.find("#CHROM") != -1:
            break

    return vcf_header


def is_empty(lines):
    """Check if the file is empty."""
    for line in lines:
        striped_line = line.strip()
        if striped_line != "" and not striped_line.startswith("#"):
            return False

    return True


def combine_headers(files_lines):
    """Combine the headers of multiple VCF files."""
    c_line = None
    h_lines = []
    printed_format = False

    # Process each file type
    for key in ['sv', 'cnv', 'repeat']:
        files_for_type = files_lines[key]
        if files_for_type is None:
            continue
            
        # Handle multiple files for the same type
        if isinstance(files_for_type, list):
            files_to_process = files_for_type
        else:
            files_to_process = [files_for_type]
            
        for lines in files_to_process:
            if lines:
                for line in lines:
                    if line.startswith('##') and line not in h_lines:
                        if line.startswith('##fileformat') and not printed_format:
                            h_lines.append(line)
                            printed_format = True
                            continue
                        elif not line.startswith('##fileformat'):
                            h_lines.append(line)
                    elif line.startswith('#CHROM') and not c_line:
                        c_line = line
    
    if c_line:
        h_lines.append(c_line)
    return h_lines


def create_unified_file(files_lines, output_path, skip_svtype):
    """
    Create the unified file.

    Prints the sv file complete (with header), then the cnv without a header
    and repeats without a header and with modified info section
    if sv file is empty prints all cnv lines (with header) and repeat
    lines without a header if cnv is also empty prints repeat lines with the
    header
    """
    # We need to combine the headers first, to avoid issues downstream when
    # sorting/indexing with bcftools and keep the output compliant.
    header = combine_headers(files_lines)

    # Then, combine the body of the files
    with open(output_path, "w") as output_h:
        # Save the combined header first
        if header:
            output_h.write("".join(header))

        # Process files in order: sv, cnv, repeat
        for ftype in ['sv', 'cnv', 'repeat']:
            if files_lines[ftype] is not None:
                write_vcf_content(ftype, files_lines, output_h, skip_svtype)


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser(description="Unify VCF files")
    parser.add_argument(
        '-o', '--outputPath',
        help='the unified output VCF path (required)',
        required=True
    )
    parser.add_argument(
        '-s', '--svPaths',
        help='SV input file paths (optional, can specify multiple)',
        required=False, 
        nargs='*',
        default=[]
    )
    parser.add_argument(
        '--svPath',
        help='Single SV input file path (optional, for backward compatibility)',
        required=False, 
        default=None
    )
    parser.add_argument(
        '-c', '--cnvPath',
        help='CNV input file path (optional)',
        required=False, default=None
    )
    parser.add_argument(
        '-r', '--repeatPath',
        help='repeats input file path (optional)',
        required=False, default=None
    )
    
    args = parser.parse_args()
    
    # Define input variables.
    output_path = args.outputPath
    sv_paths = args.svPaths if args.svPaths else ([args.svPath] if args.svPath else [])
    cnv_path = args.cnvPath
    repeat_path = args.repeatPath
    skip_svtype = True

    # Check that there is at least one input
    if not sv_paths and not cnv_path and not repeat_path:
        logger.warning("Empty/no vcf files to concatenate")
        return

    struct_lines = {
        "sv": [],
        "cnv": None,
        "repeat": None
    }

    # Process multiple SV files
    if sv_paths:
        for sv_path in sv_paths:
            if sv_path and sv_path != "OPTIONAL_FILE":
                if not os.path.isfile(sv_path):
                    logger.warning(f'SV file "{sv_path}" does not exist.')
                    continue
                
                lines = get_lines(sv_path)
                if is_empty(lines):
                    logger.warning(f'SV file "{sv_path}" is empty.')
                    continue
                
                struct_lines["sv"].append(lines)
    
    # If no valid SV files, set to None
    if not struct_lines["sv"]:
        struct_lines["sv"] = None

    # Process CNV file
    if cnv_path and cnv_path != "OPTIONAL_FILE":
        if not os.path.isfile(cnv_path):
            logger.warning(f'CNV file "{cnv_path}" does not exist.')
        else:
            lines = get_lines(cnv_path)
            if is_empty(lines):
                logger.warning('CNV file is empty.')
            else:
                struct_lines["cnv"] = lines

    # Process repeat file
    if repeat_path and repeat_path != "OPTIONAL_FILE":
        if not os.path.isfile(repeat_path):
            logger.warning(f'Repeat file "{repeat_path}" does not exist.')
        else:
            lines = get_lines(repeat_path)
            if is_empty(lines):
                logger.warning('Repeat file is empty.')
            else:
                struct_lines["repeat"] = lines

    create_unified_file(struct_lines, output_path, skip_svtype)


if __name__ == "__main__":
    main()