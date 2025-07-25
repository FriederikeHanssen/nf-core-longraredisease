#!/usr/bin/env python3
import gzip
import os
import logging
import shutil

# Gets a file and return its lines.
# This is done in order to handle both zipped and unzipped files.
def __get_lines__(file_path: str) -> []:
    try:
        with gzip.open(file_path, 'rt') as file_h:
            lines = file_h.readlines()
    except:
        with open(file_path) as file_h:
            lines = file_h.readlines()
    return lines

# adds SVTYPE=REP to the INFO section of the vcf
# this is necessary for repeat vcf for providers that are not ONT, to enable
# Geneyx application reading these lines as REP effect.
# ONT repeat files already contain SVTYPE=STR in their INFO section, which is
# also recognized by Geneyx as REP effect
def __modify_repeat_line__(repeat_line: str) -> str:
    (chrom, pos, id, ref, alt, qual, filter, info, format, _) = repeat_line.split("\t")
    info = info + ';SVTYPE=REP'
    repeat_line = "\t".join((chrom, pos, id, ref, alt, qual, filter, info, format, _))
    return repeat_line

def __modify_roh_line__(vcf_line) -> str:
        bed_data = vcf_line.strip().split('\t')
        #Bed files are 0 Based while VCF files are 1 based.
        pos_value = min([int(bed_data[1]),int(bed_data[2])]) + 1
        end_value = max([int(bed_data[1]),int(bed_data[2])])
        chromosome = bed_data[0]
        roh_score = bed_data[3]
        return f"{chromosome}\t{pos_value}\t.\tN\t<ROH>\t.\tPASS\tEND={end_value};SVTYPE=ROH;ROH_SCORE={roh_score}\tGT\t1/1\n"


# checks whether this line contains SVTYPE= in it
# this is Geneyx way to recognizing structural variants
# also checks that the gt (genotype) is not "./." (which means a no-call)
# a line is valid if it contains SVTYPE= and it had a valid genotype
def __is_valid_line__(line: str) -> bool:
    # check if line is header or null / empty
    if line.startswith('#') or not line:
        return True

    cells = line.split('\t')
    info_cell = cells[7]
    format_key = cells[8]
    format_value = cells[9]

    # if the "svtype" line is exists the line valid
    if 'SVTYPE=' in info_cell:
        return True

    gt_index = format_key.split(':').index('GT')
    gt_value = format_value.split(':')[gt_index]

    # if "svtype" not present and GT is 0 or "./." the line is invalid
    if gt_index == '0' or gt_value == './.':
        return False

    return True


# writes the content of a file into the output file
# skips the header of the file
# modifies lines when necessary (e.g
def __write_vcf_content__(ftype: str, lines_dict: dict, output_h, skip_svtype: bool):
    start_read_content_flag = False
    if ftype == 'roh':
        start_read_content_flag = True
    for vcf_line in lines_dict[ftype]:
        if start_read_content_flag:
            if ftype == 'repeat' and not skip_svtype:
                vcf_line = __modify_repeat_line__(vcf_line)
            elif ftype == 'roh':
                vcf_line = __modify_roh_line__(vcf_line)
            output_h.write(vcf_line)
        if vcf_line.find('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT') != -1:
            start_read_content_flag = True

    if not start_read_content_flag:
        logging.warning(f'No vcf header for given {ftype} file')


# gets vcf lines and returns the header lines
def __get_vcf_header_from_lines__(vcf_lines: []):
    vcf_header = ""
    #IF ROH
    for line in vcf_lines:
        vcf_header = vcf_header + line
        if line.find("#CHROM") != -1:
            break
    return vcf_header


def __is_empty__(lines: []) -> bool:
    for line in lines:
        striped_line = line.strip()
        if striped_line != '' and not striped_line.startswith('#'):
            return False

    return True

# sorts the vcf so that it will be indexed when uploaded to Geneyx, to allow viewing in the IGV
# the sorting is not done using bcftools but rather with linux sorting of the content, to avoid problems with the header
def __sort_vcf__(vcf_file_path: str):
    # create a directory to put all the temporary files in
    temp_dir_name = "vcf-temp"
    if(not os.path.exists(temp_dir_name)):
        os.mkdir(temp_dir_name)
    try:
        header_file_path = os.path.join(temp_dir_name, "header.txt")
        content_file_path = os.path.join(temp_dir_name, "content.vcf")
        sorted_content_file_path = os.path.join(temp_dir_name, "sorted_content.vcf")
        os.system(f'grep "^#" "{vcf_file_path}" > {header_file_path}')
        os.system(f'grep -v "^#" "{vcf_file_path}" > {content_file_path}')
        os.system(f'sort -k1,1V -k2,2n {content_file_path} > {sorted_content_file_path}')
        os.system(f'cat {header_file_path} {sorted_content_file_path} > "{vcf_file_path}"')
        os.system(f'bgzip "{vcf_file_path}"')
        bgzipped_vcf_path = f'{vcf_file_path}.gz'
    except:
        logging.error("could not sort unified vcf file")
    finally:
        if(os.path.exists(temp_dir_name)):
            shutil.rmtree(temp_dir_name)
        pass


# prints the sv file complete (with header), then the cnv without a header
# and repeats without a header and with modified info section
# if sv file is empty prints all cnv lines (with header) and repeat lines without a header
# if cnv is also empty prints repeat lines with the header

def __create_unified_file__(files_lines: dict, output_path: str, skip_svtype: bool):

    cnv_printed = False
    rep_printed = False

    with open(output_path, 'w+') as output_h:

        # prints into the output file all the sv file (with header) 
        # or all the cnv file (with header) if the sv file is empty 
        # or all the repeats file (with header) if the sv and cnv files are empty

        if files_lines['sv'] is not None:
            output_h.write("".join(files_lines['sv']))
        elif files_lines['cnv'] is not None:
            cnv_printed = True
            output_h.write("".join(filter(lambda l: __is_valid_line__(l), files_lines['cnv'])))
        elif files_lines['repeat'] is not None:
            rep_printed = True
            # separates the header and the content since repeat lines should be modified
            rep_header = __get_vcf_header_from_lines__(files_lines['repeat'])
            output_h.write(rep_header)
            # prints the content of the vcf without the header, and modifying the INFO section
            __write_vcf_content__('repeat', files_lines, output_h, skip_svtype)
        else:
            logging.warning('Empty/no vcf files to concatenate')
            return

        # in case sv exist, adds the CNV lines (without the header)
        # if sv doesn't exist, CNV was already written into the output file with its header
        if files_lines['cnv'] is not None and not cnv_printed:
            __write_vcf_content__('cnv', files_lines, output_h, skip_svtype)
        # in case sv or cnv exist, adds the repeats lines (without the header)
        # if sv and cnv don't exist, the repeats file was already written into the output file with its header
        if files_lines['repeat'] is not None and not rep_printed:
            __write_vcf_content__('repeat', files_lines, output_h, skip_svtype)
        if files_lines['roh'] is not None:
            __write_vcf_content__('roh', files_lines, output_h, skip_svtype)

        output_h.flush()
        output_h.close()
    
    
    # sorts the unified vcf to enable indexing in the application
    __sort_vcf__(output_path)


# creates a dictionary with the different vcf files per type and calls the unifying function
def run(output_path: str, sv_path: str = None, cnv_path: str = None, repeat_path: str = None, roh_bed_path : str = None, skip_svtype: bool = False):
    output_file_path = output_path

    struct_files = {
        'sv': sv_path,
        'cnv': cnv_path,
        'repeat': repeat_path,
        'roh' : roh_bed_path
    }
    struct_lines = {}

    for file_type in struct_files.keys():
        struct_lines[file_type] = None

        if struct_files[file_type] is None:
            logging.warning(f'Can\'t unify file type "{file_type}". The file does not exist.')
            continue

        lines = __get_lines__(struct_files[file_type])
        if __is_empty__(lines):
            logging.warning(f'Can\'t unify file type "{file_type}". The file is empty or contains only headers.')
            continue

        struct_lines[file_type] = lines

    __create_unified_file__(struct_lines, output_file_path, skip_svtype)