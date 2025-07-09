#!/usr/bin/env python3
import gzip
import os
import logging
import shutil
import datetime

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

# Create comprehensive header that includes all field definitions
def __create_comprehensive_header__(sample_name: str = "test") -> str:
    """Create a comprehensive header that includes all necessary definitions from all callers"""
    
    header_lines = [
        "##fileformat=VCFv4.2\n",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n",
        f"##source=UnifiedVCF-MultiCaller\n",
        f"##fileDate={datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n",
        
        # Complete contig definitions for all chromosomes
        "##contig=<ID=chr1>\n", "##contig=<ID=chr2>\n", "##contig=<ID=chr3>\n",
        "##contig=<ID=chr4>\n", "##contig=<ID=chr5>\n", "##contig=<ID=chr6>\n",
        "##contig=<ID=chr7>\n", "##contig=<ID=chr8>\n", "##contig=<ID=chr9>\n",
        "##contig=<ID=chr10>\n", "##contig=<ID=chr11>\n", "##contig=<ID=chr12>\n",
        "##contig=<ID=chr13>\n", "##contig=<ID=chr14>\n", "##contig=<ID=chr15>\n",
        "##contig=<ID=chr16>\n", "##contig=<ID=chr17>\n", "##contig=<ID=chr18>\n",
        "##contig=<ID=chr19>\n", "##contig=<ID=chr20>\n", "##contig=<ID=chr21>\n",
        "##contig=<ID=chr22>\n", "##contig=<ID=chrX>\n", "##contig=<ID=chrY>\n",
        "##contig=<ID=chrMT>\n",
        
        # ALT definitions for all SV types
        "##ALT=<ID=INS,Description=\"Insertion of novel sequence relative to the reference\">\n",
        "##ALT=<ID=DEL,Description=\"Deletion relative to the reference\">\n",
        "##ALT=<ID=DUP,Description=\"Region of elevated copy number relative to the reference\">\n",
        "##ALT=<ID=INV,Description=\"Inversion of reference sequence\">\n",
        "##ALT=<ID=BND,Description=\"Breakend of translocation\">\n",
        "##ALT=<ID=CNV,Description=\"Copy number variant\">\n",
        "##ALT=<ID=LOH,Description=\"Loss of heterozygosity\">\n",
        "##ALT=<ID=REP,Description=\"Tandem repeat\">\n",
        "##ALT=<ID=ROH,Description=\"Run of homozygosity\">\n",
        
        # Comprehensive INFO field definitions
        "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variant\">\n",
        "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variant\">\n",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n",
        "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n",
        "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n",
        "##INFO=<ID=CILEN,Number=2,Type=Integer,Description=\"Confidence interval around inserted/deleted material between breakends\">\n",
        "##INFO=<ID=RE,Number=1,Type=Integer,Description=\"Number of read support this record\">\n",
        "##INFO=<ID=STRAND,Number=A,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format\">\n",
        "##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"Supporting read names of SVs (comma separated)\">\n",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n",
        
        # Sniffles-specific INFO fields
        "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Number of reads supporting the structural variation\">\n",
        "##INFO=<ID=SUPPORT_INLINE,Number=1,Type=Integer,Description=\"Number of reads supporting an INS/DEL SV (non-split events only)\">\n",
        "##INFO=<ID=SUPPORT_LONG,Number=1,Type=Integer,Description=\"Number of soft-clipped reads putatively supporting the long insertion SV\">\n",
        "##INFO=<ID=STDEV_POS,Number=1,Type=Float,Description=\"Standard deviation of structural variation start position\">\n",
        "##INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description=\"Standard deviation of structural variation length\">\n",
        "##INFO=<ID=COVERAGE,Number=.,Type=Float,Description=\"Coverages near upstream, start, center, end, downstream of structural variation\">\n",
        "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count, summed up over all samples\">\n",
        "##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description=\"List of read support for all samples\">\n",
        "##INFO=<ID=CONSENSUS_SUPPORT,Number=1,Type=Integer,Description=\"Number of reads that support the generated insertion (INS) consensus sequence\">\n",
        "##INFO=<ID=NM,Number=.,Type=Float,Description=\"Mean number of query alignment length adjusted mismatches of supporting reads\">\n",
        "##INFO=<ID=PHASE,Number=.,Type=String,Description=\"Phasing information derived from supporting reads\">\n",
        "##INFO=<ID=CUTPASTE,Number=0,Type=Flag,Description=\"Genomic origin of interspersed duplication seems to be deleted\">\n",
        "##INFO=<ID=STD_SPAN,Number=1,Type=Float,Description=\"Standard deviation in span of merged SV signatures\">\n",
        "##INFO=<ID=STD_POS,Number=1,Type=Float,Description=\"Standard deviation in position of merged SV signatures\">\n",
        "##INFO=<ID=STD_POS1,Number=1,Type=Float,Description=\"Standard deviation of breakend 1 position\">\n",
        "##INFO=<ID=STD_POS2,Number=1,Type=Float,Description=\"Standard deviation of breakend 2 position\">\n",
        "##INFO=<ID=MOSAIC,Number=0,Type=Flag,Description=\"Structural variation classified as putative mosaic\">\n",
        
        # Spectre/CNV-specific INFO fields
        "##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">\n",
        
        # ROH-specific INFO fields
        "##INFO=<ID=ROH_SCORE,Number=1,Type=Float,Description=\"Run of homozygosity score\">\n",
        
        # FILTER definitions
        "##FILTER=<ID=q5,Description=\"Quality below 5\">\n",
        "##FILTER=<ID=GT,Description=\"Genotype filter\">\n",
        "##FILTER=<ID=SUPPORT_MIN,Description=\"Minimum read support filter\">\n",
        "##FILTER=<ID=STDEV_POS,Description=\"SV Breakpoint standard deviation filter\">\n",
        "##FILTER=<ID=STDEV_LEN,Description=\"SV length standard deviation filter\">\n",
        "##FILTER=<ID=COV_MIN,Description=\"Minimum coverage filter\">\n",
        "##FILTER=<ID=COV_MIN_GT,Description=\"Minimum coverage filter (missing genotype)\">\n",
        "##FILTER=<ID=COV_CHANGE,Description=\"Coverage change filter\">\n",
        "##FILTER=<ID=COV_CHANGE_INS,Description=\"Coverage change filter for INS\">\n",
        "##FILTER=<ID=COV_CHANGE_FRAC_US,Description=\"Coverage fractional change filter: upstream-start\">\n",
        "##FILTER=<ID=COV_CHANGE_FRAC_SC,Description=\"Coverage fractional change filter: start-center\">\n",
        "##FILTER=<ID=COV_CHANGE_FRAC_CE,Description=\"Coverage fractional change filter: center-end\">\n",
        "##FILTER=<ID=COV_CHANGE_FRAC_ED,Description=\"Coverage fractional change filter: end-downstream\">\n",
        "##FILTER=<ID=MOSAIC_AF,Description=\"Mosaic variant allele frequency filter\">\n",
        "##FILTER=<ID=NOT_MOSAIC_AF,Description=\"Variant allele frequency filter for non-mosaic\">\n",
        "##FILTER=<ID=ALN_NM,Description=\"Length adjusted mismatch filter\">\n",
        "##FILTER=<ID=STRAND_BND,Description=\"Strand support filter for BNDs\">\n",
        "##FILTER=<ID=STRAND,Description=\"Strand support filter for germline SVs\">\n",
        "##FILTER=<ID=STRAND_MOSAIC,Description=\"Strand support filter for mosaic SVs\">\n",
        "##FILTER=<ID=SVLEN_MIN,Description=\"SV length filter\">\n",
        "##FILTER=<ID=SVLEN_MIN_MOSAIC,Description=\"SV length filter for mosaic SVs\">\n",
        "##FILTER=<ID=hom_ref,Description=\"Genotype is homozygous reference\">\n",
        "##FILTER=<ID=not_fully_covered,Description=\"Tandem duplication is not fully covered by a single read\">\n",
        
        # Comprehensive FORMAT field definitions
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
        "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"High-quality reference reads\">\n",
        "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"High-quality variant reads\">\n",
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods rounded to the closest integer\">\n",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n",
        "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase-block, zero if none or not phased\">\n",
        "##FORMAT=<ID=ID,Number=1,Type=String,Description=\"Individual sample SV ID for multi-sample output\">\n",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n",
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth for each allele\">\n",
        
        # Spectre/CNV-specific FORMAT fields
        "##FORMAT=<ID=HO,Number=1,Type=Float,Description=\"Homozygosity\">\n",
        "##FORMAT=<ID=CD,Number=1,Type=Float,Description=\"Copy depth\">\n",
        
        # Column header
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
    ]
    
    return "".join(header_lines)

# writes the content of a file into the output file
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

# Fixed sorting function using bcftools for proper chromosome ordering
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

def __create_unified_file__(files_lines: dict, output_path: str, skip_svtype: bool):
    with open(output_path, 'w+') as output_h:
        
        # Write comprehensive header that includes all field definitions
        comprehensive_header = __create_comprehensive_header__()
        output_h.write(comprehensive_header)
        
        # Write content from each file (without headers)
        for file_type in ['sv', 'cnv', 'repeat', 'roh']:
            if files_lines[file_type] is not None:
                __write_vcf_content__(file_type, files_lines, output_h, skip_svtype)

        output_h.flush()
        output_h.close()
    
    # sorts the unified vcf to enable indexing in the application
    __sort_vcf__(output_path)

# creates a dictionary with the different vcf files per type and calls the unifying function
def run(output_path: str, sv_path: str = None, cnv_path: str = None, repeat_path: str = None, roh_bed_path: str = None, skip_svtype: bool = False):
    output_file_path = output_path

    struct_files = {
        'sv': sv_path,
        'cnv': cnv_path,
        'repeat': repeat_path,
        'roh': roh_bed_path
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