description = '''
A local rephaser build to deal with difficult to phase regions with segmental duplications, such as the region around the CYP21A2 gene.
'''
####################################################################################################
# LOAD PACKAGES                                                                                    #
####################################################################################################

import pysam
import argparse
import sys, os
import Levenshtein
import subprocess

####################################################################################################
# ARGPARSER ARGUMENTS                                                                              #
####################################################################################################

parser = argparse.ArgumentParser(prog = "NanoCAH",
                                 description = description,
                                 epilog = "Good luck with using this software! Best wishes, Laura Skak")

parser.add_argument("-i","--infile",
                    type=str,
                    required = True,
                    help = "The full or relative path to the unfiltered alignment file."
                    )

parser.add_argument("-r","--reference",
                    metavar= "REFERENCE",
                    type=str,
                    required = True,
                    help = "The full or relative path to the reference used for the unfiltered bam file.")

parser.add_argument("-s","--range",
                    metavar= "CHROM:START-END",
                    type=str,
                    help = "The range to be evaluated. Should be in the following format; CHROM:START-END.")

parser.add_argument("-b","--bed",
                    type=str,
                    help = "The full or relative path to a file in bed format containing the ranges to be evaluated.") # FIXME: not implemented

parser.add_argument("-o","--outfile",
                    type=str,
                    default = "outfile.bam",
                    help = "The desired name for alignment outfile.")

parser.add_argument("-d","--outdir",
                    type=str,
                    default = "test", # FIXME: not implemented
                    help = "The full or relative path to the alignment file containing the filtered and phased reads.")

parser.add_argument("-p","--prefix",
                    type=str,
                    default = "",
                    help = "The prefix that should be added to all outfiles.")

parser.add_argument("-a","--assembly",
                    metavar= "FASTA",
                    type=str,
                    default = "asm.fasta",
                    help = "The desired full or relative path to the output file containing the haplotype assemblies.")

parser.add_argument("--max_diff",
                    type=int,
                    default = 10,
                    help = "The maximum allowed SNP differences between blocks. Is used for low variance merging.")

parser.add_argument("--minimum_count",
                    type=int,
                    default = 3,
                    help = "The minimal read count for a SNP to be determined as a usable variant.")

parser.add_argument("--minimum_overlap",
                    type=int,
                    default = 3,
                    help = "The minimal overlapping variants for read-to-block or block-to-block comparison.")

parser.add_argument("--minimum_coverage",
                    type=int,
                    default = 10,
                    help = "The minimal number of reads for a block to be considered.")

parser.add_argument("--minimum_depth",
                    type=int,
                    default = 4,
                    help = "The minimal depth per site to be considered high coverage.")

parser.add_argument("--minimum_variance",
                    type=int,
                    default = 0.50,
                    help = "The maximal allowed variance in a usable SNP.")

parser.add_argument("--minimum_identity",
                    type=int,
                    default = 0.25,
                    help = "The minimum identity required to merge two high confidence blocks.")

parser.add_argument("--maximal_iterations",
                    type=int,
                    default = 5,
                    help = "The maximal allowed variance in a usable SNP.")

parser.add_argument("--maximal_block_overlap",
                    type=int,
                    default = 1000000,
                    help = "The maximal allowed overlap between two blocks to be merged.")

parser.add_argument("--minimum_pct_overlap",
                    type=int,
                    default = 0.5,
                    help = "The minimum percentage overlap for two blocks to merge or a read to be added to a block.")

parser.add_argument("--show_removed_reads",
                    action='store_true',
                    help = "Keep all removed reads in the final bam.")

parser.add_argument("--merge_overlap",
                    action='store_true',
                    help = "If two blocks are found to contain a large insertion, they will be merged to one read. The default is merging blocks with large deletions, but seperating blocks with large insertions in a left and right haplotype. Not implemented yet.") # FIXME: is not implemented

parser.add_argument("--output_raw_assemblies",
                    action='store_true',
                    help = "") # FIXME: is not implemented

parser.add_argument("--output_intermediate_states",
                    action='store_true',
                    help = "This is primarily used for testing, but can sometimes be informative to understand the process of NanoCAH read filtering. It will produce five intermediate alignment outputs; after initial grouping, after no-variance merging, after low variance merging, after read selection and after haplotype creation.") # FIXME: is not implemented

parser.add_argument("-v", "--verbosity",
                    type = int,
                    default = 2,
                    help = "The level of verbosity. Level 0 to 4.")

args = parser.parse_args()


####################################################################################################
# SET CONSTANTS                                                                                    #
####################################################################################################

infile = args.infile
# bam_file = args.outfile
# asm_file = args.assembly

if not os.path.exists(args.outdir):

    os.makedirs(args.outdir)

if args.output_intermediate_states and not os.path.exists(f'{args.outdir}/test_outputs'):

    os.makedirs(f'{args.outdir}/test_outputs')

samfile = pysam.AlignmentFile(infile, "rb")
fastafile = pysam.FastaFile(args.reference)

use_SAs = False # FIXME: Virker ikke, når True
minimum_base_quality = 0 # FIXME: måske ikke nødvendig
minimum_mapping_quality = 40 # FIXME: not implemented
minimum_count = args.minimum_count
maxdiff = args.max_diff
minimum_overlap = args.minimum_overlap
minimum_depth = args.minimum_depth
minimum_coverage = args.minimum_coverage
max_cleanup_iterations = args.maximal_iterations
max_variance = args.minimum_variance
minimum_identity = args.minimum_identity
minimum_pct_overlap = args.minimum_pct_overlap
expected_polyploidy = 2 # FIXME: MISSING ARG. Not implemented
show_removed = args.show_removed_reads
verbosity = args.verbosity
min_overlap = 1000 # FIXME: MISSING ARG.
# max_overlap = 1000000 # FIXME: MISSING ARG.
overlap_jump = 500 # FIXME: MISSING ARG.
min_overlap_identity = 0.75
max_error_length = 10000 # FIXME: MISSING ARG.
max_variance = 0.30 # FIXME: MISSING ARG.
minimum_pct_common_reads = 0.05 # FIXME: MISSING ARG.

seperator = "******************************************************************"

####################################################################################################
# FUNCTIONS                                                                                        #
####################################################################################################

def create_range_list():

    if args.bed != None:

        if args.range != None:

            print(f'Both a range and a bed file is provided. Only the ranges in bedfile will be evaulated.', file=sys.stderr)

        range_list = []
        
        with open(args.bed, 'r') as file:

            for line in file.readlines():

                row = line.strip().split("\t")

                row[1] = int(row[1])
                row[2] = int(row[2])

                range_list.append(tuple(row))

    elif args.range != None:

        chrom = args.range.split(":")[0]
        start = int(args.range.split(":")[1].split("-")[0])
        end = int(args.range.split(":")[1].split("-")[1])

        range_list = [(chrom, start, end, "range")]

    else:

        print(f'Neither a range and a bed file is provided. Use the -b or -s flags and run again.', file=sys.stderr)
    
    return range_list

def print_subprocess_title(text):
    
    if verbosity > 2:

        print(f"\n{seperator}\n{text}\n{seperator}\n", file=sys.stderr)

    elif verbosity > 0:

        print(f"\n\n{text}\n", file=sys.stderr)

def get_reference_sequence(read_type, read, start_pos, read_dict, keep_insertions = False):

    if keep_insertions:

        sequence = list(read_dict["p"][read]["-"][0].query_sequence)

    else:

        sequence = read_dict["p"][read]["-"][0].query_sequence

    cigar_tuples = read_dict[read_type][read][start_pos][0].cigartuples

    primary_left_flank = read_dict["p"][read]["-"][0].cigartuples[0]
    primary_right_flank = read_dict["p"][read]["-"][0].cigartuples[-1]

    if keep_insertions:

        if primary_left_flank[0] == 5:

            sequence = ["N"]*primary_left_flank[1]+sequence

        if primary_right_flank[0] == 5:

            sequence = sequence+["N"]*primary_right_flank[1]

        pass

    else:

        if primary_left_flank[0] == 5:

            sequence = "N"*primary_left_flank[1]+sequence

        if primary_right_flank[0] == 5:

            sequence = sequence+"N"*primary_right_flank[1]

    running_pos = 0

    for tuple in cigar_tuples:

        if tuple[0] == 4:

            if running_pos == 0:

                sequence = sequence[int(tuple[1]):]

            else:

                sequence = sequence[:running_pos]

            running_pos -= int(tuple[1])

        elif tuple[0] == 2:

            if keep_insertions:

                sequence = sequence[:running_pos] + ["D"]*int(tuple[1]) + sequence[running_pos:]

            else:

                sequence = sequence[:running_pos] + "D"*int(tuple[1]) + sequence[running_pos:]

        elif tuple[0] == 1:

            if keep_insertions:

                sequence = sequence[:running_pos-1] + ["".join(sequence[running_pos-1:running_pos-1+tuple[1]+1])] + sequence[running_pos+tuple[1]:]

            else:

                sequence = sequence[:running_pos-1] + "I" + sequence[running_pos+tuple[1]:]

            running_pos -= int(tuple[1])

        running_pos += int(tuple[1])

    return sequence

def add_secondary_reads_to_read_dict(read_dict, variant_dict): # Calls get_reference_sequence

    unmatched_secondary_dict = {}

    for read in list(read_dict["s"].keys()):

        if read not in read_dict["p"]: # FIXME: Deleted secondary reads with no primary partner, which is bad

            if read not in unmatched_secondary_dict:

                unmatched_secondary_dict[read] = {}

            for start_pos in read_dict["s"][read]:

                read_dict["s"][read][start_pos][3] = "unmatched_secondary"

                unmatched_secondary_dict[read][start_pos] = read_dict["s"][read][start_pos]

            del read_dict["s"][read]

            continue

        for start_pos in read_dict["s"][read]:

            chrom = read_dict["s"][read][start_pos][0].reference_name

            sequence = get_reference_sequence("s", read, start_pos, read_dict)

            if read_dict["s"][read][start_pos][0].reference_end-read_dict["s"][read][start_pos][0].reference_start != len(sequence): # FIXME: test, hvis den kører er der noget galt

                print("If this prints, there is something wrong with add_secondary_reads_to_read_dict.")

            read_variant_dict = {}

            for i in range(len(sequence)):

                base = sequence[i]

                if base == "N":

                    continue

                ref_pos = start_pos+i

                if ref_pos in variant_dict[chrom]:

                    if base not in variant_dict[chrom][ref_pos]:

                        variant_dict[chrom][ref_pos][base] = []

                    variant_dict[chrom][ref_pos][base].append(("s", read, start_pos))

                    read_variant_dict[ref_pos] = base

                    continue

                if base in ["D", "I"]:

                    continue

                ref_base = fastafile.fetch(read_dict["s"][read][start_pos][0].reference_name, ref_pos, ref_pos+1).upper()

                if base != ref_base:

                    if ref_pos not in variant_dict[chrom]:

                        variant_dict[chrom][ref_pos] = {}

                    if base not in variant_dict[chrom][ref_pos]:

                        variant_dict[chrom][ref_pos][base] = []

                    variant_dict[chrom][ref_pos][base].append(("s", read, start_pos))

                    read_variant_dict[ref_pos] = base

            read_dict["s"][read][start_pos][1] = read_variant_dict
    
    return unmatched_secondary_dict

def add_read_to_base_dict(base_dict, key, read_type, value):

    if key not in base_dict:

        base_dict[key] = []

    base_dict[key].append((read_type, value, "-"))

def add_primary_reads_to_read_dict(read_dict, variant_dict): # Calls add_read_to_base_dict

    range_list = create_range_list()

    for chrom, start, end, name in range_list:

        iter = samfile.pileup(chrom, start, end)

        for pileupcolumn in iter:

            pos = pileupcolumn.pos

            if pos in variant_dict[chrom]:

                ref_base = fastafile.fetch(chrom, pos, pos+1).upper()

                base_dict = variant_dict[chrom][pos]

                for pileupread in pileupcolumn.pileups:

                    if pileupread.alignment.has_tag("SA") and not use_SAs: # FIXME: Don't know how well it handles supplementary alignments

                        continue

                    if pileupread.alignment.is_secondary:

                        continue

                    if pileupread.is_del:

                        add_read_to_base_dict(base_dict, "D", "p", pileupread.alignment.query_name)

                        read_dict["p"][pileupread.alignment.query_name]["-"][1][pos] = "D"

                    elif pileupread.indel > 0:

                        add_read_to_base_dict(base_dict, "I", "p", pileupread.alignment.query_name)

                        read_dict["p"][pileupread.alignment.query_name]["-"][1][pos] = "I"

                    elif ref_base != pileupread.alignment.query_sequence[pileupread.query_position]:

                        alt_base = pileupread.alignment.query_sequence[pileupread.query_position]

                        add_read_to_base_dict(base_dict, alt_base, "p", pileupread.alignment.query_name)

                        read_dict["p"][pileupread.alignment.query_name]["-"][1][pos] = alt_base

                    else:

                        add_read_to_base_dict(base_dict, ref_base, "p", pileupread.alignment.query_name)
                        read_dict["p"][pileupread.alignment.query_name]["-"][1][pos] = ref_base

def add_supplementary_reads_to_supplementary_read_dict(read_dict, variant_dict):

    range_list = create_range_list()

    for chrom, start, end, name in range_list:
            
        iter = samfile.pileup(chrom, start, end)

        for pileupcolumn in iter:

            pos = pileupcolumn.pos

            if pos in variant_dict:

                ref_base = fastafile.fetch(chrom, pos, pos+1).upper()

                for pileupread in pileupcolumn.pileups:

                    if not pileupread.alignment.has_tag("SA"):

                        continue

                    if pileupread.is_del:

                        read_dict[pileupread.alignment.query_name][pileupread.alignment.reference_start][1][pos] = "D"

                    elif pileupread.indel > 0:

                        read_dict[pileupread.alignment.query_name][pileupread.alignment.reference_start][1][pos] = "I"

                    elif ref_base != pileupread.alignment.query_sequence[pileupread.query_position]:

                        alt_base = pileupread.alignment.query_sequence[pileupread.query_position]

                        read_dict[pileupread.alignment.query_name][pileupread.alignment.reference_start][1][pos] = alt_base

                    else:

                        read_dict[pileupread.alignment.query_name][pileupread.alignment.reference_start][1][pos] = ref_base

def create_dicts(samfile): # Calls add_secondary_reads_to_read_dict and add_primary_reads_to_read_dict

    print_subprocess_title("Making the read dict and variant dict.")

    # Variant dict will be a dictionary containing all positions (keys) that hav at least one read that differs from the reference. The values are the the pileupcolumn object.
    variant_dict = {}

    # The read dict will contain all reads that overlap with any of the positions in the variant dict. The values are the pileupread.alignment object.
    read_dict = {"p": {}, "s": {}}

    # The read dict will contain all reads that overlap with any of the positions in the variant dict. The values are the pileupread.alignment object.
    supplementary_read_dict = {}

    # Keeps track of the depth to calculate if a SNP is over- or under covered compared to the median.
    depth_dict = {}

    range_list = create_range_list()

    for chrom, start, end, name in range_list:

        if chrom not in variant_dict:

            variant_dict[chrom] = {}
            depth_dict[chrom] = {}

        iter = samfile.pileup(chrom, start, end, stepper="nofilter")

        for pileupcolumn in iter:

            pileupcolumn.set_min_base_quality(0)

            pos = pileupcolumn.pos

            depth_dict[chrom][pos] = pileupcolumn.n

            # Checks if position is within range
            if not (start <= pos and pos <= end):

                continue

            ref_base = fastafile.fetch(chrom, pos, pos+1).upper()

            for pileupread in pileupcolumn.pileups:

                if pileupread.alignment.has_tag("SA"): # FIXME: Don't know how well it handles supplementary alignments

                    if pileupread.alignment.query_name not in supplementary_read_dict:

                        supplementary_read_dict[pileupread.alignment.query_name] = {}

                    if pileupread.alignment.reference_start not in supplementary_read_dict[pileupread.alignment.query_name]:

                        supplementary_read_dict[pileupread.alignment.query_name][pileupread.alignment.reference_start] = [pileupread.alignment, {}, {}, "supplementary", name]

                    continue

                if pileupread.alignment.is_secondary:

                    if pileupread.alignment.query_name not in read_dict["s"]:

                        read_dict["s"][pileupread.alignment.query_name] = {}

                    if pileupread.alignment.reference_start not in read_dict["s"][pileupread.alignment.query_name]:

                        read_dict["s"][pileupread.alignment.query_name][pileupread.alignment.reference_start] = [pileupread.alignment, {}, {}, "ungrouped", name]

                    continue

                if pileupread.alignment.query_name not in read_dict["p"]:

                    read_dict["p"][pileupread.alignment.query_name] = {}

                    read_dict["p"][pileupread.alignment.query_name]["-"] = [pileupread.alignment, {}, {}, "ungrouped", name]

                if pileupread.is_refskip: # FIXME: test

                    print("If this prints, there is a problem. A read is refskip and handling of this is not implemented yet.", pileupread.is_refskip)

                    continue

                if pileupread.query_position == None:

                    continue

                if ref_base != pileupread.alignment.query_sequence[pileupread.query_position] and pileupread.alignment.query_qualities[pileupread.query_position] >= minimum_base_quality:

                    variant_dict[chrom][pos] = {}

    unmatched_secondary_dict = add_secondary_reads_to_read_dict(read_dict, variant_dict)

    add_primary_reads_to_read_dict(read_dict, variant_dict)

    add_supplementary_reads_to_supplementary_read_dict(supplementary_read_dict, variant_dict)

    for chrom in variant_dict:

        variant_dict[chrom] = dict(sorted(variant_dict[chrom].items()))

    for read_type in read_dict:
        for read in read_dict[read_type]:
            for start_pos in read_dict[read_type][read]:

                read_dict[read_type][read][start_pos][2] = read_dict[read_type][read][start_pos][1].copy()

    if verbosity > 1:

        print(f'{sum([len(variant_dict[chrom]) for chrom in variant_dict])} variable positions were found.\n', file=sys.stderr)

    return variant_dict, read_dict, depth_dict, supplementary_read_dict, unmatched_secondary_dict


def remove_variant(chrom, pos, haploblock_dict, read_dict):

    for base in variant_dict[chrom][pos]:

        for read_type, read, start_pos in variant_dict[chrom][pos][base]:

            if pos in read_dict[read_type][read][start_pos][1]:

                del read_dict[read_type][read][start_pos][1][pos]

            if pos in read_dict[read_type][read][start_pos][2]:

                del read_dict[read_type][read][start_pos][2][pos]

    del variant_dict[chrom][pos]

    for hap in haploblock_dict:

        if pos in haploblock_dict[hap][0]:

            del haploblock_dict[hap][0][pos]

def remove_bad_variants(variant_dict): # Calls remove_variant

    print_subprocess_title("Removing bad variants.")

    for chrom in variant_dict:
        for pos in list(variant_dict[chrom].keys()):

            depth_dict[chrom][pos] = sum([len(variant_dict[chrom][pos][base]) for base in variant_dict[chrom][pos]]) # FIXME: new addition, because there was a problem with removing variant because of high supplementary read counts. 

            ref_base = fastafile.fetch(chrom, pos, pos+1).upper()

            if len(variant_dict[chrom][pos]) == 1:

                remove_variant(chrom, pos, {}, read_dict)

                continue

            if sum(count >= minimum_count for count in [len(variant_dict[chrom][pos][base]) for base in variant_dict[chrom][pos]]) <= 1:

                remove_variant(chrom, pos, {}, read_dict)

                continue

            if len(variant_dict[chrom][pos]) > 2:

                for base in list(variant_dict[chrom][pos].keys()):

                    if len(variant_dict[chrom][pos][base]) < minimum_count:

                        for read_type, read, start_pos in variant_dict[chrom][pos][base]:

                            if pos in read_dict[read_type][read][start_pos][1]:

                                del read_dict[read_type][read][start_pos][1][pos]
                                del read_dict[read_type][read][start_pos][2][pos]

                        del variant_dict[chrom][pos][base]

                if len(variant_dict[chrom][pos]) < 2:

                    remove_variant(chrom, pos, {}, read_dict)

                    continue

            # The following section catches if the remaining variants are only the reference base and a insertion or deletion. These variants are not informative.

            if sum([base in ["I", "D", ref_base] for base in variant_dict[chrom][pos]]) == len(variant_dict[chrom][pos]):

                remove_variant(chrom, pos, {}, read_dict)

                continue

            if sum([len(variant_dict[chrom][pos][base]) if base != "deleted" else 0 for base in variant_dict[chrom][pos]])/depth_dict[chrom][pos] < 1-max_variance:

                remove_variant(chrom, pos, {}, read_dict)

                continue

    if verbosity > 1:

        print(f'{sum([len(variant_dict[chrom]) for chrom in variant_dict])} usable variants were found.\n', file=sys.stderr)

    return variant_dict


def remove_pos_from_haploblocks(read_type, read, start_pos, pos, haploblock_dict):

    # Makes sure that the variant is removed from the haploblock, if this read is the only one supporting it.

    block = read_dict[read_type][read][start_pos][3]

    if block in haploblock_dict and pos in haploblock_dict[block][0]:

        if sum([pos in read_dict[read_type][read][start_pos][1] for read_type, read, start_pos in haploblock_dict[block][1]]) == 0:

            del haploblock_dict[block][0][pos]

            if verbosity > 3:

                print(f'\tpos {pos} removed from {block}.', file=sys.stderr)

def remove_read_from_variant_dict(pos, read_type, read, start_pos, haploblock_dict): # Calls remove_pos_from_haploblocks and remove_variant

    base = read_dict[read_type][read][start_pos][1][pos]
    chrom = read_dict[read_type][read][start_pos][0].reference_name

    for i in range(len(variant_dict[chrom][pos][base])):

        if variant_dict[chrom][pos][base][i] == (read_type, read, start_pos):

            del variant_dict[chrom][pos][base][i]

            # To make sure the read variants are still deleted if the variant or variant base is later removed they are saved in a "deleted" category

            if "deleted" not in variant_dict[chrom][pos]:

                variant_dict[chrom][pos]["deleted"] = []

            variant_dict[chrom][pos]["deleted"].append((read_type, read, start_pos))

            break

    # If the base count is less than the minimum count the base is removed from variant dict and reads in read dict.

    if len(variant_dict[chrom][pos][base]) < minimum_count:

        for read_type_1, read_1, start_pos_1 in variant_dict[chrom][pos][base]:

            del read_dict[read_type_1][read_1][start_pos_1][1][pos]
            del read_dict[read_type_1][read_1][start_pos_1][2][pos]

            remove_pos_from_haploblocks(read_type_1, read_1, start_pos_1, pos, haploblock_dict)

        del variant_dict[chrom][pos][base]

        new_deleted_list = []

        for read_type_1, read_1, start_pos_1 in variant_dict[chrom][pos]["deleted"]:

            if read_dict[read_type_1][read_1][start_pos_1][2][pos] == base:

                del read_dict[read_type_1][read_1][start_pos_1][2][pos]

                remove_pos_from_haploblocks(read_type_1, read_1, start_pos_1, pos, haploblock_dict)

            else:

                new_deleted_list.append((read_type_1, read_1, start_pos_1))

        if len(new_deleted_list) == 0:

            del variant_dict[chrom][pos]["deleted"]

        else:

            variant_dict[chrom][pos]["deleted"] = new_deleted_list

        for block in haploblock_dict:

            if pos in haploblock_dict[block][0] and haploblock_dict[block][0][pos] == base:

                del haploblock_dict[block][0][pos]

                if verbosity > 3:

                    print(f'\tpos {pos} removed from {block}', file=sys.stderr)

    if sum([base in ["deleted", "I", "D", fastafile.fetch(chrom, pos, pos+1).upper()] for base in variant_dict[chrom][pos]]) == len(variant_dict[chrom][pos]) or sum([base not in ["deleted"] for base in variant_dict[chrom][pos]]) == 1:

        remove_variant(chrom, pos, haploblock_dict, read_dict)

    elif "deleted" in variant_dict[chrom][pos] and len(variant_dict[chrom][pos]["deleted"])/sum([len(variant_dict[chrom][pos][base]) for base in variant_dict[chrom][pos]]) > max_variance: # If the deleted base count exceeds max_variance it is likely a bad variant and should therefore be removed.

        remove_variant(chrom, pos, haploblock_dict, read_dict)
    
    elif sum([len(variant_dict[chrom][pos][base]) if base != "deleted" else 0 for base in variant_dict[chrom][pos]])/depth_dict[chrom][pos] < 1-max_variance:

        remove_variant(chrom, pos, haploblock_dict, read_dict)

    if pos in read_dict[read_type][read][start_pos][1]:

        del read_dict[read_type][read][start_pos][1][pos]

        remove_pos_from_haploblocks(read_type, read, start_pos, pos, haploblock_dict)

def make_fit_list(other_variant_dict, chrom, range, haploblock_dict, max_diff, minimum_overlap):

    fit_list = []

    if len(other_variant_dict) == 0:

        return fit_list 

    for hap in haploblock_dict:

        haploblock_variant_dict = haploblock_dict[hap][0]

        if len(haploblock_dict[hap][0]) == 0:

            continue

        # if chrom != haploblock_dict[hap][2]:

        #     continue

        if range != haploblock_dict[hap][3]:

            continue

        if not (min(haploblock_variant_dict) <= max(other_variant_dict) and min(other_variant_dict) <= max(haploblock_variant_dict)):
        # if min(list(haploblock_variant_dict.keys())) <= max(list(other_variant_dict.keys())) and min(list(other_variant_dict.keys())) <= max(list(haploblock_variant_dict.keys())):

            continue

        overlap = 0
        identity = 0
        mismatch_list = []

        for pos in list(haploblock_variant_dict.keys()):

            if pos not in variant_dict[chrom]: # FIXME: test

                print(f'If this prints there is a problem. pos {pos} not in variant dict. {hap}.')

            if pos in other_variant_dict:

                overlap += 1

                if haploblock_variant_dict[pos] == other_variant_dict[pos]:

                    identity += 1

                else:

                    mismatch_list.append((pos, haploblock_variant_dict[pos], other_variant_dict[pos]))
            
                if len(mismatch_list) > max_diff: # This saves time by termination the forloop when the SNP difference is more than the maximal SNP difference

                    break

        if overlap == 0:

            continue

        if overlap >= minimum_overlap+len(mismatch_list) and len(mismatch_list) <= max_diff:

            fit_list.append((hap, overlap, identity, mismatch_list))

    return fit_list

def add_reads_to_haplotype_dict(): # Calls make_fit_list and remove_read_from_variant_dict

    print_subprocess_title("Adding reads to haploblock dict.")

    haploblock_dict = {}

    hap_count = 1

    prev_reads = []

    for chrom in variant_dict:
        for variant_pos in list(variant_dict[chrom].keys()):

            if variant_pos not in variant_dict[chrom]:

                continue

            base_dict = variant_dict[chrom][variant_pos]

            for base in base_dict:

                read_list = base_dict[base]

                for read_type, read, start_pos in read_list:

                    if f'{read}_{read_type}_{start_pos}' in prev_reads:

                        continue

                    else:

                        prev_reads.append(f'{read}_{read_type}_{start_pos}')

                    if len(read_dict[read_type][read][start_pos][2]) == 0:

                        continue

                    read_variant_dict = read_dict[read_type][read][start_pos][2]

                    fit_list = make_fit_list(read_variant_dict, read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, 0, minimum_overlap)

                    if len(fit_list) == 0:

                        hap_name = f'block_{hap_count}'

                        haploblock_dict[hap_name] = [read_variant_dict.copy(), [(read_type, read, start_pos)], read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4]]

                        hap_count += 1

                        read_dict[read_type][read][start_pos][3] = hap_name

                    elif len(fit_list) == 1:

                        hap_name = fit_list[0][0]

                        haploblock_dict[hap_name][0].update(read_variant_dict.copy())

                        haploblock_dict[hap_name][1].append((read_type, read, start_pos))

                        # for pos, _, _ in fit_list[0][3]: # FIXME: should not be necessary

                        #     remove_read_from_variant_dict(pos, read_type, read, start_pos, haploblock_dict)

                        read_dict[read_type][read][start_pos][3] = hap_name

                    else:

                        max_overlap = 0
                        max_length = 0
                        to_hap = None

                        # Finds the haplotype with the largest overlap.
                        for hap, overlap, _, _ in fit_list:

                            if overlap > max_overlap:

                                max_overlap = overlap
                                max_length = len(haploblock_dict[hap][1])
                                to_hap = hap

                            elif overlap == max_overlap and len(haploblock_dict[hap][1]) > max_length:

                                max_length = len(haploblock_dict[hap][1])
                                to_hap = hap

                        haploblock_dict[to_hap][0].update(read_variant_dict.copy())
                        haploblock_dict[to_hap][1].append((read_type, read, start_pos))

                        read_dict[read_type][read][start_pos][3] = to_hap

    # Makes the ungrouped list

    not_unique_list = []

    for read_type in read_dict:
        for read in read_dict[read_type]:
            for start_pos in read_dict[read_type][read]:

                if len(read_dict[read_type][read][start_pos][2]) == 0:

                    read_dict[read_type][read][start_pos][3] = "unphasable"

                    continue

                if read_dict[read_type][read][start_pos][3] == "ungrouped":

                    not_unique_list.append((read_type, read, start_pos))

    if verbosity > 1:

        print(f'{len(haploblock_dict)} blocks were created during initial grouping.\n', file=sys.stderr)

    if args.output_intermediate_states:

        make_test_output(f'{args.outdir}/test_outputs/{args.prefix}_secondary_filtered_1_of_5_post_initial_grouping.bam')

    return haploblock_dict, not_unique_list, hap_count


def reintroduce_read_to_variant_dict(read_type, read, start_pos):

    chrom = read_dict[read_type][read][start_pos][0].reference_name

    read_dict[read_type][read][start_pos][3] = "ungrouped"

    for pos in read_dict[read_type][read][start_pos][2]:

        if "deleted" in variant_dict[chrom][pos]:

            for i in range(len(variant_dict[chrom][pos]["deleted"])):

                read_vec = variant_dict[chrom][pos]["deleted"][i]

                if read_vec == (read_type, read, start_pos):

                    del variant_dict[chrom][pos]["deleted"][i]

                    if len(variant_dict[chrom][pos]["deleted"]) == 0:

                        del variant_dict[chrom][pos]["deleted"]

                    base = read_dict[read_type][read][start_pos][2][pos]

                    variant_dict[chrom][pos][base].append((read_type, read, start_pos))

                    break

    read_dict[read_type][read][start_pos][1] = read_dict[read_type][read][start_pos][2].copy()

def split_blocks(minimum_depth, haploblock_dict, not_unique_list): # Calls reintroduce_read_to_variant_dict

    if verbosity > 3:

        print(f'Blocks are split based on a required minimum depth of {minimum_depth}.', file=sys.stderr)


    for block in list(haploblock_dict.keys()):

        blocks = []

        start = None
        end = None
        first = True

        if len(haploblock_dict[block][1]) < minimum_depth:

            if verbosity > 3:

                print(f'{block} is removed.', file=sys.stderr)

            for read_type, read, start_pos in haploblock_dict[block][1]:

                reintroduce_read_to_variant_dict(read_type, read, start_pos)

                not_unique_list.append((read_type, read, start_pos))

            del haploblock_dict[block]

            continue

        if len(haploblock_dict[block][1]) == 1: # There's no splitting one read, so might as well continue

            continue

        pos_list = sorted(list(haploblock_dict[block][0].keys()))

        for i in range(len(pos_list)):

            pos = pos_list[i]

            new_reads = []

            for read_type, read, start_pos in haploblock_dict[block][1]:

                if pos in read_dict[read_type][read][start_pos][1]:

                    new_reads.append(f'{read_type}_{read}_{start_pos}')

            if first and len(new_reads) >= minimum_depth:

                start = pos
                end = pos

                all_reads = new_reads

                first = False

                continue

            elif first and len(new_reads) < minimum_depth:

                continue

            found_count = sum(read in all_reads for read in new_reads)

            if found_count >= minimum_depth:

                end = pos

                for read in new_reads:

                    if read not in all_reads:

                        all_reads.append(read)

            elif start != None and end != None:

                blocks.append([start, end, {}, [], haploblock_dict[block][2], haploblock_dict[block][3], all_reads])

                if len(new_reads) >= minimum_depth:

                    all_reads = new_reads
                    start = pos
                    end = pos

                else:

                    start = None
                    end = None
                    first = True

        if start != None and end != None:

            blocks.append([start, end, {}, [], haploblock_dict[block][2], haploblock_dict[block][3], all_reads])

        print(f'   {block}', [(block[0], block[1], len(block[2]), len(block[3])) for block in blocks], len(haploblock_dict[block][1]), len(haploblock_dict[block][0]), file=sys.stderr)

        if len(blocks) > 1:

            new_blocks = []

            i = 0

            while i < len(blocks)-1:

                found = False

                for j in range(i+1, len(blocks)):

                    if sum([read in blocks[j][6] for read in blocks[i][6]]) >= minimum_depth:

                        all_reads = []

                        for k in range(i, j+1):

                            all_reads += blocks[k][6]

                        blocks[j] = [blocks[i][0], blocks[j][1], {}, [], haploblock_dict[block][2], haploblock_dict[block][3], list(set(all_reads))]

                        found = True

                        break

                if not found:

                    new_blocks.append(blocks[i])
                    i += 1

                else:

                    i = j

            new_blocks.append(blocks[-1])

            blocks = new_blocks

        print(f'   {block}', [(block[0], block[1], len(block[2]), len(block[3])) for block in blocks], len(haploblock_dict[block][1]), len(haploblock_dict[block][0]), file=sys.stderr)

        for read_type, read, start_pos in haploblock_dict[block][1]:

            if len(read_dict[read_type][read][start_pos][2]) == 0:

                read_dict[read_type][read][start_pos][3] = "unphasable"

                continue

            read_dict[read_type][read][start_pos][3] = "ungrouped"

            if len(read_dict[read_type][read][start_pos][1]) == 0:

                continue

            mid = sum(list(read_dict[read_type][read][start_pos][1])) / len(read_dict[read_type][read][start_pos][1])

            found = False

            for i in range(len(blocks)):

                start = blocks[i][0]
                end = blocks[i][1]

                if start <= mid and mid <= end:

                    if found: # FIXME: test

                        print("If this print it is because something is wrong with split blocks... A read is found in more than one block...")

                    blocks[i][2].update(read_dict[read_type][read][start_pos][1].copy())
                    blocks[i][3].append((read_type, read, start_pos))

                    found = True

            if found:

                read_dict[read_type][read][start_pos][3] = block

        new_blocks = []

        for i in range(len(blocks)):

            if len(blocks[i][3]) >= minimum_depth:

                new_blocks.append(blocks[i])
            
            elif len(blocks[i][3]) != 0:

                for read_type, read, start_pos in blocks[i][3]:

                    read_dict[read_type][read][start_pos][3] = "ungrouped"

        blocks = new_blocks

        if verbosity > 4:

            print(f'   {block}', [(block[0], block[1], len(block[2]), len(block[3])) for block in blocks], len(haploblock_dict[block][1]), len(haploblock_dict[block][0]), file=sys.stderr)

        for read_type, read, start_pos in haploblock_dict[block][1]:

            if read_dict[read_type][read][start_pos][3] == "ungrouped":

                reintroduce_read_to_variant_dict(read_type, read, start_pos)

                not_unique_list.append((read_type, read, start_pos))

        if len(blocks) == 0:

            del haploblock_dict[block]

        elif len(blocks) == 1:

            haploblock_dict[block] = [blocks[0][2], blocks[0][3], blocks[0][4], blocks[0][5]]

        else:

            if verbosity > 3:

                print(f'   {block} is split into {len(blocks)} subsections.', file=sys.stderr)

            for i in range(len(blocks)):

                haploblock_dict[f'{block}_{i+1}'] = [blocks[i][2], blocks[i][3], blocks[0][4], blocks[0][5]]

                for read_type, read, start_pos in blocks[i][3]:

                    read_dict[read_type][read][start_pos][3] = f'{block}_{i+1}'

            del haploblock_dict[block]

def compare_haploblocks(block1, block2, haploblock_dict):

    overlap = 0
    identity = 0
    mismatch_list = []

    for pos in haploblock_dict[block1][0]:

        if pos in haploblock_dict[block2][0]:

            overlap += 1

            if haploblock_dict[block1][0][pos] == haploblock_dict[block2][0][pos]:

                identity += 1

            else:

                mismatch_list.append((pos, haploblock_dict[block1][0][pos], haploblock_dict[block2][0][pos]))

    return overlap, identity, mismatch_list

def merge_haploblocks(block1, block2, mismatch_list, haploblock_dict): # Calls remove_read_from_variant_dict

    chrom = haploblock_dict[block1][2]

    for pos, _, _ in mismatch_list:

        if sum([pos in read_dict[read_type][read][start_pos][1] for read_type, read, start_pos in haploblock_dict[block2][1]]) > sum([pos in read_dict[read_type][read][start_pos][1] for read_type, read, start_pos in haploblock_dict[block1][1]]):

            hap = block1

        else:

            hap = block2

        read_count = 0

        for read_type, read, start_pos in haploblock_dict[hap][1]:

            if pos in read_dict[read_type][read][start_pos][1]:

                read_count += 1

                remove_read_from_variant_dict(pos, read_type, read, start_pos, haploblock_dict)

        if pos in haploblock_dict[hap][0]:

            del haploblock_dict[hap][0][pos]

        if verbosity > 3 and read_count > 0:

            if pos in variant_dict[chrom]:

                print(f"\tpos {pos} removed from the {read_count} {hap} read(s). Base dict is now {[(base, len(variant_dict[chrom][pos][base])) for base in variant_dict[chrom][pos]]}", file=sys.stderr)

            else:

                print(f"\tpos {pos} removed from the {read_count} {hap} read(s). This caused the base count to fall below the minimum count and the pos was removed.", file=sys.stderr)

    haploblock_dict[block1][0].update(haploblock_dict[block2][0].copy())
    haploblock_dict[block1][1] += haploblock_dict[block2][1]

    for read_type, read, start_pos in haploblock_dict[block2][1]:

        read_dict[read_type][read][start_pos][3] = block1

    del haploblock_dict[block2]

def low_variance_merge_haplotypes(max_diff, minimum_overlap, haploblock_dict): # Calls merge_haploblocks

    for block1 in list(haploblock_dict):

        if block1 not in haploblock_dict:

            continue

        if len(haploblock_dict[block1][1]) < minimum_coverage and max_diff != 0: # Makes sure that all blocks can be merged if there is no SNP differences between them.

            continue

        for block2 in list(haploblock_dict):

            if block1 == block2:

                continue

            if haploblock_dict[block1][2] != haploblock_dict[block2][2]:

                continue

            if len(haploblock_dict[block1][0]) == 0 or len(haploblock_dict[block2][0]) == 0:

                continue
        
            if not (min(haploblock_dict[block1][0]) <= max(haploblock_dict[block2][0]) and min(haploblock_dict[block2][0]) <= max(haploblock_dict[block1][0])):

                continue

            overlap, identity, mismatch_list = compare_haploblocks(block1, block2, haploblock_dict)

            if (overlap == identity and overlap >= minimum_overlap+len(mismatch_list)) or (len(mismatch_list) <= max_diff and overlap >= minimum_overlap + len(mismatch_list)):

                if sum([sum([base1 not in ["D", "I"] and base2 not in ["D", "I"] and pos in read_dict[read_type][read][start_pos][1] for read_type, read, start_pos in haploblock_dict[block2][1]]) > minimum_count and sum([pos in read_dict[read_type][read][start_pos][1] for read_type, read, start_pos in haploblock_dict[block1][1]]) > minimum_count for pos, base1, base2 in mismatch_list]) > 0:

                    # There is a high confidence pos in both blocks, so it is likely two different haplotypes.

                    continue

                if verbosity > 3:

                    print(f"   {block1} and {block2} is merged in low variance. Overlap is {overlap} and mismatch count is {len(mismatch_list)}. {mismatch_list}", file=sys.stderr)

                merge_haploblocks(block1, block2, mismatch_list, haploblock_dict)

def split_and_merge_primary_and_secondary_reads(max_diff, minimum_overlap, minimum_depth, haploblock_dict, not_unique_list): # Calls split_blocks and low_variance_merge_haplotypes

    # Splitting primary and secondary reads.

    temp_haploblock_dict = {}

    for block in list(haploblock_dict.keys()):

        for read_type, read, start_pos in haploblock_dict[block][1]:

            if read_type not in temp_haploblock_dict:

                temp_haploblock_dict[read_type] = {}

            new_block_name = f'{block}_{read_type}'

            if new_block_name not in temp_haploblock_dict[read_type]:

                temp_haploblock_dict[read_type][new_block_name] = [{}, [], haploblock_dict[block][2], haploblock_dict[block][3]]

            temp_haploblock_dict[read_type][new_block_name][0].update(read_dict[read_type][read][start_pos][1].copy())
            temp_haploblock_dict[read_type][new_block_name][1].append((read_type, read, start_pos))

            read_dict[read_type][read][start_pos][3] = new_block_name

    new_hapblock_dict = {}
    count = 0

    for group in temp_haploblock_dict:

        # Because the removed variants are only removed from the blocks in one read type group, they have to be removed after.

        for block_name in temp_haploblock_dict[group]:

            for pos in list(temp_haploblock_dict[group][block_name][0].keys()):

                if pos not in variant_dict[temp_haploblock_dict[group][block_name][2]]:

                    del temp_haploblock_dict[group][block_name][0][pos]

        split_blocks(minimum_depth, temp_haploblock_dict[group], not_unique_list)

        while True:

            haploblock_dict_pre_filter = len(temp_haploblock_dict[group])
            ungrouped_list_pre_filter = len(not_unique_list)

            low_variance_merge_haplotypes(0, minimum_overlap, temp_haploblock_dict[group])
            low_variance_merge_haplotypes(max_diff, minimum_overlap, temp_haploblock_dict[group])

            if haploblock_dict_pre_filter == len(temp_haploblock_dict[group]) and ungrouped_list_pre_filter == len(not_unique_list):

                break

        for block in temp_haploblock_dict[group]:

            count += 1
            block_name = f'block_{count}'

            new_hapblock_dict[block_name] = temp_haploblock_dict[group][block]

            for read_type, read, start_pos in temp_haploblock_dict[group][block][1]:

                read_dict[read_type][read][start_pos][3] = block_name

    # Because the removed variants are only removed from the blocks in one read type group, they have to be removed after.

    for block_name in list(new_hapblock_dict.keys()):

        if len(new_hapblock_dict[block_name][1]) < minimum_depth:

            not_unique_list += new_hapblock_dict[block_name][1]

            for read_type, read, start_pos in new_hapblock_dict[block_name][1]:

                reintroduce_read_to_variant_dict(read_type, read, start_pos)

            del new_hapblock_dict[block_name]

            continue

        chrom = new_hapblock_dict[block_name][2]

        for pos in list(new_hapblock_dict[block_name][0].keys()):

            if pos not in variant_dict[chrom]:

                del new_hapblock_dict[block_name][0][pos]

    return new_hapblock_dict, not_unique_list

def merge_low_coverage_blocks(haploblock_dict): # Calls make_fit_list

    for block in list(haploblock_dict.keys()):

        if block not in haploblock_dict:

            continue

        if len(haploblock_dict[block][0]) >= minimum_overlap:

            continue

        # The following makes sure that reads that have less SNPs than the minimum overlap has the chance to merge, if the identity is 100%.

        while True:

            merged = False

            fit_list = make_fit_list(haploblock_dict[block][0], haploblock_dict[block][2], haploblock_dict[block][3], haploblock_dict, 0, len(haploblock_dict[block][1]))

            if len(fit_list) > 1:

                for i in range(len(fit_list)):

                    block1 = block
                    block2 = fit_list[i][0]

                    if block1 != block2:

                        if verbosity > 3:

                            print(f"Low overlap merge between {block1} and {block2}. Overlap was {fit_list[i][1]} and identity was {fit_list[i][2]/fit_list[i][1]}.", file=sys.stderr)

                        haploblock_dict[block1][0].update(haploblock_dict[block2][0].copy())
                        haploblock_dict[block1][1] += haploblock_dict[block2][1]

                        for read_type, read, start_pos in haploblock_dict[block2][1]:

                            read_dict[read_type][read][start_pos][3] = block1

                        del haploblock_dict[block2]

                        merged = True

                        break

            if not merged:

                break

            if len(haploblock_dict[block][1]) < minimum_coverage: # FIXME: What is that for?...

                continue

def perform_no_variance_merging(haploblock_dict, not_unique_list): # Calls merge_low_coverage_blocks and split_and_merge_primary_and_secondary_reads

    print_subprocess_title("Performing no variance merging within mapping group.")

    if verbosity > 3:

        print(f'Performing low coverage merging.\n', file=sys.stderr)

    merge_low_coverage_blocks(haploblock_dict)

    if verbosity > 3:

        print(f'\nSplitting into mapping group (primary or secondary).\n', file=sys.stderr)

    haploblock_dict, not_unique_list = split_and_merge_primary_and_secondary_reads(0, 1, 1, haploblock_dict, not_unique_list)

    if verbosity > 3:

        print(f'', file=sys.stderr)

    if verbosity > 1:

        print(f'{len(haploblock_dict)} blocks remain after no variance merging.\n', file=sys.stderr)

    if args.output_intermediate_states:

        make_test_output(f'{args.outdir}/test_outputs/{args.prefix}_secondary_filtered_2_of_5_post_no_variance_merging.bam')

    return haploblock_dict, not_unique_list


def add_read_to_haploblock(read_type, read, start_pos, block, mismatch_list, haploblock_dict):

    chrom = read_dict[read_type][read][start_pos][0].reference_name

    for pos, _, _ in mismatch_list:

        if pos in read_dict[read_type][read][start_pos][1]:

            remove_read_from_variant_dict(pos, read_type, read, start_pos, haploblock_dict)

            if verbosity > 3:

                if pos not in variant_dict[chrom]:

                    print(f'\tpos {pos} removed for read to fit to {block}. This caused the base count to fall below the minimum count and the pos was removed.', file=sys.stderr)

                elif verbosity > 4:

                    print(f'\tpos {pos} removed for read to fit to {block}. Base dict is now {[(base, len(variant_dict[chrom][pos][base])) for base in variant_dict[chrom][pos]]}.', file=sys.stderr)

    haploblock_dict[block][0].update(read_dict[read_type][read][start_pos][1].copy())
    haploblock_dict[block][1].append((read_type, read, start_pos))

    read_dict[read_type][read][start_pos][3] = block

def reevaluate_unassigned_reads(max_diff, minimum_overlap, haploblock_dict, not_unique_list): # Calls make_fit_list, compare_haploblocks, merge_haploblocks and add_read_to_haploblock

    new_not_unique_list = []

    for read_type, read, start_pos in not_unique_list:

        read_dict[read_type][read][start_pos][1] = read_dict[read_type][read][start_pos][2].copy()

        read_variant_dict = read_dict[read_type][read][start_pos][2]

        # Removes read from ungrouped list, if the SNP count is 0

        if len(read_variant_dict) == 0:

            read_dict[read_type][read][start_pos][3] = "unphasable"

            continue

        fit_list = make_fit_list(read_variant_dict, read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, max_diff, minimum_overlap)

        if verbosity > 4:

            print(fit_list, len(read_variant_dict), (read_type, read, start_pos), file=sys.stderr)

        if len(fit_list) > 1:

            for i in range(len(fit_list)):
                for j in range(i+1, len(fit_list)):

                    block1 = fit_list[i][0]
                    block2 = fit_list[j][0]

                    if block1 not in haploblock_dict or block2 not in haploblock_dict:

                        continue

                    overlap, identity, mismatch_list = compare_haploblocks(block1, block2, haploblock_dict)

                    if (len(mismatch_list) <= max_diff and 0 < overlap) or (len(fit_list) == 2 and len(mismatch_list) <= max_diff):

                        if read_dict[read_type][read][start_pos][0].has_tag("SA") and overlap == 0:

                            if verbosity > 2:

                                print(f"The overlap between {block1} and {block2} is {overlap}, the identity is {identity} and a read overlaps both, but the read is a supplementary read, so no blocks are merged.", file=sys.stderr)

                            continue

                        if sum([sum([base1 not in ["D", "I"] and base2 not in ["D", "I"] and pos in read_dict[read_type][read][start_pos][1] for read_type, read, start_pos in haploblock_dict[block2][1]]) > minimum_count and sum([pos in read_dict[read_type][read][start_pos][1] for read_type, read, start_pos in haploblock_dict[block1][1]]) > minimum_count for pos, base1, base2 in mismatch_list]) > 0:

                            continue

                        if verbosity > 2:

                            print(f"The overlap between {block1} and {block2} is {overlap}, the identity is {identity} and a read overlaps both, so the two blocks are merged.", file=sys.stderr)

                        merge_haploblocks(block1, block2, mismatch_list, haploblock_dict)

            read_dict[read_type][read][start_pos][1] = read_dict[read_type][read][start_pos][2].copy()

            read_variant_dict = read_dict[read_type][read][start_pos][2]

            fit_list = make_fit_list(read_variant_dict, read_dict[read_type][read][start_pos][0].reference_name, haploblock_dict, max_diff, minimum_overlap)

        if len(fit_list) == 1:

            if verbosity > 4:

                print(f"   Read added to {block}", file=sys.stderr)

            add_read_to_haploblock(read_type, read, start_pos, fit_list[0][0], fit_list[0][3], haploblock_dict)

        elif len(fit_list) == 2:

            block1 = fit_list[0][0]
            block2 = fit_list[1][0]

            overlap, identity, mismatch_list = compare_haploblocks(block1, block2, haploblock_dict)

            if overlap == 0:

                # FIXME: could possibly work with the subfunctions...

                if verbosity > 2:

                    print(f"The overlap between {block1} and {block2} is {overlap} and a read overlaps both, so the two blocks are merged.", file=sys.stderr)

                merge_haploblocks(block1, block2, mismatch_list, haploblock_dict)

                add_read_to_haploblock(read_type, read, start_pos, block1, fit_list[0][3]+fit_list[1][3], haploblock_dict)

        elif len(fit_list) > 1:

            if read_dict[read_type][read][start_pos][0].has_tag("SA"):

                read_dict[read_type][read][start_pos][3] = "useful supplementary"

                new_not_unique_list.append((read_type, read, start_pos))

                continue

            max_block = None
            max_overlap = 0

            for i in range(len(fit_list)):

                if len(fit_list[i][3]) == 0 and max_overlap < fit_list[i][1]:

                    max_block = fit_list[i]
                    max_overlap = fit_list[i][1]

            if max_block != None:

                block = max_block[0]

                add_read_to_haploblock(read_type, read, start_pos, block, max_block[3], haploblock_dict)

            else:

                new_not_unique_list.append((read_type, read, start_pos))

        else:

            new_not_unique_list.append((read_type, read, start_pos))

    return new_not_unique_list

def merge_low_overlap_blocks(maxdiff, haploblock_dict):

    for block in list(haploblock_dict.keys()):

        if block not in haploblock_dict:

            continue

        if len(haploblock_dict[block][1]) < minimum_coverage:

            continue

        diff = 0

        while diff <= maxdiff:

            block_reduction = False

            fit_list = make_fit_list(haploblock_dict[block][0], haploblock_dict[block][2], haploblock_dict[block][3], haploblock_dict, diff, 1)

            if len(fit_list) > 0: # Checks if one of the mismatched SNPs is a high confidence SNP, in which case the read would not be moved.

                new_fit_list = []

                for i in range(len(fit_list)):

                    block1 = block

                    block2 = fit_list[i][0]

                    mismatch_list = fit_list[i][3]

                    if block1 != block2 and fit_list[i][2]/fit_list[i][1] > 0.5 and len(mismatch_list) > 0:

                        if sum([sum([base1 not in ["D", "I"] and base2 not in ["D", "I"] and pos in read_dict[read_type][read][start_pos][1] for read_type, read, start_pos in haploblock_dict[block2][1]]) > minimum_count and sum([pos in read_dict[read_type][read][start_pos][1] for read_type, read, start_pos in haploblock_dict[block1][1]]) > minimum_count for pos, base1, base2 in mismatch_list]) > 0:

                            continue

                        else:

                            new_fit_list.append(fit_list[i])

                    elif block1 != block2 and fit_list[i][2]/fit_list[i][1] > 0.5:

                        new_fit_list.append(fit_list[i])

                fit_list = new_fit_list

            if len(fit_list) > 0:

                fit_list = sorted(fit_list, key = lambda x:x[1], reverse = True)

                if verbosity > 3:

                    print(f"   {block} and {fit_list[0][0]} is merged in low overlap merging. Overlap is {fit_list[0][1]} and mismatch count is {len(fit_list[0][3])}.", file=sys.stderr)

                merge_haploblocks(block, fit_list[0][0], fit_list[0][3], haploblock_dict)

                block_reduction = True

            if not block_reduction:

                diff += 1

def reevaluate_remaining_ungrouped_reads(haploblock_dict, not_unique_list, merge_list = []): # Calls make_fit_list, compare_haploblocks, remove_read_from_variant_dict, merge_haploblocks and add_read_to_haploblock

    if verbosity > 3:

        print(f'', file=sys.stderr)

    if verbosity > 2:

        print("Reevaluating ungrouped reads.", file=sys.stderr)

    if verbosity > 3:

        print(f'', file=sys.stderr)

    diff = 0

    while diff <= maxdiff and len(not_unique_list) > 0:

        new_not_unique_list = []

        for read_type, read, start_pos in not_unique_list:

            read_dict[read_type][read][start_pos][1] = read_dict[read_type][read][start_pos][2].copy()

            read_variant_dict = read_dict[read_type][read][start_pos][2]

            if len(read_variant_dict) == 0:

                read_dict[read_type][read][start_pos][3] = "unphasable"

                continue

            if read_dict[read_type][read][start_pos][0].has_tag("SA") and not use_SAs: # FIXME: Don't know how well it handles supplementary alignments

                continue

            fit_list = make_fit_list(read_variant_dict, read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, diff, len(read_variant_dict)*minimum_pct_overlap)

            if verbosity > 4:

                print("This is from the remaining ungrouped reads:", fit_list, file=sys.stderr)

            if diff == 0 and len(fit_list) > 1:

                while True:

                    merged = False

                    for i in range(0, len(fit_list)):
                        for j in range(i+1, len(fit_list)):

                            block1 = fit_list[i][0]
                            block2 = fit_list[j][0]

                            overlap, identity, mismatch_list = compare_haploblocks(block1, block2, haploblock_dict)

                            if overlap == identity and overlap > 0:

                                haploblock_dict[block1][0].update(haploblock_dict[block2][0])
                                haploblock_dict[block1][1] += haploblock_dict[block2][1]

                                for read_type_2, read_2, start_pos_2 in haploblock_dict[block2][1]:

                                    read_dict[read_type_2][read_2][start_pos_2][3] = block1

                                del haploblock_dict[block2]

                                if verbosity > 3:

                                    print(f'   {block1} and {block2} is merged in reevaluation of ungrouped reads. Overlap is {overlap} and mismatch count is {len(mismatch_list)}.', file=sys.stderr)

                                if len(merge_list) > 0:

                                    merge_list = update_merge_list(block1, block2, haploblock_dict, merge_list)

                                merged = True

                            if merged:

                                break

                        if merged:

                            break

                    fit_list = make_fit_list(read_variant_dict, read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, diff, len(read_variant_dict)*minimum_pct_overlap)

                    if not merged or len(fit_list) == 1:

                        break

            if len(fit_list) == 1 or (diff == 0 and len(fit_list) > 1):

                if len(fit_list) > 1:

                    block = sorted([(hap, sum([sum([pos in read_dict[read_type_2][read_2][start_pos_2][2] for read_type_2, read_2, start_pos_2 in haploblock_dict[hap][1]]) for pos in read_dict[read_type][read][start_pos][2]])) for hap, _, _, _ in fit_list], key = lambda x:x[1])[-1][0]

                else:

                    block = fit_list[0][0]

                if verbosity > 3:

                    print(f'   Read added to {block}. Overlap is {fit_list[0][1]} and the mismatch count is {len(fit_list[0][3])}', file=sys.stderr)

                add_read_to_haploblock(read_type, read, start_pos, block, fit_list[0][3], haploblock_dict)

            elif len(fit_list) > 1:

                fit_list = sorted(fit_list, key = lambda x:x[1], reverse = True)

                if verbosity > 3:

                    print(f'   Read added to {fit_list[0][0]}. Overlap is {fit_list[0][1]} and the mismatch count is {len(fit_list[0][3])}', file=sys.stderr)

                add_read_to_haploblock(read_type, read, start_pos, fit_list[0][0], fit_list[0][3], haploblock_dict)

            else:

                new_not_unique_list.append((read_type, read, start_pos))

        if len(not_unique_list) == len(new_not_unique_list):

            diff += 1

        not_unique_list = new_not_unique_list

    # new_not_unique_list = []

    # for read_type, read, start_pos in not_unique_list:

    #     read_dict[read_type][read][start_pos][1] = read_dict[read_type][read][start_pos][2].copy()

    #     if len(read_dict[read_type][read][start_pos][2]) == 0:

    #         read_dict[read_type][read][start_pos][3] == "unphasable"

    #         continue

    #     fit_list = make_fit_list(read_dict[read_type][read][start_pos][2], read_dict[read_type][read][start_pos][0].reference_name, haploblock_dict, len(read_dict[read_type][read][start_pos][2]), 1)

    #     print(len(read_dict[read_type][read][start_pos][2]), [len(haploblock_dict[block][0]) for block in haploblock_dict])

    #     if verbosity > 4:

    #         print("This is from the remaining ungrouped reads, the last merging:", fit_list, file=sys.stderr)

    #     fit_list = sorted(fit_list, key = lambda x:x[2], reverse = True)

    #     if len(fit_list) > 0:

    #         if verbosity > 3:

    #             print(f'   Read added to {block}. Overlap is {fit_list[0][1]} and the mismatch count is {len(fit_list[0][3])}', file=sys.stderr)

    #         add_read_to_haploblock(read_type, read, start_pos, fit_list[0][0], fit_list[0][3], haploblock_dict)
        
    #     else:
                
    #         new_not_unique_list.append((read_type, read, start_pos))

    # not_unique_list = new_not_unique_list

    return not_unique_list, merge_list

def remove_read_pos_from_variant_dict(mismatch_list, read_type, read, start_pos, haploblock_dict): # Calls reintroduce_read_to_variant_dict, remove_pos_from_haploblocks and remove_variant

    chrom = read_dict[read_type][read][start_pos][0].reference_name

    pre_filter = read_dict[read_type][read][start_pos][3]

    reintroduce_read_to_variant_dict(read_type, read, start_pos)

    read_dict[read_type][read][start_pos][3] = pre_filter

    for pos, _, _ in mismatch_list:

        base = read_dict[read_type][read][start_pos][2][pos]

        for i in range(len(variant_dict[chrom][pos][base])):

            if variant_dict[chrom][pos][base][i] == (read_type, read, start_pos):

                del variant_dict[chrom][pos][base][i]

                break

        # If the base count is less than the minimum count the base is removed from variant dict and reads in read dict.

        if len(variant_dict[chrom][pos][base]) < minimum_count:

            for read_type_1, read_1, start_pos_1 in variant_dict[chrom][pos][base]:

                del read_dict[read_type_1][read_1][start_pos_1][1][pos]
                del read_dict[read_type_1][read_1][start_pos_1][2][pos]

                remove_pos_from_haploblocks(read_type_1, read_1, start_pos_1, pos, haploblock_dict)

            del variant_dict[chrom][pos][base]

            if verbosity > 3:

                print(f"The base {base} has been removed from pos {pos}.", file=sys.stderr)

            if "deleted" in variant_dict[chrom][pos]:

                new_deleted_list = []

                for read_type_1, read_1, start_pos_1 in variant_dict[chrom][pos]["deleted"]:

                    if read_dict[read_type_1][read_1][start_pos_1][2][pos] == base:

                        del read_dict[read_type_1][read_1][start_pos_1][2][pos]

                        remove_pos_from_haploblocks(read_type_1, read_1, start_pos_1, pos, haploblock_dict)

                    else:

                        new_deleted_list.append((read_type_1, read_1, start_pos_1))

                if len(new_deleted_list) == 0:

                    del variant_dict[chrom][pos]["deleted"]

                else:

                    variant_dict[chrom][pos]["deleted"] = new_deleted_list

        if "deleted" in variant_dict[chrom][pos]:

            for i in range(len(variant_dict[chrom][pos]["deleted"])):

                if variant_dict[chrom][pos]["deleted"][i] == (read_type, read, start_pos):

                    del variant_dict[chrom][pos]["deleted"][i]

                    break

            if len(variant_dict[chrom][pos]["deleted"]) == 0:

                del variant_dict[chrom][pos]["deleted"]

        del read_dict[read_type][read][start_pos][2][pos]

        if sum([base in ["deleted", "I", "D", fastafile.fetch(chrom, pos, pos+1).upper()] for base in variant_dict[chrom][pos]]) == len(variant_dict[chrom][pos]) or sum([base not in ["deleted"] for base in variant_dict[chrom][pos]]) == 1:

            if verbosity > 4:
                
                print(f'{pos} is removed.', file=sys.stderr)

            remove_variant(chrom, pos, haploblock_dict, read_dict)

        elif "deleted" in variant_dict[chrom][pos] and len(variant_dict[chrom][pos]["deleted"])/sum([len(variant_dict[chrom][pos][base]) for base in variant_dict[chrom][pos]]) > max_variance: # If the deleted base count exceeds max_variance it is likely a bad variant and should therefore be removed.

            if verbosity > 4:
                
                print(f'{pos} is removed.', file=sys.stderr)

            remove_variant(chrom, pos, haploblock_dict, read_dict)

        elif sum([len(variant_dict[chrom][pos][base]) if base != "deleted" else 0 for base in variant_dict[chrom][pos]])/depth_dict[chrom][pos] < 1-max_variance:

            if verbosity > 4:
                
                print(f'{pos} is removed.', file=sys.stderr)

            remove_variant(chrom, pos, haploblock_dict, read_dict)

        if verbosity > 3:

            if pos not in variant_dict[chrom]:

                print(f"\tpos {pos} removed from read. This caused the base count to fall below the minimum count and the pos was removed.", file=sys.stderr)

    read_dict[read_type][read][start_pos][1] = read_dict[read_type][read][start_pos][2].copy()

def cleanup_haploblock_variant_dicts(haploblock_dict): # FIXME: NOT USED

    # Clean up haploblock variant dicts

    for block in haploblock_dict:

        haploblock_variant_dict = haploblock_dict[block][0]

        chrom = haploblock_dict[block][2]

        for pos in list(haploblock_variant_dict.keys()):

            if pos not in variant_dict[chrom]: # FIXME: why necessary?

                del haploblock_dict[block][0][pos]

                continue

            count_dict = {}

            for base in variant_dict[chrom][pos]:

                if base != "deleted":

                    count_dict[base] = 0

                for read_type, read, start_pos in variant_dict[chrom][pos][base]:

                    if read_dict[read_type][read][start_pos][3] == block:

                        if base == "deleted" and pos in read_dict[read_type][read][start_pos][2]:

                            count_dict[read_dict[read_type][read][start_pos][2][pos]] += 1

                        else:

                            count_dict[base] += 1

            max_base = max(count_dict, key = count_dict.get)

            if haploblock_variant_dict[pos] != max_base and len(haploblock_variant_dict[pos]) == 1:

                if verbosity > 3:

                    print("changing variant", block, pos, haploblock_variant_dict[pos], max_base, count_dict, file=sys.stderr)

                haploblock_variant_dict[pos] = max_base

def cleanup_variants(haploblock_dict): # Calls remove_variant, remove_read_pos_from_variant_dict, make_fit_list and add_read_to_haploblock

    if verbosity > 3:

        print(f'', file=sys.stderr)

    if verbosity > 2:

        print("Removing unusable read variants.", file=sys.stderr)

    if verbosity > 3:

        print(f'', file=sys.stderr)

    #cleanup_haploblock_variant_dicts(haploblock_dict)

    for chrom in variant_dict:
            
        while True:

            pre_filter_variant = len(variant_dict[chrom])

            for pos in list(variant_dict[chrom].keys()):

                new_dict = {}

                for base in variant_dict[chrom][pos]:

                    new_dict[base] = {}

                    for read_type, read, start_pos in variant_dict[chrom][pos][base]:

                        if read_dict[read_type][read][start_pos][3] not in new_dict[base]:

                            new_dict[base][read_dict[read_type][read][start_pos][3]] = 0

                        new_dict[base][read_dict[read_type][read][start_pos][3]] += 1

                if verbosity > 4:
                    
                    print(pos, [(base, len(variant_dict[chrom][pos][base])) for base in variant_dict[chrom][pos]], new_dict, file=sys.stderr)

                if sum([(len(new_dict[base]) == 1 and "ungrouped" in new_dict[base]) or base == "deleted" for base in new_dict]) >= len(new_dict)-1:

                    if verbosity > 4:
                        
                        print(f'{pos} is removed.', file=sys.stderr)

                    remove_variant(chrom, pos, haploblock_dict, read_dict)


                else: # FIXME: understand usability

                    for base in new_dict:

                        if len(new_dict[base]) == 1 and "ungrouped" in new_dict[base]:

                            for read_type, read, start_pos in variant_dict[chrom][pos][base]:

                                remove_read_pos_from_variant_dict([(pos, None, None)], read_type, read, start_pos, haploblock_dict)

                                if pos not in variant_dict[chrom] or base not in variant_dict[chrom][pos]:

                                    break

                            for block in haploblock_dict:

                                if pos in haploblock_dict[block][0]:

                                    del haploblock_dict[block][0][pos]

                        if pos not in variant_dict[chrom]:

                            break

            new_ungrouped_list = []

            for read_type in read_dict:
                for read in read_dict[read_type]:
                    for start_pos in read_dict[read_type][read]:

                        if read_dict[read_type][read][start_pos][3] in ["unchosen", "removed"]:

                            continue

                        if read_dict[read_type][read][start_pos][0].has_tag("SA"):

                            continue

                        if read_dict[read_type][read][start_pos][3] not in haploblock_dict and read_dict[read_type][read][start_pos][3] != "unphasable": # FIXME: somewhere reads are lost. This fixes that mess-up. But the real problem should be found...

                            read_dict[read_type][read][start_pos][3] = "ungrouped"
                            new_ungrouped_list.append((read_type, read, start_pos))

                        if len(read_dict[read_type][read][start_pos][2]) == 0:

                            continue

                        diff = 0

                        while diff <= maxdiff:

                            fit_list = make_fit_list(read_dict[read_type][read][start_pos][2], read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, diff, len(read_dict[read_type][read][start_pos][2])*minimum_pct_overlap)

                            if len(fit_list) == 1 and len(fit_list[0][3]) > 0:

                                if verbosity > 4:

                                    print("Variant cleanup:", read_type, read, start_pos, read_dict[read_type][read][start_pos][3], fit_list, file=sys.stderr)

                                pre_group = read_dict[read_type][read][start_pos][3]

                                remove_read_pos_from_variant_dict(fit_list[0][3], read_type, read, start_pos, haploblock_dict)

                                read_dict[read_type][read][start_pos][3] = pre_group

                                break

                            elif len(fit_list) > 1 and sum([len(fit_list[i][3]) > 0 for i in range(len(fit_list))]) == len(fit_list):

                                if verbosity > 4:
                                    
                                    print("Variant cleanup:", read_type, read, start_pos, read_dict[read_type][read][start_pos][3], fit_list, file=sys.stderr)

                                mismatch_list = []

                                for pos, _, _ in fit_list[0][3]:

                                    for i in range(1, len(fit_list)):

                                        found = False

                                        for pos2, _, _ in fit_list[i][3]:

                                            if pos == pos2:

                                                found = True

                                        if not found:

                                            break

                                    if found:

                                        mismatch_list.append((pos, None, None))

                                if len(mismatch_list) > 0:

                                    if verbosity > 4:

                                        print(f"The common variants are: {mismatch_list}", file=sys.stderr)

                                    pre_group = read_dict[read_type][read][start_pos][3]

                                    remove_read_pos_from_variant_dict(mismatch_list, read_type, read, start_pos, haploblock_dict)

                                    read_dict[read_type][read][start_pos][3] = pre_group
                                
                                break
                            
                            else:

                                diff += 1

                        if read_dict[read_type][read][start_pos][3] == "ungrouped":

                            new_ungrouped_list.append((read_type, read, start_pos))

            if pre_filter_variant == len(variant_dict[chrom]):

                break

    return new_ungrouped_list

def remove_mismatched_variants_from_read(mismatch_list, read_type, read, start_pos, haploblock_dict): # Calls remove_read_from_variant_dict

    chrom = read_dict[read_type][read][start_pos][0].reference_name

    for pos, _, _ in mismatch_list:

        if pos in read_dict[read_type][read][start_pos][2]:

            remove_read_from_variant_dict(pos, read_type, read, start_pos, haploblock_dict)

            if verbosity > 3:

                if pos not in variant_dict[chrom]:

                    print(f"\tpos {pos} removed from read. This caused the base count to fall below the minimum count and the pos was removed.", file=sys.stderr)

                elif verbosity > 4:

                    print(f"\tpos {pos} removed from read. Base dict is now {[(base, len(variant_dict[chrom][pos][base])) for base in variant_dict[chrom][pos]]}", file=sys.stderr)

def cleanup_blocks(minimum_depth, maxdiff, haploblock_dict, not_unique_list): # Calls make_fit_list, reintroduce_read_to_variant_dict and remove_mismatched_variants_from_read

    blocks_by_coverage = sorted([(block, len(haploblock_dict[block][1])) for block in haploblock_dict], key = lambda x:x[1])

    for block2, _ in blocks_by_coverage:

        # if verbosity > 4:

        #     print(f'The old coverage for {block2} is {len(haploblock_dict[block2][1])}.', file=sys.stderr)

        chrom = haploblock_dict[block2][2]

        new_read_variant_dict = {}
        new_read_list = []

        for read_type, read, start_pos in haploblock_dict[block2][1]:

            read_variant_dict = read_dict[read_type][read][start_pos][2]

            if len(read_variant_dict) == 0:

                read_dict[read_type][read][start_pos][3] = "unphasable"

                continue

            diff = 0
            found = False

            while diff <= maxdiff:

                fit_list = make_fit_list(read_variant_dict, read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, diff, 1)

                if len(fit_list) == 0:

                    diff += 1

                else:

                    fit_list = sorted([(fit_list[i], sum([sum([pos in read_dict[read_type_2][read_2][start_pos_2][2] and read_dict[read_type_2][read_2][start_pos_2][2][pos] == read_dict[read_type][read][start_pos][2][pos] for read_type_2, read_2, start_pos_2 in haploblock_dict[fit_list[i][0]][1]]) for pos in read_dict[read_type][read][start_pos][2]])) for i in range(len(fit_list))], key = lambda x:x[1], reverse = True)

                    fit_list = [fit_list[i][0] for i in range(len(fit_list))]

                    block = fit_list[0][0]

                    reintroduce_read_to_variant_dict(read_type, read, start_pos)

                    read_dict[read_type][read][start_pos][3] = block

                    if diff > 0:

                        remove_mismatched_variants_from_read(fit_list[0][3], read_type, read, start_pos, haploblock_dict)

                    if block == block2:

                        new_read_variant_dict.update(read_dict[read_type][read][start_pos][1].copy())
                        new_read_list.append((read_type, read, start_pos))

                    else:

                        haploblock_dict[block][0].update(read_dict[read_type][read][start_pos][1].copy())
                        haploblock_dict[block][1].append((read_type, read, start_pos))

                    found = True

                    break

            if not found:

                reintroduce_read_to_variant_dict(read_type, read, start_pos)

                not_unique_list.append((read_type, read, start_pos))

        # if verbosity > 4:

        #     print(f'The new coverage for {block2} is {len(new_read_list)}.', file=sys.stderr)

        if len(new_read_list) < minimum_depth:

            del haploblock_dict[block2]

            if verbosity > 3:

                print(f'   {block2} is removed due to read count falling below {minimum_depth}.', file=sys.stderr)

            if len(new_read_list) > 0:

                for read_type, read, start_pos in new_read_list:

                    reintroduce_read_to_variant_dict(read_type, read, start_pos)

                not_unique_list += new_read_list

        else:

            for pos in list(new_read_variant_dict.keys()):

                if pos not in variant_dict[chrom]:

                    del new_read_variant_dict[pos]

            haploblock_dict[block2] = [new_read_variant_dict, new_read_list, chrom, haploblock_dict[block2][3]]

    split_blocks(minimum_depth, haploblock_dict, not_unique_list)

    new_haploblock_dict = {}

    count = 1

    for block in haploblock_dict:

        # This part makes sure no removed variants are carried into the haploblocks

        for pos in list(haploblock_dict[block][0].keys()):

            if pos not in variant_dict[chrom]:

                del haploblock_dict[block][0][pos]

        new_block_name = f'block_{count}'

        new_haploblock_dict[new_block_name] = haploblock_dict[block]

        for read_type, read, start_pos in new_haploblock_dict[new_block_name][1]:

            read_dict[read_type][read][start_pos][3] = new_block_name

        count += 1

    return new_haploblock_dict, not_unique_list

def make_unphasable_consensus(haploblock_dict, not_unique_list, split_haplotypes = True): # make_fit_list, reintroduce_read_to_variant_dict, split_blocks, low_variance_merge_haplotypes

    if verbosity > 3:

        print(f'', file=sys.stderr)

    if verbosity > 2:

        print("Cleanup of ambigously mapped reads.", file=sys.stderr)

    if verbosity > 3:

        print(f'', file=sys.stderr)

    for block in haploblock_dict:

        new_read_variant_dict = {}
        new_read_list = []

        for read_type, read, start_pos in haploblock_dict[block][1]:

            read_variant_dict = read_dict[read_type][read][start_pos][2]

            if len(read_variant_dict) == 0:

                if read_dict[read_type][read][start_pos][3] != "unchosen":

                    read_dict[read_type][read][start_pos][3] = "unphasable"

                continue

            if read_dict[read_type][read][start_pos][3] == "ungrouped":

                continue

            fit_list = make_fit_list(read_variant_dict, read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, 0, len(read_variant_dict)*minimum_pct_overlap)

            if len(fit_list) == 1:

                new_read_variant_dict.update(read_variant_dict.copy())
                new_read_list.append((read_type, read, start_pos))

                read_dict[read_type][read][start_pos][3] = block

            else:

                if verbosity > 4:

                    print("Unphasable consensus:", read_type, read, start_pos, read_dict[read_type][read][start_pos][3], fit_list, file=sys.stderr)

                not_unique_list.append((read_type, read, start_pos))
                reintroduce_read_to_variant_dict(read_type, read, start_pos)

        haploblock_dict[block] = [new_read_variant_dict, new_read_list, haploblock_dict[block][2], haploblock_dict[block][3]]

    if split_haplotypes:

        split_blocks(1, haploblock_dict, not_unique_list)

    return not_unique_list

def split_and_merge_unique_and_multimapped_reads(max_diff, minimum_overlap, minimum_depth, haploblock_dict, not_unique_list): # Calls split_blocks and low_variance_merge_haplotypes

    # Splitting primary and secondary reads.

    temp_haploblock_dict = {"u": {}, "m": {}}

    for block in list(haploblock_dict.keys()):

        for read_type, read, start_pos in haploblock_dict[block][1]:

            if read in read_dict["p"] and read in read_dict["s"]:

                read_group = "m"
            
            else:

                read_group = "u"

            if read_group not in temp_haploblock_dict:

                temp_haploblock_dict[read_group] = {}

            new_block_name = f'{block}_{read_group}'

            if new_block_name not in temp_haploblock_dict[read_group]:

                temp_haploblock_dict[read_group][new_block_name] = [{}, []]

            temp_haploblock_dict[read_group][new_block_name][0].update(read_dict[read_type][read][start_pos][1].copy())
            temp_haploblock_dict[read_group][new_block_name][1].append((read_type, read, start_pos))

            read_dict[read_type][read][start_pos][3] = new_block_name

    new_hapblock_dict = {}
    count = 0

    for group in temp_haploblock_dict:

        # Because the removed variants are only removed from the blocks in one read type group, they have to be removed after.

        for block_name in temp_haploblock_dict[group]:

            for pos in list(temp_haploblock_dict[group][block_name][0].keys()):

                if pos not in variant_dict:

                    del temp_haploblock_dict[group][block_name][0][pos]

        split_blocks(minimum_depth, temp_haploblock_dict[group], not_unique_list)

        while True:

            haploblock_dict_pre_filter = len(temp_haploblock_dict[group])
            ungrouped_list_pre_filter = len(not_unique_list)

            low_variance_merge_haplotypes(0, minimum_overlap, temp_haploblock_dict[group])
            low_variance_merge_haplotypes(max_diff, minimum_overlap, temp_haploblock_dict[group])

            if haploblock_dict_pre_filter == len(temp_haploblock_dict[group]) and ungrouped_list_pre_filter == len(not_unique_list):

                break

        for block in temp_haploblock_dict[group]:

            count += 1
            block_name = f'block_{count}'

            new_hapblock_dict[block_name] = temp_haploblock_dict[group][block]

            for read_type, read, start_pos in temp_haploblock_dict[group][block][1]:

                read_dict[read_type][read][start_pos][3] = block_name

    # Because the removed variants are only removed from the blocks in one read type group, they have to be removed after.

    for block_name in list(new_hapblock_dict.keys()):

        if len(new_hapblock_dict[block_name][1]) < minimum_depth:

            not_unique_list += new_hapblock_dict[block_name][1]

            for read_type, read, start_pos in new_hapblock_dict[block_name][1]:

                reintroduce_read_to_variant_dict(read_type, read, start_pos)

            del new_hapblock_dict[block_name]

            continue

        for pos in list(new_hapblock_dict[block_name][0].keys()):

            if pos not in variant_dict:

                del new_hapblock_dict[block_name][0][pos]

    return new_hapblock_dict, not_unique_list

def final_merging(haploblock_dict, not_unique_list):

    # This function makes sure we merge all blocks that do not have more than one match, where the matches do not match themselves. That would indicate an area, where we cannot be sure if a block fits one or the other, and det blocks should therefore be kept seperate. 

    merge_dict = {}

    for block in haploblock_dict:

        fit_list = make_fit_list(haploblock_dict[block][0], haploblock_dict[block][2], haploblock_dict[block][3], haploblock_dict, 0, 1)

        for i in range(len(fit_list)):

            if fit_list[i][0] != block:

                if block not in merge_dict:

                    merge_dict[block] = []

                merge_dict[block].append(fit_list[i][0])

    for block in list(merge_dict.keys()):

        if block not in merge_dict:

            continue

        if len(merge_dict[block]) != 0:

            merge_list = merge_dict[block].copy()

            i = 0

            while i < len(merge_list):

                for block2 in merge_dict[merge_list[i]]:

                    if block2 not in merge_list:

                        merge_list.append(block2)

                i += 1

            for i in range(len(merge_list)):
                for j in range(i+1, len(merge_list)):

                    overlap, identity, mismatch_list = compare_haploblocks(merge_list[i], merge_list[j], haploblock_dict)

                    if len(mismatch_list) > 0:

                        if merge_list[i] in merge_dict[block]:

                            merge_dict[block].remove(merge_list[i])

                        if merge_list[j] in merge_dict[block]:

                            merge_dict[block].remove(merge_list[j])

                        if block in merge_dict[merge_list[i]]:

                            merge_dict[merge_list[i]].remove(block)

                        if block in merge_dict[merge_list[j]]:

                            merge_dict[merge_list[j]].remove(block)

        if len(merge_dict[block]) == 0:

            del merge_dict[block]

    for block in merge_dict:

        if block not in haploblock_dict:

            continue

        merge_list = merge_dict[block].copy()

        i = 0

        while i < len(merge_list):

            for block2 in merge_dict[merge_list[i]]:

                if block2 not in merge_list and block2 != block:

                    merge_list.append(block2)

            i += 1

        for i in range(len(merge_list)):

            merge_haploblocks(block, merge_list[i], [], haploblock_dict)

    new_not_unique_list = []

    for read_type, read, start_pos in not_unique_list:

        fit_list = make_fit_list(read_dict[read_type][read][start_pos][2], read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, 0, 1)

        if verbosity > 4:

            print("Final merging:", read_type, read, start_pos, read_dict[read_type][read][start_pos][3], fit_list, file=sys.stderr)

        if len(fit_list) == 1:

            add_read_to_haploblock(read_type, read, start_pos, fit_list[0][0], fit_list[0][3], haploblock_dict)
        
        elif len(fit_list) == 2:

            identity, overlap, mismatch_list = compare_haploblocks(fit_list[0][0], fit_list[1][0], haploblock_dict)

            if verbosity > 4:

                print(identity, overlap, len(mismatch_list), file=sys.stderr)

            if identity == overlap:

                merge_haploblocks(fit_list[0][0], fit_list[1][0], mismatch_list, haploblock_dict)

                add_read_to_haploblock(read_type, read, start_pos, fit_list[0][0], fit_list[0][3], haploblock_dict)

                if verbosity > 3:

                    print(f'   {fit_list[0][0]} and {fit_list[1][0]} is merged in final merging due to an ungrouped reads overlapping both. Overlap is {overlap} and mismatch count is {len(mismatch_list)}.', file=sys.stderr)
            
            else:

                new_not_unique_list.append((read_type, read, start_pos))

        else:

            new_not_unique_list.append((read_type, read, start_pos))
    
    not_unique_list = new_not_unique_list

    new_haploblock_dict = {}

    count = 1

    for block in haploblock_dict:

        new_block_name = f'block_{count}'

        new_haploblock_dict[new_block_name] = haploblock_dict[block]

        for read_type, read, start_pos in haploblock_dict[block][1]:

            read_dict[read_type][read][start_pos][3] = new_block_name

        count += 1

    return new_haploblock_dict, not_unique_list

def perform_low_variance_merging(haploblock_dict, not_unique_list): # Calls low_variance_merge_haplotypes, reevaluate_unassigned_reads, split_and_merge_primary_and_secondary_reads, merge_low_overlap_blocks, reevaluate_remaining_ungrouped_reads and perform_ambigous_block_cleanup

    print_subprocess_title("Performing low variance merging.")

    count = 0
    diff = 0

    if verbosity > 2:

        print(f"Starting merging with the maximum difference: {diff}.\n", file=sys.stderr)

    while True:

        count += 1

        variant_dict_pre_filter = sum([len(variant_dict[chrom]) for chrom in variant_dict])
        haploblock_dict_pre_filter = len(haploblock_dict)

        low_variance_merge_haplotypes(0, minimum_overlap, haploblock_dict)
        low_variance_merge_haplotypes(diff, minimum_overlap, haploblock_dict)

        if haploblock_dict_pre_filter == len(haploblock_dict) and variant_dict_pre_filter == sum([len(variant_dict[chrom]) for chrom in variant_dict]):

            if verbosity > 3:

                print(f'', file=sys.stderr)

            if verbosity > 2:

                print(f'After {count} iteration(s) of low variance merging the block count is no longer decreasing, so low overlap merging is performed.', file=sys.stderr)

            if verbosity > 3:

                print(f'', file=sys.stderr)

            merge_low_overlap_blocks(diff, haploblock_dict)

            if haploblock_dict_pre_filter != len(haploblock_dict) or variant_dict_pre_filter != sum([len(variant_dict[chrom]) for chrom in variant_dict]):

                if verbosity > 3:

                    print(f'', file=sys.stderr)

                if verbosity > 2:

                    print(f'Continues low variance merging due to low overlap merging causing removal of {variant_dict_pre_filter-sum([len(variant_dict[chrom]) for chrom in variant_dict])} additional variant(s) and/or {haploblock_dict_pre_filter-len(haploblock_dict)} additional blocks.', file=sys.stderr)

                if verbosity > 3:

                    print(f'', file=sys.stderr)

                continue

            if verbosity > 3:

                print(f'Performing ambigous mapping cleanup.\n', file=sys.stderr)

            haploblock_dict, not_unique_list = cleanup_blocks(1, diff, haploblock_dict, not_unique_list)

            low_variance_merge_haplotypes(0, minimum_overlap, haploblock_dict)
            low_variance_merge_haplotypes(diff, minimum_overlap, haploblock_dict)

            if variant_dict_pre_filter != sum([len(variant_dict[chrom]) for chrom in variant_dict]):

                if verbosity > 3:

                    print(f'', file=sys.stderr)

                if verbosity > 2:

                    print(f'Continues low variance merging due to ambigous mapping cleanup causing removal of {variant_dict_pre_filter-sum([len(variant_dict[chrom]) for chrom in variant_dict])} additional variant(s).', file=sys.stderr)

                if verbosity > 3:

                    print(f'', file=sys.stderr)

                continue

            if maxdiff > diff:

                if verbosity > 3:

                    print(f'', file=sys.stderr)

                if verbosity > 2:

                    print(f"Low overlap merging and ambigous mapping cleanup did not remove additional variant(s), so the allowed difference is increased.", file=sys.stderr)

                diff += 1

                if verbosity > 2:

                    print(f"\nIncreasing the maximum difference to: {diff}.\n", file=sys.stderr)

            else:

                if verbosity > 3:

                    print(f'', file=sys.stderr)

                if verbosity > 2:

                    print(f"Low overlap merging and ambigous mapping cleanup did not remove additional variant(s), so ambigous mapping block cleanup with required minimum depth of {minimum_depth} is initiated.", file=sys.stderr)

                if verbosity > 3:

                    print(f'', file=sys.stderr)

                variant_dict_pre_filter = sum([len(variant_dict[chrom]) for chrom in variant_dict])

                haploblock_dict, not_unique_list = cleanup_blocks(minimum_depth, maxdiff, haploblock_dict, not_unique_list)

                while True:

                    haploblock_dict_pre_filter_2 = len(haploblock_dict)
                    ungrouped_list_pre_filter_2 = len(not_unique_list)

                    low_variance_merge_haplotypes(diff, minimum_overlap, haploblock_dict)

                    merge_low_overlap_blocks(maxdiff, haploblock_dict)

                    not_unique_list, _ = reevaluate_remaining_ungrouped_reads(haploblock_dict, not_unique_list)

                    if haploblock_dict_pre_filter_2 == len(haploblock_dict) and ungrouped_list_pre_filter_2 == len(not_unique_list) or variant_dict_pre_filter != sum([len(variant_dict[chrom]) for chrom in variant_dict]):

                        break

                if variant_dict_pre_filter != sum([len(variant_dict[chrom]) for chrom in variant_dict]):

                    if verbosity > 2:

                        print(f'Continues low variance merging due to ambigous mapping cleanup causing removal of {variant_dict_pre_filter-sum([len(variant_dict[chrom]) for chrom in variant_dict])} additional variant(s).', file=sys.stderr)

                    continue

                haploblock_dict, not_unique_list = split_and_merge_primary_and_secondary_reads(maxdiff, 1, minimum_depth, haploblock_dict, not_unique_list)

                haploblock_dict, not_unique_list = cleanup_blocks(1, maxdiff, haploblock_dict, not_unique_list)
                
                while True:

                    haploblock_dict_pre_filter_2 = len(haploblock_dict)
                    ungrouped_list_pre_filter_2 = len(not_unique_list)

                    not_unique_list, _ = reevaluate_remaining_ungrouped_reads(haploblock_dict, not_unique_list)

                    merge_low_overlap_blocks(maxdiff, haploblock_dict)

                    low_variance_merge_haplotypes(diff, minimum_overlap, haploblock_dict)

                    if haploblock_dict_pre_filter_2 == len(haploblock_dict) and ungrouped_list_pre_filter_2 == len(not_unique_list) or variant_dict_pre_filter != sum([len(variant_dict[chrom]) for chrom in variant_dict]):

                        break

                if variant_dict_pre_filter != sum([len(variant_dict[chrom]) for chrom in variant_dict]):

                    if verbosity > 2:

                        print(f'Continues low variance merging due to ambigous mapping cleanup causing removal of {variant_dict_pre_filter-sum([len(variant_dict[chrom]) for chrom in variant_dict])} additional variant(s).', file=sys.stderr)

                    continue

                not_unique_list = cleanup_variants(haploblock_dict)

                not_unique_list = make_unphasable_consensus(haploblock_dict, not_unique_list, False)

                i = 1

                while i < max_cleanup_iterations:

                    haploblock_dict_pre_filter = len(haploblock_dict)
                    ungrouped_list_pre_filter = len(not_unique_list)

                    split_blocks(minimum_depth, haploblock_dict, not_unique_list)

                    while True:

                        haploblock_dict_pre_filter_2 = len(haploblock_dict)
                        ungrouped_list_pre_filter_2 = len(not_unique_list)

                        merge_low_overlap_blocks(0, haploblock_dict)

                        not_unique_list, _ = reevaluate_remaining_ungrouped_reads(haploblock_dict, not_unique_list)

                        if haploblock_dict_pre_filter_2 == len(haploblock_dict) and ungrouped_list_pre_filter_2 == len(not_unique_list):

                            break

                    if variant_dict_pre_filter != sum([len(variant_dict[chrom]) for chrom in variant_dict]):

                        break
                        
                    not_unique_list = cleanup_variants(haploblock_dict)

                    not_unique_list = make_unphasable_consensus(haploblock_dict, not_unique_list, False)

                    if haploblock_dict_pre_filter == len(haploblock_dict) and ungrouped_list_pre_filter == len(not_unique_list):

                        break

                    i += 1

                if variant_dict_pre_filter != sum([len(variant_dict[chrom]) for chrom in variant_dict]):

                    if verbosity > 2:

                        print(f'Continues low variance merging due to ambigous mapping cleanup causing removal of {variant_dict_pre_filter-sum([len(variant_dict[chrom]) for chrom in variant_dict])} additional variant(s).', file=sys.stderr)

                    continue

                split_blocks(minimum_depth, haploblock_dict, not_unique_list)

                while True:

                    haploblock_dict_pre_filter = len(haploblock_dict)
                    ungrouped_list_pre_filter = len(not_unique_list)

                    haploblock_dict, not_unique_list = final_merging(haploblock_dict, not_unique_list)

                    if haploblock_dict_pre_filter == len(haploblock_dict) and ungrouped_list_pre_filter == len(not_unique_list):

                        break

                if variant_dict_pre_filter != sum([len(variant_dict[chrom]) for chrom in variant_dict]):

                    if verbosity > 2:

                        print(f'Continues low variance merging due to ambigous mapping cleanup causing removal of {variant_dict_pre_filter-sum([len(variant_dict[chrom]) for chrom in variant_dict])} additional variant(s).', file=sys.stderr)

                    continue

                if verbosity > 3:

                    print(f'', file=sys.stderr)

                if verbosity > 2:

                    print(f'After {count-1} iteration(s) the block count is no longer decreasing, so the process is stopped.', file=sys.stderr)

                break

        else:

            if verbosity > 3:

                print(f'', file=sys.stderr)

            if verbosity > 2:

                print(f'After {count} iteration(s) of low variance merging the block count is {len(haploblock_dict)}. {variant_dict_pre_filter-sum([len(variant_dict[chrom]) for chrom in variant_dict])} variants were removed.', file=sys.stderr)

            if verbosity > 3:

                print(f'', file=sys.stderr)

    if verbosity > 2:

        print("\n", file=sys.stderr)

    if verbosity > 1:

        print(f'After low variance merging the block count is {len(haploblock_dict)} and {sum([len(variant_dict[chrom]) for chrom in variant_dict])} variants remain usable.\n', file=sys.stderr)

    if verbosity > 4:

        print(len(not_unique_list), sum([read_dict[read_type][read][start_pos][3] == "ungrouped" and len(read_dict[read_type][read][start_pos][2]) != 0 for read_type in read_dict for read in read_dict[read_type] for start_pos in read_dict[read_type][read]]),"\n", file=sys.stderr)
        test_haploblock_dict(haploblock_dict)
        test_not_unique_list(not_unique_list, haploblock_dict)
        test_variant_dict()

    if args.output_intermediate_states:

        make_test_output(f'{args.outdir}/test_outputs/{args.prefix}_secondary_filtered_3_of_5_post_low_variance_merging.bam')

    return haploblock_dict, not_unique_list


def make_hap_dict(haploblock_dict):

    hap_dict = {}

    for block in haploblock_dict:

        read_list = []
        read_variant_dict = {}
        hap_dict[block] = [{}, read_variant_dict, read_list]

        for read_type, read, start_pos in haploblock_dict[block][1]:

            if read_type == "s":

                sec_block = read_dict["p"][read]["-"][3]

                if sec_block not in hap_dict[block][0]:

                    hap_dict[block][0][sec_block] = [0, {}, []]

                hap_dict[block][0][sec_block][0] += 1
                hap_dict[block][0][sec_block][1].update(read_dict[read_type][read][start_pos][2].copy())
                hap_dict[block][0][sec_block][2].append((read_type, read, start_pos))

            else:

                read_list.append((read_type, read, start_pos))
                read_variant_dict.update(read_dict[read_type][read][start_pos][2].copy())

    return hap_dict

def cleanup_hap_dict(hap_dict, block): # Calls remove_read_from_variant_dict

    new_read_list = []
    new_read_dict = {}

    for read_type, read, start_pos in hap_dict[block][2]:

        if read_dict[read_type][read][start_pos][3] not in ["removed", "ungrouped"]:

            new_read_dict.update(read_dict[read_type][read][start_pos][1].copy())
            new_read_list.append((read_type, read, start_pos))

    hap_dict[block] = [hap_dict[block][0], new_read_dict, new_read_list]

    for sec_block in hap_dict[block][0]:

        new_count = 0
        new_read_list = []
        new_read_dict = {}

        for read_type, read, start_pos in hap_dict[block][0][sec_block][2]:

            if read_dict[read_type][read][start_pos][3] not in ["removed", "ungrouped"]:

                new_count += 1
                new_read_dict.update(read_dict[read_type][read][start_pos][1].copy())
                new_read_list.append((read_type, read, start_pos))

            elif read_dict[read_type][read][start_pos][3] == "removed": # FIXME: is this necessary?

                for pos in list(read_dict[read_type][read][start_pos][1].keys()):

                    remove_read_from_variant_dict(pos, read_type, read, start_pos, haploblock_dict)

        hap_dict[block][0][sec_block] = [new_count, new_read_dict, new_read_list]

def choose_primary_or_secondary_reads(hap_dict): # Calls cleanup_hap_dict

    if verbosity > 2:

        print("Choosing between the primary or a secondary alignment for each read.\n", file=sys.stderr)

    merges = []
    compared_pairs = []

    for hap in list(haploblock_dict.keys()):

        block1 = hap

        if len(hap_dict[block1][2]) < minimum_coverage:

            continue

        if verbosity > 2:

            print(f"Evaluating {block1}", file=sys.stderr)

        if verbosity > 4:

            print("\t", block1, len(hap_dict[block1][2]), [(hap, hap_dict[block1][0][hap][0]) for hap in hap_dict[block1][0]], file=sys.stderr)

        for block2 in list(hap_dict[block1][0].keys()):

            if block2 in ["unphasable", "ungrouped"]:

                continue

            block1 = hap

            if f'{block1}_{block2}' in compared_pairs or f'{block2}_{block1}' in compared_pairs:

                continue

            else:

                compared_pairs.append(f'{block1}_{block2}')

            if verbosity > 2:

                if block1 == block2:

                    print(f" Both the primary and secondary reads map to {block1}", file=sys.stderr)

                else:

                    print(f" Comparing {block1} and {block2}", file=sys.stderr)

            if verbosity > 4:

                if block1 == block2:

                    print("\t", block1, len(hap_dict[block1][2]), [(hap, hap_dict[block1][0][hap][0]) for hap in hap_dict[block1][0]])

                else:

                    print("\t", block1, len(hap_dict[block1][2]), [(hap, hap_dict[block1][0][hap][0]) for hap in hap_dict[block1][0]], file=sys.stderr)
                    print("\t", block2, len(hap_dict[block2][2]), [(hap, hap_dict[block2][0][hap][0]) for hap in hap_dict[block2][0]], file=sys.stderr)

            if block1 not in hap_dict[block2][0]:

                hap_dict[block2][0][block1] = [0, {}, []]

            if block1 == block2:

                removed_reads = []

                block_min = min(haploblock_dict[block1][0])
                block_max = max(haploblock_dict[block1][0])

                for _, read, start_pos in hap_dict[block1][0][block1][2]:

                    if min(min(read_dict["p"][read]["-"][2])-block_min, block_max-max(read_dict["p"][read]["-"][2])) < min(min(read_dict["s"][read][start_pos][2])-block_min, block_max-max(read_dict["s"][read][start_pos][2])):

                        removed_reads.append(f'{"p"}_{read}_{"-"}')

                    else:

                        removed_reads.append(f'{"s"}_{read}_{start_pos}')

                new_read_dict = {}
                new_read_list = []

                for read_type, read, start_pos in hap_dict[block1][0][block1][2]:

                    if f'{read_type}_{read}_{start_pos}' in removed_reads:

                        read_dict[read_type][read][start_pos][3] = "ungrouped"

                    else:

                        new_read_dict.update(read_dict[read_type][read][start_pos][2].copy())
                        new_read_list.append((read_type, read, start_pos))

                hap_dict[block1][0][block1][0] = len(new_read_list)
                hap_dict[block1][0][block1][1] = new_read_dict
                hap_dict[block1][0][block1][2] = new_read_list

                if "ungrouped" in hap_dict[block1][0]:

                    for read_type, read, start_pos in hap_dict[block1][0]["ungrouped"][2]:

                        read_dict[read_type][read][start_pos][3] = "ungrouped"

                    hap_dict[block1][0]["ungrouped"] = [0, {}, []]

                new_read_dict = {}
                new_read_list = []

                for read_type, read, start_pos in hap_dict[block1][2]:

                    if f'{read_type}_{read}_{start_pos}' in removed_reads:

                        read_dict[read_type][read][start_pos][3] = "ungrouped"

                    elif read in read_dict["s"] and "ungrouped" in [read_dict["s"][read][start_pos][3] for start_pos in read_dict["s"][read]]: # FIXME: test. This is added because if the secondary read is ungrouped because it ambigously maps to more than one read, there is no way to know if the secondary read is correct. Therefore the primary read is ungrouped.

                        read_dict[read_type][read][start_pos][3] = "ungrouped"

                    else:

                        new_read_dict.update(read_dict[read_type][read][start_pos][2].copy())
                        new_read_list.append((read_type, read, start_pos))

                hap_dict[block1][1] = new_read_dict
                hap_dict[block1][2] = new_read_list

            elif len(hap_dict[block2][2])-hap_dict[block1][0][block2][0] < minimum_coverage or len(hap_dict[block1][2])-hap_dict[block2][0][block1][0] < minimum_coverage:

                if len(hap_dict[block2][2])-hap_dict[block1][0][block2][0] > len(hap_dict[block1][2])-hap_dict[block2][0][block1][0]:

                    block1, block2 = block2, block1

                if verbosity > 2:

                    print(f'  {block2} secondary reads are removed because they map to {block1}')
                    print(f'  {block2} primary reads are removed because their secondary reads map to {block1}', file=sys.stderr)

                if block2 in hap_dict[block1][0]:

                    new_read_list = []
                    new_read_dict = {}

                    for read_type, read, start_pos in hap_dict[block2][2]:

                        found = False

                        for _, read_2, _ in hap_dict[block1][0][block2][2]:

                            if read == read_2:

                                found = True

                                break

                        if not found:

                            new_read_list.append((read_type, read, start_pos))
                            new_read_dict.update(read_dict[read_type][read][start_pos][1].copy())

                        else:

                            read_dict[read_type][read][start_pos][3] = "removed"

                    hap_dict[block2] = [hap_dict[block2][0], new_read_dict, new_read_list]

                if block1 in hap_dict[block2][0]:

                    for read_type, read, start_pos in hap_dict[block2][0][block1][2]:

                        read_dict[read_type][read][start_pos][3] = "removed"

                    hap_dict[block2][0][block1] = [hap_dict[block2][0][block1][0], {}, []]

            else:

                if max(list(hap_dict[block1][1].keys())) < min(list(hap_dict[block2][1].keys())) or max(list(hap_dict[block2][1].keys())) < min(list(hap_dict[block1][1].keys())):

                    merge_status = "no overlap"

                elif min(list(hap_dict[block1][1].keys())) < min(list(hap_dict[block2][1].keys())) and max(list(hap_dict[block2][1].keys())) < max(list(hap_dict[block1][1].keys())) or min(list(hap_dict[block2][1].keys())) < min(list(hap_dict[block1][1].keys())) and max(list(hap_dict[block1][1].keys())) < max(list(hap_dict[block2][1].keys())):

                    merge_status = "within"

                elif (min(list(hap_dict[block1][1].keys())) < min(list(hap_dict[block2][1].keys())) and min(list(hap_dict[block2][1].keys())) < max(list(hap_dict[block1][1].keys())) and max(list(hap_dict[block1][1].keys())) < max(list(hap_dict[block2][1].keys()))) or (min(list(hap_dict[block2][1].keys())) < min(list(hap_dict[block1][1].keys())) and min(list(hap_dict[block1][1].keys())) < max(list(hap_dict[block2][1].keys())) and max(list(hap_dict[block2][1].keys())) < max(list(hap_dict[block1][1].keys()))):

                    merge_status = "overlap"

                else:

                    merge_status = "error"

                if merge_status in ["within", "error"]:

                    if merge_status == "error" and len(haploblock_dict[block2][0]) > len(haploblock_dict[block1][0]) or merge_status == "within" and min(haploblock_dict[block2][0]) < min(haploblock_dict[block1][0]) and max(haploblock_dict[block1][0]) < max(haploblock_dict[block2][0]):

                        block1, block2 = block2, block1

                    read_list = [read for _, read, _ in hap_dict[block1][0][block2][2]]

                    new_read_dict = {}
                    new_read_list = []

                    for read_type, read, start_pos in hap_dict[block2][2]:

                        if read not in read_list:

                            new_read_dict.update(read_dict[read_type][read][start_pos][2].copy())
                            new_read_list.append((read_type, read, start_pos))

                        else:

                            read_dict[read_type][read][start_pos][3] = "ungrouped"

                    for read_type, read, start_pos in hap_dict[block2][0][block1][2]:

                        read_dict[read_type][read][start_pos][3] = "ungrouped"

                    hap_dict[block2][0][block1] = [hap_dict[block2][0][block1][0], {}, []]

                    hap_dict[block2][1] = new_read_dict
                    hap_dict[block2][2] = new_read_list
                
                else:

                    for read_type, read, start_pos in hap_dict[block1][0][block2][2]+hap_dict[block2][0][block1][2]:

                        read_dict[read_type][read][start_pos][3] = "ungrouped"

                # To figure out if two blocks should be merged, we evaluate how many of the possible paired reads mapped to the blocks are shared between them.

                # to make sure that reads are not counted twice, the secondary reads with primary reads in either of the blocks are not counted.
                total_paired_reads = sum([read in read_dict["s"] for _, read, _ in hap_dict[block1][2]+hap_dict[block2][2]]) + sum([hap_dict[block1][0][sub_block][0] if sub_block not in [block1, block2] else 0 for sub_block in hap_dict[block1][0]]) + sum([hap_dict[block2][0][sub_block][0] if sub_block not in [block1, block2] else 0 for sub_block in hap_dict[block2][0]])

                # The total shared are all paired reads that are not mapped to a different block. Therefore reads, where their partner is ungrouped or unphasable will count towards shared reads.
                total_shared_reads = sum([hap_dict[block1][0][sub_block][0] if sub_block in [block2, "ungrouped", "unphasable"] else 0 for sub_block in hap_dict[block1][0]]) + sum([hap_dict[block2][0][sub_block][0] if sub_block in [block1, "ungrouped", "unphasable"] else 0 for sub_block in hap_dict[block2][0]]) + sum([read_type == "s" and read_dict["p"][read]["-"][3] in [block1, block2] for read_type, read, _ in not_unique_list])

                if verbosity > 4:
                    
                    print(total_paired_reads, total_shared_reads, total_shared_reads/total_paired_reads, file=sys.stderr)

                if total_shared_reads/total_paired_reads > minimum_pct_common_reads:

                    if merge_status != "error":

                        merges.append([block1, block2, merge_status, hap_dict[block1][0][block2][0]+hap_dict[block2][0][block1][0]])

                    if verbosity > 2:

                        if merge_status == "error":

                            print(f'  There is seemingly an error.', file=sys.stderr)

                        elif merge_status == "within":

                            if min(list(hap_dict[block1][1].keys())) < min(list(hap_dict[block2][1].keys())) and max(list(hap_dict[block2][1].keys())) < max(list(hap_dict[block1][1].keys())):

                                print(f'  {block1} and {block2} seemingly come from the same haplotype and will therefore eventually be merged. {block2} is {merge_status} {block1}.', file=sys.stderr)

                            elif min(list(hap_dict[block2][1].keys())) < min(list(hap_dict[block1][1].keys())) and max(list(hap_dict[block1][1].keys())) < max(list(hap_dict[block2][1].keys())):

                                print(f'  {block2} and {block1} seemingly come from the same haplotype and will therefore eventually be merged. {block1} is {merge_status} {block2}.', file=sys.stderr)

                        else:

                            if hap_dict[block1][0][block2][0] > 0:

                                print(f'  {block2} secondary reads are ungrouped because they also map to {block1}', file=sys.stderr)

                            if hap_dict[block2][0][block1][0] > 0:

                                print(f'  {block1} secondary reads are ungrouped because they also map to {block2}', file=sys.stderr)

                            if merge_status == "overlap":

                                indel = "insertion"

                            else:

                                indel = "deletion"

                            print(f'  {block1} and {block2} seemingly come from the same haplotype and will therefore eventually be merged. They have {merge_status} and are therefore most likely a large {indel}.', file=sys.stderr)

                else:

                    if verbosity > 2:

                        print(f'  {block1} and {block2} seemingly come from the same haplotype, but only {total_shared_reads}/{total_paired_reads} = {(total_shared_reads/total_paired_reads)*100}% of the possible paired reads connect the two blocks, so they should not be merged.', file=sys.stderr)
                
                hap_dict[block1][0][block2][2] = []
                hap_dict[block2][0][block1][2] = []

            if verbosity > 4:

                if block1 == block2:

                    print("\t", block1, len(hap_dict[block1][2]), [(hap, hap_dict[block1][0][hap][0]) for hap in hap_dict[block1][0]], file=sys.stderr)

                else:

                    print("\t", block1, len(hap_dict[block1][2]), [(hap, hap_dict[block1][0][hap][0]) for hap in hap_dict[block1][0]], file=sys.stderr)
                    print("\t", block2, len(hap_dict[block2][2]), [(hap, hap_dict[block2][0][hap][0]) for hap in hap_dict[block2][0]], file=sys.stderr)

    return hap_dict, merges

def cleanup_low_depth_reads(hap_dict, not_unique_list): # Calls reintroduce_read_to_variant_dict and cleanup_hap_dict

    # Makes sure there are no lone reads

    if verbosity > 2:

        print("\nCleaning up lone reads.", file=sys.stderr)

    for block in list(hap_dict.keys()):

        # Continues if there are no reads attributed to the block
        if len(hap_dict[block][2]) + sum([len(hap_dict[block][0][sec_block][2]) for sec_block in hap_dict[block][0]]) == 0:

            continue

        if verbosity > 2:

            print(f'Evaluating {block}', file=sys.stderr)

        if verbosity > 4:

            print("\t", block, len(hap_dict[block][2]), [(hap, hap_dict[block][0][hap][0]) for hap in hap_dict[block][0]], file=sys.stderr)


        # Makes a read list and read variant dictionary with all reads in the block

        read_variant_dict = hap_dict[block][1].copy()
        read_list = hap_dict[block][2].copy()

        for sec_block in list(hap_dict[block][0].keys()):

            if hap_dict[block][0][sec_block][0] == 0:

                continue

            # Removed secondary reads mapped to another block i that block is below the minimum primary read depth.
            if sec_block != "ungrouped" and len(hap_dict[block][2]) < minimum_coverage:

                if verbosity > 2:

                    print(f'   {hap_dict[block][0][sec_block][0]} reads were removed because their primary read is in {sec_block}.', file=sys.stderr)

                for read_type, read, start_pos in hap_dict[block][0][sec_block][2]:

                    read_dict[read_type][read][start_pos][3] = "removed"

                hap_dict[block][0][sec_block] = [{}, {}, []]

                continue

            read_variant_dict.update(hap_dict[block][0][sec_block][1].copy())
            read_list += hap_dict[block][0][sec_block][2]

        # If the primary read depth is small then the reads are moved to ungrouped.

        if len(hap_dict[block][2]) < minimum_coverage:

            if verbosity > 2:

                print(f'   {len(read_list)} reads were classed as removed because they belong to a removed block.', file=sys.stderr)

            for read_type, read, start_pos in read_list:

                reintroduce_read_to_variant_dict(read_type, read, start_pos)

                read_dict[read_type][read][start_pos][3] = "removed"

            # Makes it empty so that the block is deleted in a later step

            hap_dict[block] = [{}, {}, []]

            if verbosity > 4:

                print("\t", block, len(hap_dict[block][2]), [(hap, hap_dict[block][0][hap][0]) for hap in hap_dict[block][0]], file=sys.stderr)

            continue

        reads_to_remove = []
        old_count_list = []

        for chrom in variant_dict:
            for pos in variant_dict[chrom]:

                if pos not in read_variant_dict:

                    continue

                # Checks if there are more than the minimal read count mapped to that position
                if sum([pos in read_dict[read_type][read][start_pos][1] for read_type, read, start_pos in read_list]) >= minimum_depth:

                    continue

                # If less reads than the minimal read count is mapped to the position, then the following section figures out which reads this is.

                count_list = []

                for read_type, read, start_pos in read_list:

                    if sum([(read_type, read, start_pos) == read_vec for read_vec in old_count_list]) > 0:

                        continue # FIXME: føles lidt som dobbelkonfekt...

                    if pos in read_dict[read_type][read][start_pos][1]:

                        count_list.append((read_type, read, start_pos))

                old_count_list += count_list

                # This would mean that all reads in this position has been evaluated
                if len(count_list) == 0:

                    continue

                # Checks if the read, which overlaps with a position the less than the minimum depth has other postitions which is above the minimum depth
                for read_type, read, start_pos in count_list:

                    min_pos = min(read_dict[read_type][read][start_pos][1])
                    max_pos = max(read_dict[read_type][read][start_pos][1])

                    min_count = 0
                    max_count = 0

                    for read_type_2, read_2, start_pos_2 in read_list:

                        if min_pos in read_dict[read_type_2][read_2][start_pos_2][1]:

                            min_count += 1

                            if min_count == minimum_depth:

                                break

                        if max_pos in read_dict[read_type_2][read_2][start_pos_2][1]:

                            max_count += 1

                            if max_count == minimum_depth:

                                break

                    if min_count == minimum_depth or max_count == minimum_depth:

                        continue

                    else:

                        reads_to_remove.append((read_type, read, start_pos))

        if len(reads_to_remove) == 0:

            if verbosity > 3:

                print(f'   No reads cover a region with a depth less than {minimum_depth}.', file=sys.stderr)

            continue

        elif verbosity > 2:

            print(f'   {len(reads_to_remove)} reads cover a region with a depth less than {minimum_depth}. These reads are therefore moved to the ungrouped list.', file=sys.stderr)

        if verbosity > 3:

            print(f'   The following reads were removed from {block}: {reads_to_remove}', file=sys.stderr)

        # The following moves the low depth reads to ungrouped.

        for read_type, read, start_pos in reads_to_remove:

            if read_type == "p" and read not in read_dict["s"]:

                reintroduce_read_to_variant_dict(read_type, read, start_pos)

                not_unique_list.append((read_type, read, start_pos))

            else:

                read_dict["p"][read]["-"][3] = "ungrouped"
                not_unique_list.append(("p", read, "-"))

                for start_pos in read_dict["s"][read]:

                    read_dict["s"][read][start_pos][3] = "ungrouped"
                    not_unique_list.append(("s", read, start_pos))

        if verbosity > 4:

            print("\t", block, len(hap_dict[block][2]), [(hap, hap_dict[block][0][hap][0]) for hap in hap_dict[block][0]], file=sys.stderr)

def update_block_info(block, hap_dict, haploblock_dict):
    
    new_read_dict = hap_dict[block][1].copy()
    new_read_list = hap_dict[block][2].copy()

    for sec_block in hap_dict[block][0]:

        new_read_dict.update(hap_dict[block][0][sec_block][1].copy())
        new_read_list += hap_dict[block][0][sec_block][2]

    if block not in haploblock_dict:

        haploblock_dict[block] = [None, None]

    haploblock_dict[block][0] = new_read_dict
    haploblock_dict[block][1] = new_read_list

    return new_read_dict, new_read_list

def update_haploblock_dict(hap_dict, merges, haploblock_dict):

    if verbosity > 2:

        print("\nUpdating haploblock dict.\n", file=sys.stderr)

    for block in hap_dict:

        if len(hap_dict[block][2]) == 0:

            if verbosity > 2:

                print(f'{block} is removed.', file=sys.stderr)

            if block in haploblock_dict:

                del haploblock_dict[block]

            continue

        update_block_info(block, hap_dict, haploblock_dict)

        if verbosity > 2:

            print(f'{block} is updated.', file=sys.stderr)

    if verbosity > 4:
        
        print(merges, file=sys.stderr)

    i = 0

    while i < len(merges):

        removed = False

        block1 = merges[i][0]
        block2 = merges[i][1]
        merge_status = merges[i][2]

        if block1 not in haploblock_dict or block2 not in haploblock_dict:

            del merges[i]

            continue

        if merge_status == "within": # FIXME: should it not be handled somehow?

            i += 1

            continue

        if merge_status == "overlap":

            # This section captures the situation where bad merging due to a homogeneous region has merged secondary reads where they do not belong. 

            overlap, identity, mismatch_list = compare_haploblocks(block1, block2, haploblock_dict)

            if sum([base1 not in ["I", "D", "no overlap", "overlap", "within"] and base2 not in ["I", "D", "no overlap", "overlap", "within"] and sum([pos in read_dict[read_type][read][start_pos][2] for read_type, read, start_pos in haploblock_dict[block1][1]]) > minimum_coverage and sum([pos in read_dict[read_type][read][start_pos][2] for read_type, read, start_pos in haploblock_dict[block2][1]]) > minimum_coverage for pos, base1, base2 in mismatch_list]) == 0:

                if verbosity > 3:

                    print(f"   {block1} and {block2} is merged in haploblock updating. Overlap is {overlap} and mismatch count is {len(mismatch_list)}.", file=sys.stderr)

                merge_haploblocks(block1, block2, mismatch_list, haploblock_dict)

                del merges[i]

                continue

        for j in range(i+1, len(merges)):

            block1 = merges[j][0]
            block2 = merges[j][1]
            merge_status = merges[j][2]

            if block1 not in haploblock_dict or block2 not in haploblock_dict:

                if block1 in haploblock_dict:

                    for read_type, read, start_pos in hap_dict[block1][0][block2][2]:

                        read_dict[read_type][read][start_pos][3] = "ungrouped"

                    del hap_dict[block1][0][block2]

                    update_block_info(block1, hap_dict, haploblock_dict)

                if block2 in haploblock_dict:

                    for read_type, read, start_pos in hap_dict[block2][0][block1][2]:

                        read_dict[read_type][read][start_pos][3] = "ungrouped"

                    del hap_dict[block2][0][block1]

                    update_block_info(block2, hap_dict, haploblock_dict)

                del merges[j]

                break

            if merge_status == "within":

                continue # FIXME: should it not be handled somehow?

            if block1 in merges[i] or block2 in merges[i]:

                removed = True

                if merges[i][3] < merges[j][3]:

                    merges[i] = merges[j]

                block1, block2, _, _ = merges[j]

                for read_type, read, start_pos in hap_dict[block1][0][block2][2] + hap_dict[block2][0][block1][2]:

                    read_dict[read_type][read][start_pos][3] = "ungrouped"

                hap_dict[block1][0][block2] = [0, {}, []]
                hap_dict[block2][0][block1] = [0, {}, []]

                update_block_info(block1, hap_dict, haploblock_dict)
                update_block_info(block2, hap_dict, haploblock_dict)

                del merges[j]

                break

        if not removed:

            i += 1

    if verbosity > 4:
        
        print(merges, file=sys.stderr)

    merge_list = []

    for hap1, hap2, merge_status, _ in merges:

        if hap1 not in haploblock_dict or hap2 not in haploblock_dict:

            continue

        if min(haploblock_dict[hap1][0]) < min(haploblock_dict[hap2][0]) and max(haploblock_dict[hap2][0]) < max(haploblock_dict[hap1][0]) or min(haploblock_dict[hap1][0]) < min(haploblock_dict[hap2][0]) and max(haploblock_dict[hap1][0]) < max(haploblock_dict[hap2][0]):

            merge_list.append([hap1, hap2, merge_status])

        elif min(haploblock_dict[hap2][0]) < min(haploblock_dict[hap1][0]) and max(haploblock_dict[hap1][0]) < max(haploblock_dict[hap2][0]) or min(haploblock_dict[hap2][0]) < min(haploblock_dict[hap1][0]) and max(haploblock_dict[hap2][0]) < max(haploblock_dict[hap1][0]):

            merge_list.append([hap2, hap1, merge_status])

        else:

            print("Something has gone wrong in update_haploblock_dict")

    if verbosity > 4:
        
        print(merge_list, file=sys.stderr)

    if verbosity > 2:

        print("", file=sys.stderr)

    return merge_list

def update_merge_list(hap1, hap2, haploblock_dict, merged_list):

    if verbosity > 4:

        print("Before merge updating:", merged_list, file=sys.stderr)

    i = 0

    while i < len(merged_list):

        merged = False

        merge_vec = merged_list[i]

        if hap1 in merge_vec and hap2 in merge_vec:

            del merged_list[i]

            continue

        elif hap2 in merge_vec:

            if hap2 == merge_vec[0]:

                merged_list[i][0] = hap1

            elif hap2 == merge_vec[1]:

                merged_list[i][1] = hap1

            merged = True

        elif hap1 in merge_vec:

            merged = True

        if merged:

            block1 = merged_list[i][0]
            block2 = merged_list[i][1]

            if max(list(haploblock_dict[block1][0].keys())) < min(list(haploblock_dict[block2][0].keys())) or max(list(haploblock_dict[block2][0].keys())) < min(list(haploblock_dict[block1][0].keys())):

                merge_status = "no overlap"

                if merged_list[i][2] == "within" and max(list(haploblock_dict[block2][0].keys())) < min(list(haploblock_dict[block1][0].keys())):

                    merged_list[i][0], merged_list[i][1] = merged_list[i][1], merged_list[i][0]

            elif min(list(haploblock_dict[block1][0].keys())) < min(list(haploblock_dict[block2][0].keys())) and max(list(haploblock_dict[block2][0].keys())) < max(list(haploblock_dict[block1][0].keys())) or min(list(haploblock_dict[block2][0].keys())) < min(list(haploblock_dict[block1][0].keys())) and max(list(haploblock_dict[block1][0].keys())) < max(list(haploblock_dict[block2][0].keys())):

                merge_status = "within"

                if min(list(haploblock_dict[block2][0].keys())) < min(list(haploblock_dict[block1][0].keys())) and max(list(haploblock_dict[block1][0].keys())) < max(list(haploblock_dict[block2][0].keys())):

                    merged_list[i][0], merged_list[i][1] = merged_list[i][1], merged_list[i][0]

            elif (min(list(haploblock_dict[block1][0].keys())) < min(list(haploblock_dict[block2][0].keys())) and min(list(haploblock_dict[block2][0].keys())) < max(list(haploblock_dict[block1][0].keys())) and max(list(haploblock_dict[block1][0].keys())) < max(list(haploblock_dict[block2][0].keys()))) or (min(list(haploblock_dict[block2][0].keys())) < min(list(haploblock_dict[block1][0].keys())) and min(list(haploblock_dict[block1][0].keys())) < max(list(haploblock_dict[block2][0].keys())) and max(list(haploblock_dict[block2][0].keys())) < max(list(haploblock_dict[block1][0].keys()))):

                merge_status = "overlap"

                if merged_list[i][2] == "within" and (min(list(haploblock_dict[block2][0].keys())) < min(list(haploblock_dict[block1][0].keys())) and min(list(haploblock_dict[block1][0].keys())) < max(list(haploblock_dict[block2][0].keys())) and max(list(haploblock_dict[block2][0].keys())) < max(list(haploblock_dict[block1][0].keys()))):

                    merged_list[i][0], merged_list[i][1] = merged_list[i][1], merged_list[i][0]

            else:

                merge_status = "error"

            if merge_status == "error":

                pass

            else:

                merged_list[i][2] = merge_status

        i += 1

    merged_list = set([",".join(merge) for merge in merged_list])

    merged_list = [merge.split(",") for merge in merged_list]

    if verbosity > 4:

        print("After merge updating:", merged_list, file=sys.stderr)

    return merged_list

def stitch_blocks(haploblock_dict, maxdiff, merged_list): # Calls compare_haploblocks, merge_haploblocks and update_merge_list

    if verbosity > 3:

        print(f'', file=sys.stderr)

    if verbosity > 2:

        print("Merging blocks.", file=sys.stderr)

    if verbosity > 3:

        print(f'', file=sys.stderr)

    hap_names = list(haploblock_dict.keys())

    for i in range(len(hap_names)):

        hap = hap_names[i]

        for j in range(i+1, len(hap_names)):

            hap1 = hap

            hap2 = hap_names[j]

            if hap not in haploblock_dict or hap2 not in haploblock_dict:

                continue

            overlap, identity, mismatch_list = compare_haploblocks(hap1, hap2, haploblock_dict)

            if verbosity > 4 and overlap > 0:

                print(hap1, hap2, overlap, identity, len(mismatch_list), file=sys.stderr)

            if overlap > 0 and len(mismatch_list) <= maxdiff and (minimum_identity <= identity/overlap or sum([base1 not in ["I", "D", "no overlap", "overlap", "within"] and base2 not in ["I", "D", "no overlap", "overlap", "within"] and sum([pos in read_dict[read_type][read][start_pos][2] for read_type, read, start_pos in haploblock_dict[hap1][1]]) > minimum_coverage and sum([pos in read_dict[read_type][read][start_pos][2] for read_type, read, start_pos in haploblock_dict[hap2][1]]) > minimum_coverage for pos, base1, base2 in mismatch_list]) == 0):

                if len(mismatch_list) > 0:

                    # print(mismatch_list, sum([sum([pos in haploblock_dict[hap][0] and len(haploblock_dict[hap][0][pos]) == 1 for hap in haploblock_dict]) <= expected_polyploidy for pos, _, _ in mismatch_list]))

                    if sum([sum([pos in haploblock_dict[hap][0] and len(haploblock_dict[hap][0][pos]) == 1 for hap in haploblock_dict]) <= expected_polyploidy for pos, _, _ in mismatch_list]): # Checks if the merging would result fewer haplotypes than what is expected. If yes, the merge is not performed. # TODO: not implemented probably

                        continue

                if verbosity > 3:

                    print(f'   {hap1} and {hap2} is merged, because they had an overlap of {overlap} and a mismatch count of {len(mismatch_list)}.', file=sys.stderr)

                merge_haploblocks(hap1, hap2, mismatch_list, haploblock_dict)
                merged_list = update_merge_list(hap1, hap2, haploblock_dict, merged_list)

    return merged_list

def add_to_fit_lists(read_type, read, start_pos, fit_lists, maxdiff, haploblock_dict): # Calls make_fit_list and remove_read_pos_from_variant_dict

    read_variant_dict = read_dict[read_type][read][start_pos][2]

    if maxdiff < len(read_variant_dict)*minimum_pct_overlap:

        fit_list = make_fit_list(read_variant_dict, read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, maxdiff, len(read_variant_dict)*minimum_pct_overlap)

    else:

        fit_list = make_fit_list(read_variant_dict, read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, len(read_variant_dict)*minimum_pct_overlap, len(read_variant_dict)*minimum_pct_overlap)

    if len(fit_list) > 0:

        fit_lists.append((start_pos, fit_list))

        if len(fit_list) == 1 and len(fit_list[0][3]) > 0:

            if verbosity > 3:

                print(f'The following variants were removed from {read_type} {start_pos}: {fit_list[0][3]} ')

            pre_filter = read_dict[read_type][read][start_pos][3]

            remove_read_pos_from_variant_dict(fit_list[0][3], read_type, read, start_pos, haploblock_dict)

            read_dict[read_type][read][start_pos][3] = pre_filter

def reevaluate_ungrouped_primary_and_secondary_reads(haploblock_dict, merged_list): # Calls reintroduce_read_to_variant_dict and add_to_fit_lists

    if verbosity > 3:

        print(f'', file=sys.stderr)

    if verbosity > 2:

        print("Removing ungrouped reads with a primary or secondary read in a block.", file=sys.stderr)

    if verbosity > 3:

        print(f'', file=sys.stderr)

    unchosen_list = []
    new_not_unique_list = []

    for read in read_dict["p"]:

        if read in read_dict["s"]:

            if verbosity > 4:

                print(read, (len(read_dict["p"][read]["-"][2]), read_dict["p"][read]["-"][3]), [(len(read_dict["s"][read][start_pos][2]), read_dict["s"][read][start_pos][3]) for start_pos in read_dict["s"][read]], file=sys.stderr)

            groups = set([read_dict["p"][read]["-"][3]] + [read_dict["s"][read][start_pos][3] for start_pos in read_dict["s"][read]])

            # Makes sure all reads have usable SNPs. If one doesn't it could give a bad result.

            if len(read_dict["p"][read]["-"][2]) == 0 or sum([len(read_dict["s"][read][start_pos][2]) == 0 for start_pos in read_dict["s"][read]]) > 0:

                if read_dict["p"][read]["-"][3] != "removed":

                    reintroduce_read_to_variant_dict("p", read, "-")

                    if len([read_dict["p"][read]["-"][3]] + [read_dict["s"][read][start_pos][3] for start_pos in read_dict["s"][read]]) == 2 and "removed" in groups:

                        if len(read_dict["p"][read]["-"][2]) == 0:

                            read_dict["p"][read]["-"][3] = "unphasable"

                        else:

                            new_not_unique_list.append(("p", read, "-"))

                    else:

                        read_dict["p"][read]["-"][3] = "unchosen"

                for start_pos in read_dict["s"][read]:

                    if read_dict["s"][read][start_pos][3] != "removed":

                        reintroduce_read_to_variant_dict("s", read, start_pos)

                        if len([read_dict["p"][read]["-"][3]] + [read_dict["s"][read][start_pos][3] for start_pos in read_dict["s"][read]]) == 2 and "removed" in groups:

                            if len(read_dict["s"][read][start_pos][2]) == 0:

                                read_dict["s"][read][start_pos][3] = "unphasable"

                            else:

                                new_not_unique_list.append(("s", read, start_pos))

                        else:

                            read_dict["s"][read][start_pos][3] = "unchosen"

                groups = set([read_dict["p"][read]["-"][3]] + [read_dict["s"][read][start_pos][3] for start_pos in read_dict["s"][read]])

                if "unchosen" in groups:

                    unchosen_list.append(read)

            elif "removed" in groups:

                possible_reads = []

                if read_dict["p"][read]["-"][3] != "removed":

                    possible_reads.append(("p", read, "-"))

                for start_pos in read_dict["s"][read]:

                    if read_dict["s"][read][start_pos][3] != "removed":

                        possible_reads.append(("s", read, start_pos))

                if len(possible_reads) == 1:

                    read_type, read, start_pos = possible_reads[0]

                    read_variant_dict = read_dict[read_type][read][start_pos][2]

                    fit_list = make_fit_list(read_variant_dict, read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, 0, 1)

                    if verbosity > 4:

                        print(fit_list, file=sys.stderr)

                    if len(fit_list) == 1:

                        block = fit_list[0][0]

                        if read_dict[read_type][read][start_pos][3] == "ungrouped":

                            reintroduce_read_to_variant_dict(read_type, read, start_pos)
                            add_read_to_haploblock(read_type, read, start_pos, block, fit_list[0][3], haploblock_dict)

                        elif read_dict[read_type][read][start_pos][3] != block: # FIXME: test

                            print("If this prints, something is wrong in reevaluate_ungrouped_primary_and_secondary_reads.")

                    else:

                        reintroduce_read_to_variant_dict(read_type, read, start_pos)
                        new_not_unique_list.append((read_type, read, start_pos))

                else:

                    # FIXME: not sure what we might miss here... it is not very common to have more than 2 mappings

                    for read_type, read, start_pos in possible_reads:

                        read_dict[read_type][read][start_pos][3] = "unchosen"

            else:

                fit_lists = []

                read_type = "p"
                start_pos = "-"

                add_to_fit_lists(read_type, read, start_pos, fit_lists, 0, haploblock_dict)

                for start_pos in read_dict["s"][read]:

                    read_type = "s"

                    add_to_fit_lists(read_type, read, start_pos, fit_lists, 0, haploblock_dict)

                if verbosity > 4:

                    print(fit_lists, file=sys.stderr)

                if len(fit_lists) == 1 and len(fit_lists[0][1]) == 1:

                    start_pos = fit_lists[0][0]

                    read_type = "p"

                    if start_pos == "-":

                        block = fit_lists[0][1][0][0]

                        if block != read_dict[read_type][read][start_pos][3]:

                            reintroduce_read_to_variant_dict(read_type, read, start_pos)
                            add_read_to_haploblock(read_type, read, start_pos, block, fit_lists[0][1][0][3], haploblock_dict)

                    else:

                        read_dict[read_type][read]["-"][3] = "removed"

                    read_type = "s"

                    for start_pos_2 in read_dict[read_type][read]:

                        if start_pos == start_pos_2:

                            block = fit_lists[0][1][0][0]

                            if read_dict[read_type][read][start_pos][3] != block:

                                reintroduce_read_to_variant_dict(read_type, read, start_pos)
                                add_read_to_haploblock(read_type, read, start_pos, block, fit_lists[0][1][0][3], haploblock_dict)

                        else:

                            read_dict[read_type][read][start_pos_2][3] = "removed"

                elif len(fit_lists) == 1:

                    if fit_lists[0][0] == "-":

                        for start_pos_2 in read_dict["s"][read]:

                            read_dict["s"][read][start_pos_2][3] = "removed"

                        read_dict["p"][read]["-"][3] = "ungrouped"

                        reintroduce_read_to_variant_dict("p", read, "-")
                        new_not_unique_list.append(("p", read, "-"))

                    else:

                        for start_pos_2 in read_dict["s"][read]:

                            if start_pos_2 != fit_lists[0][0]:

                                read_dict["s"][read][start_pos_2][3] = "removed"

                            else:

                                read_dict["s"][read][start_pos_2][3] = "ungrouped"

                                reintroduce_read_to_variant_dict("s", read, start_pos_2)
                                new_not_unique_list.append(("s", read, start_pos_2))

                        read_dict["p"][read]["-"][3] = "removed"
                
                # FIXME: Added because I had a problem with multimapped reads that should be kept are disappearing... Seeing if this solves the problem.
                elif len(fit_lists) == 2 and sum([len(fit_lists[i][1]) == 1 for i in range(len(fit_lists))]) == 2 and sum([sorted([fit_lists[0][1][0][0], fit_lists[1][1][0][0]]) == sorted(merged_list[i][:2]) and merged_list[i][2] != "within" for i in range(len(merged_list))]) == 1:

                    # FIXME: This means that the multimapped read is only mapping to the two sides of blocks that should be merged. At the same time they are not deemed to be two blocks, where one is within, because that may cause problems. It should be researched if this is in fact true or just an incorrect assumption.

                    for i in range(len(fit_lists)):

                        if fit_lists[i][0] == "-":

                            read_type = "p"
                            start_pos = fit_lists[i][0]

                        else:

                            read_type = "s"
                            start_pos = fit_lists[i][0]

                        if read_dict[read_type][read][start_pos][3] == "ungrouped":

                            block = fit_lists[i][1][0][0]

                            reintroduce_read_to_variant_dict(read_type, read, start_pos)
                            add_read_to_haploblock(read_type, read, start_pos, block, fit_lists[i][1][0][3], haploblock_dict)

                # elif sum([len(fit_lists[i][1]) == 1 for i in range(len(fit_lists))]) == 2 and fit_lists[0][1][0][0] != fit_lists[1][1][0][0]: # FIXME: created a problem with a within group, which stays

                #     for i in range(len(fit_lists)):

                #         if fit_lists[i][0] == "-":

                #             read_type = "p"
                #             start_pos = fit_lists[i][0]

                #         else:

                #             read_type = "s"
                #             start_pos = fit_lists[i][0]

                #         if read_dict[read_type][read][start_pos][3] == "ungrouped":

                #             block = fit_lists[i][1][0][0]

                #             reintroduce_read_to_variant_dict(read_type, read, start_pos)
                #             add_read_to_haploblock(read_type, read, start_pos, block, fit_lists[i][1][0][3], haploblock_dict)

                else:

                    read_dict["p"][read]["-"][3] = "unchosen"

                    for start_pos_2 in read_dict["s"][read]:

                        read_dict["s"][read][start_pos_2][3] = "unchosen"

                    unchosen_list.append(read)

            groups = set([read_dict["p"][read]["-"][3]] + [read_dict["s"][read][start_pos][3] for start_pos in read_dict["s"][read]])

            if verbosity > 4:

                print([read_dict["p"][read]["-"][3]], [read_dict["s"][read][start_pos][3] for start_pos in read_dict["s"][read]], file=sys.stderr)

                print(f'', file=sys.stderr)

        elif read_dict["p"][read]["-"][3] in ["ungrouped", "removed"]:

            reintroduce_read_to_variant_dict("p", read, "-")
            new_not_unique_list.append(("p", read, "-"))

    print([(hap, set([read_dict[read_type][read][start_pos][3] for read_type, read, start_pos in haploblock_dict[hap][1]]), sum([read_dict[read_type][read][start_pos][3] == "unchosen" for read_type, read, start_pos in haploblock_dict[hap][1]]), sum([read_dict[read_type][read][start_pos][3] == hap for read_type, read, start_pos in haploblock_dict[hap][1]]), len(haploblock_dict[hap][1])) for hap in haploblock_dict], file=sys.stderr)

    for block in list(haploblock_dict.keys()):

        if sum([read_dict[read_type][read][start_pos][3] == block for read_type, read, start_pos in haploblock_dict[block][1]]) < minimum_coverage:

            if verbosity > 3:

                print(f'{block} is removed due to low coverage', file=sys.stderr)

            for read_type, read, start_pos in haploblock_dict[block][1]:

                if read_dict[read_type][read][start_pos][3] == block:

                    new_not_unique_list.append((read_type, read, start_pos))

            del haploblock_dict[block]

            merged_list = update_merge_list(block, block, haploblock_dict, merged_list)

    return unchosen_list, new_not_unique_list, merged_list

def update_merges(merge_list, haploblock_dict):

    if verbosity > 3:

        print(f'', file=sys.stderr)

    if verbosity > 2:

        print("Modify the blocks that should be merged.", file=sys.stderr)

    if verbosity > 3:

        print(f'', file=sys.stderr)

    for i in range(len(merge_list)):

        merge_status = merge_list[i][-1]

        chrom = haploblock_dict[merge_list[i][0]][2]

        if merge_status == "within":

            within_hap = merge_list[i][1]

            within_max = max(haploblock_dict[within_hap][0])
            within_min = min(haploblock_dict[within_hap][0])

            for pos in variant_dict[chrom]:

                if pos < within_min or within_max < pos:

                    haploblock_dict[within_hap][0][pos] = merge_status

        else:

            left_hap = merge_list[i][0]
            right_hap = merge_list[i][-2]

            for j in range(1, len(merge_list[i])-2):

                middle_hap = merge_list[i][j]

                middle_max = max(haploblock_dict[middle_hap][0])
                middle_min = min(haploblock_dict[middle_hap][0])

                for pos in variant_dict[chrom]:

                    if pos < middle_min or middle_max < pos:

                        haploblock_dict[middle_hap][0][pos] = merge_status

            left_max = max(haploblock_dict[merge_list[i][0]][0])
            right_min = min(haploblock_dict[merge_list[i][-2]][0])

            for pos in variant_dict[chrom]:

                if left_max < pos:

                    haploblock_dict[left_hap][0][pos] = merge_status

                if pos < right_min:

                    haploblock_dict[right_hap][0][pos] = merge_status

    return merge_list

def evaluate_unchosen_list(unchosen_list, haploblock_dict, not_unique_list): # Calls add_to_fit_lists

    if verbosity > 3:

        print(f'', file=sys.stderr)

    if verbosity > 2:

        print("Reevaluate unchosen reads.", file=sys.stderr)

    if verbosity > 3:

        print(f'', file=sys.stderr)

    new_unchosen_list = []

    for read in unchosen_list:

        if verbosity > 4:

            print(read, [(read_dict["p"][read]["-"][3],len(read_dict["p"][read]["-"][2]))], [(read_dict["s"][read][start_pos][3], len(read_dict["s"][read][start_pos][2])) for start_pos in read_dict["s"][read]],  file=sys.stderr)

        if len(read_dict["p"][read]["-"][2]) == 0 or sum([len(read_dict["s"][read][start_pos][2]) == 0 for start_pos in read_dict["s"][read]]) > 0:

            continue

        diff = 0

        while True:

            fit_lists = []

            read_type = "p"
            start_pos = "-"

            add_to_fit_lists(read_type, read, start_pos, fit_lists, diff, haploblock_dict)

            for start_pos in read_dict["s"][read]:

                read_type = "s"

                add_to_fit_lists(read_type, read, start_pos, fit_lists, diff, haploblock_dict)

            if len(fit_lists) > 0 or diff == maxdiff:

                break

            else:

                diff += 1
        
        if verbosity > 4:

            print("Evaluate unchosen reads:", fit_lists, file=sys.stderr)

        if len(fit_lists) == 1 and len(fit_lists[0][1]) == 1:

            start_pos = fit_lists[0][0]

            if start_pos == "-":

                block = fit_lists[0][1][0][0]

                haploblock_dict[block][0].update(read_dict["p"][read][start_pos][2].copy())
                haploblock_dict[block][1].append(("p",read,start_pos))
                read_dict["p"][read][start_pos][3] = block

            else:

                read_dict["p"][read]["-"][3] = "removed"

            for start_pos_2 in read_dict["s"][read]:

                if start_pos == start_pos_2:

                    block = fit_lists[0][1][0][0]

                    haploblock_dict[block][0].update(read_dict["s"][read][start_pos][2].copy())
                    haploblock_dict[block][1].append(("s",read,start_pos))
                    read_dict["s"][read][start_pos][3] = block

                else:

                    read_dict["s"][read][start_pos_2][3] = "removed"

        elif len(fit_lists) == 1:

            start_pos = fit_lists[0][0]

            if start_pos == "-":

                for start_pos_2 in read_dict["s"][read]:

                    read_dict["s"][read][start_pos_2][3] = "removed"

                read_dict["p"][read]["-"][3] = "ungrouped"
                not_unique_list.append(("p",read,start_pos))

            else:

                read_dict["p"][read]["-"][3] = "removed"

                for start_pos_2 in read_dict["s"][read]:

                    if start_pos != start_pos_2:

                        read_dict["s"][read][start_pos_2][3] = "removed"

                    else:

                        read_dict["s"][read][start_pos_2][3] = "ungrouped"
                        not_unique_list.append(("s",read,start_pos))
            
        else:

            new_unchosen_list.append(read)

        if verbosity > 4:

            print([read_dict["p"][read]["-"][3]], [read_dict["s"][read][start_pos][3] for start_pos in read_dict["s"][read]], file=sys.stderr)

            print("", file=sys.stderr)

    return new_unchosen_list

def remove_removed_read_variants(haploblock_dict): # Calls remove_variant and cleanup_variants

    if verbosity > 3:

        print(f'', file=sys.stderr)

    if verbosity > 2:

        print(f'Removing unusable variants after read selection', file=sys.stderr)

    if verbosity > 3:

        print(f'', file=sys.stderr)

    for chrom in variant_dict:
        for pos in list(variant_dict[chrom].keys()):

            new_dict = {}

            for base in variant_dict[chrom][pos]:

                new_dict[base] = {}

                for read_type, read, start_pos in variant_dict[chrom][pos][base]:

                    if read_dict[read_type][read][start_pos][3] not in new_dict[base]:

                        new_dict[base][read_dict[read_type][read][start_pos][3]] = 0

                    new_dict[base][read_dict[read_type][read][start_pos][3]] += 1

            if verbosity > 4:
                
                print(chrom, pos, [(base, len(variant_dict[chrom][pos][base])) for base in variant_dict[chrom][pos]], new_dict, file=sys.stderr)

            remove_count = 0

            for base in new_dict:

                if sum([group in ["removed"] or new_dict[base][group] < minimum_count for group in new_dict[base]]) == len(new_dict[base]):

                    remove_count += 1

            if len(new_dict)-1 <= remove_count or len(new_dict)-1 <= sum([sum(list(new_dict[base].values())) < minimum_coverage for base in new_dict]):

                remove_variant(chrom, pos, haploblock_dict, read_dict)

                if verbosity > 3:

                    print(f'   {pos} is removed since there are no variation in the remaining reads.',  file=sys.stderr)

    cleanup_variants(haploblock_dict)

    for block in list(haploblock_dict.keys()):

        if len(haploblock_dict[block][0]) == 0:

            for read_type, read, start_pos in haploblock_dict[block][1]:

                read_dict[read_type][read][start_pos][3] = "unphasable"

            del haploblock_dict[block]

            if verbosity > 3:

                print(f"   {block} is removed and all reads are ungrouped.",  file=sys.stderr)

    new_ungrouped_list = []

    for read_type in read_dict:
        for read in read_dict[read_type]:
            for start_pos in read_dict[read_type][read]:

                if read_dict[read_type][read][start_pos][3] == "ungrouped":

                    new_ungrouped_list.append((read_type, read, start_pos))

    return new_ungrouped_list

def reevaluating_ungrouped_reads(not_unique_list):

    if verbosity > 3:

        print(f'', file=sys.stderr)

    if verbosity > 2:

        print("Reevaluating ungrouped reads.", file=sys.stderr)

    if verbosity > 3:

        print(f'', file=sys.stderr)

    diff = 0

    while diff <= maxdiff and len(not_unique_list) > 0:

        new_not_unique_list = []

        for read_type, read, start_pos in not_unique_list:

            read_dict[read_type][read][start_pos][1] = read_dict[read_type][read][start_pos][2].copy()

            read_variant_dict = read_dict[read_type][read][start_pos][2]

            if len(read_variant_dict) == 0:

                read_dict[read_type][read][start_pos][3] = "unphasable"

                continue

            fit_list = make_fit_list(read_variant_dict, read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, diff, len(read_variant_dict)*minimum_pct_overlap)

            if verbosity > 4:

                print("This is from the remaining ungrouped reads:", fit_list, file=sys.stderr)

            if len(fit_list) == 1:

                block = fit_list[0][0]

                if verbosity > 3:

                    print(f'   Read added to {block}. Overlap is {fit_list[0][1]} and the mismatch count is {len(fit_list[0][3])}', file=sys.stderr)

                add_read_to_haploblock(read_type, read, start_pos, block, fit_list[0][3], haploblock_dict)

            else:

                new_not_unique_list.append((read_type, read, start_pos))

        if len(not_unique_list) == len(new_not_unique_list):

            diff += 1

        not_unique_list = new_not_unique_list

    return not_unique_list

def reevaluate_secondary_reads(not_unique_list, haploblock_dict): # Calls make_hap_dict, choose_primary_or_secondary_reads, cleanup_low_depth_reads, update_haploblock_dict, stitch_blocks, reevaluate_ungrouped_primary_and_secondary_reads, update_merges, reevaluate_remaining_ungrouped_reads, evaluate_unchosen_list, remove_removed_read_variants

    print_subprocess_title("Performing multiple mapped read selection.")

    hap_dict = make_hap_dict(haploblock_dict)

    hap_dict, merges = choose_primary_or_secondary_reads(hap_dict)

    # Makes sure there are no lone reads

    cleanup_low_depth_reads(hap_dict, not_unique_list)

    # Update the haploblock_dict with the hap_dict information

    merges = update_haploblock_dict(hap_dict, merges, haploblock_dict)

    merges = stitch_blocks(haploblock_dict, maxdiff, merges)

    unchosen_list, not_unique_list, merges = reevaluate_ungrouped_primary_and_secondary_reads(haploblock_dict, merges)

    merge_list = update_merges(merges, haploblock_dict)

    while True:

        haploblock_dict_pre_filter = len(haploblock_dict)
        ungrouped_list_pre_filter = len(not_unique_list)
        variant_dict_pre_filter = sum([len(variant_dict[chrom]) for chrom in variant_dict])

        not_unique_list = reevaluating_ungrouped_reads(not_unique_list)

        unchosen_list = evaluate_unchosen_list(unchosen_list, haploblock_dict, not_unique_list)

        # TODO: Noget går galt med haploblock varianter
        not_unique_list = remove_removed_read_variants(haploblock_dict)

        merge_list = stitch_blocks(haploblock_dict, 0, merge_list)

        if haploblock_dict_pre_filter == len(haploblock_dict) and ungrouped_list_pre_filter == len(not_unique_list) and variant_dict_pre_filter == sum([len(variant_dict[chrom]) for chrom in variant_dict]):

            break

    # while True:

    #     haploblock_dict_pre_filter = len(haploblock_dict)
    #     ungrouped_list_pre_filter = len(not_unique_list)

    #     merge_list = stitch_blocks(haploblock_dict, maxdiff, merge_list)

    #     not_unique_list = reevaluating_ungrouped_reads(not_unique_list)

    #     unchosen_list = evaluate_unchosen_list(unchosen_list, haploblock_dict, not_unique_list)

    #     not_unique_list = make_unphasable_consensus(haploblock_dict, not_unique_list, False)

    #     if haploblock_dict_pre_filter == len(haploblock_dict) and ungrouped_list_pre_filter == len(not_unique_list):

    #         break

    if verbosity > 2:

        print("", file=sys.stderr)

    if verbosity > 1:

        print(f'After read selection the block count is {len(haploblock_dict)} and {sum([len(variant_dict[chrom]) for chrom in variant_dict])} usable variants remain.\n', file=sys.stderr)

    if args.output_intermediate_states:

        make_test_output(f'{args.outdir}/test_outputs/{args.prefix}_secondary_filtered_4_of_5_post_read_selection.bam')

    return hap_dict, merge_list, not_unique_list, unchosen_list


def group_supplementary_reads(haploblock_dict):

    if verbosity > 3:

        print(f'', file=sys.stderr)

    if verbosity > 2:

        print("Grouping supplementary reads.", file=sys.stderr)

    if verbosity > 3:

        print(f'', file=sys.stderr)

    for read in supplementary_read_dict:
        for start_pos in supplementary_read_dict[read]:

            chrom = supplementary_read_dict[read][start_pos][0].reference_name

            supplementary_read_dict[read][start_pos][2] = {}

            for pos in supplementary_read_dict[read][start_pos][1]:

                if pos in variant_dict[chrom]:

                    supplementary_read_dict[read][start_pos][2][pos] = supplementary_read_dict[read][start_pos][1][pos]

            read_variant_dict = supplementary_read_dict[read][start_pos][2]

            if len(read_variant_dict) == 0:

                continue

            fit_list = make_fit_list(read_variant_dict, supplementary_read_dict[read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, maxdiff, 1)

            if verbosity > 4:

                print(read, start_pos, len(read_variant_dict), fit_list, file=sys.stderr)

            if len(fit_list) == 1:

                hap = fit_list[0][0]

                supplementary_read_dict[read][start_pos][3] = hap

def create_haplotypes(haploblock_dict, not_unique_list, merge_list): # Calls reintroduce_read_to_variant_dict, remove_removed_read_variants, reevaluate_remaining_ungrouped_reads, make_unphasable_consensus

    print_subprocess_title("Merging blocks to create haplotypes.")

    for block in list(haploblock_dict.keys()):

        if len(haploblock_dict[block][1]) < minimum_coverage and sum([block in merge_vec for merge_vec in merge_list]) == 0:

            for read_type, read, start_pos in haploblock_dict[block][1]:

                reintroduce_read_to_variant_dict(read_type, read, start_pos)
                not_unique_list.append((read_type, read, start_pos))

            if verbosity > 2:

                print(f"{block} removed due to low coverage.", file=sys.stderr)

            del haploblock_dict[block]

    haplotype_dict = {}

    for chrom in variant_dict:
        for pos in variant_dict[chrom]:

            if sum([pos in haploblock_dict[hap][0] and hap not in haplotype_dict for hap in haploblock_dict]) < 1:

                continue

            haps = []
            hap_names = []

            for hap in haploblock_dict:

                if pos in haploblock_dict[hap][0]:

                    haps.append(hap)

                    if hap in haplotype_dict:

                        hap_names.append(haplotype_dict[hap])

            if verbosity > 4:
                
                print(pos, haps, hap_names, file=sys.stderr)

            count = 1

            for hap in haps:

                if hap not in haplotype_dict:

                    found = False

                    while not found:

                        if count not in hap_names:

                            haplotype_dict[hap] = count
                            hap_names.append(count)

                            found = True

                        else:

                            count += 1

            if verbosity > 4:
                
                print(pos, haplotype_dict, file=sys.stderr)

    if verbosity > 4:
        
        print(haplotype_dict, file=sys.stderr)

    count = 0

    seen_blocks = []

    for block1, block2, merge_status in merge_list:

        if block1 in seen_blocks or block2 in seen_blocks: # FIXME: test

            print("two merges contain the same block...") # FIXME: It need to be able to handle this

        block1_name = haplotype_dict[block1]
        block2_name = haplotype_dict[block2]

        count += 1

        for key, value in haplotype_dict.items():

            if value == block1_name and value == block2_name: # FIXME: probably unnecessary

                print("two blocks have the same name...")

            if value == block1_name:

                seen_blocks.append(key)

                if merge_status == "within":

                    haplotype_dict[key] = f'Hap{count}'

                else:

                    haplotype_dict[key] = f'Hap{count}_left'

            if value == block2_name:

                seen_blocks.append(key)

                if merge_status == "within":

                    haplotype_dict[key] = f'Hap{count}_within'

                else:

                    haplotype_dict[key] = f'Hap{count}_right'

    if verbosity > 4:
        
        print(haplotype_dict, file=sys.stderr)

    for block in haplotype_dict:

        if block not in seen_blocks:

            count += 1

            block_name = haplotype_dict[block]

            for key, value in haplotype_dict.items():

                if value == block_name:

                    seen_blocks.append(key)

                    haplotype_dict[key] = f'Hap{count}'

    if verbosity > 4:
        
        print(haplotype_dict, file=sys.stderr)

    for block in haplotype_dict:

        hap_name = haplotype_dict[block]

        if hap_name not in haploblock_dict:

            haploblock_dict[hap_name] = haploblock_dict[block]

        else:

            haploblock_dict[hap_name][0].update(haploblock_dict[block][0].copy())
            haploblock_dict[hap_name][1] += haploblock_dict[block][1]

        for read_type, read, start_pos in haploblock_dict[block][1]:

            read_dict[read_type][read][start_pos][3] = hap_name

        del haploblock_dict[block]

    while True:

        variant_dict_pre_filter = len(haploblock_dict)
        ungrouped_list_pre_filter = len(not_unique_list)

        not_unique_list = remove_removed_read_variants(haploblock_dict)

        not_unique_list, merge_list = reevaluate_remaining_ungrouped_reads(haploblock_dict, not_unique_list, merge_list)

        if variant_dict_pre_filter == len(haploblock_dict) and ungrouped_list_pre_filter == len(not_unique_list):

            break

    not_unique_list = make_unphasable_consensus(haploblock_dict, not_unique_list, False)

    not_unique_list = remove_removed_read_variants(haploblock_dict)

    group_supplementary_reads(haploblock_dict)

    if verbosity > 2:

        print(f'', file=sys.stderr)

    if verbosity > 1:

        print(f'After creation of haplotypes the block count is {len(haploblock_dict)}.\n', file=sys.stderr)

    if verbosity > 4:

        print(len(not_unique_list), sum([read_dict[read_type][read][start_pos][3] == "ungrouped" and len(read_dict[read_type][read][start_pos][2]) != 0 for read_type in read_dict for read in read_dict[read_type] for start_pos in read_dict[read_type][read]]), file=sys.stderr)
        test_haploblock_dict(haploblock_dict)
        test_not_unique_list(not_unique_list, haploblock_dict)
        test_variant_dict()

    if args.output_intermediate_states:

        make_test_output(f'{args.outdir}/test_outputs/{args.prefix}_secondary_filtered_5_of_5_post_haplotype_creation.bam')

    return haplotype_dict, not_unique_list


def add_base_to_haplo_variant_dict(pos, alt_base, hap_name, haplo_variant_dict):

    if hap_name not in haplo_variant_dict[pos]:

        haplo_variant_dict[pos][hap_name] = {}

    if alt_base not in haplo_variant_dict[pos][hap_name]:

        haplo_variant_dict[pos][hap_name][alt_base] = 0

    haplo_variant_dict[pos][hap_name][alt_base] += 1

def create_haplo_variant_dict(): # Calls get_reference_sequence, make_fit_list and add_base_to_haplo_variant_dict

    for chrom, start, end, _ in create_range_list():

        hap_unphasale_range = {}

        for hap in haploblock_dict:

            if "left" in hap:

                hap_unphasale_range[hap] = [start, 0]

            elif "right" in hap:

                hap_unphasale_range[hap] = [float("inf"), end]

            elif "within" in hap:

                hap_unphasale_range[hap] = [float("inf"), 0]

            else:

                hap_unphasale_range[hap] = [start, end]

                continue

            for pos in haploblock_dict[hap][0]:

                if len(haploblock_dict[hap][0][pos]) == 1:

                    if pos < hap_unphasale_range[hap][0]:

                        hap_unphasale_range[hap][0] = pos

                    if hap_unphasale_range[hap][1] < pos:

                        hap_unphasale_range[hap][1] = pos

        iter = samfile.pileup(chrom, start, end, stepper="nofilter")

        haplo_variant_dict = {}
        haplo_asm = {}
        secondary_dict = {"s": {}, "p": {}}

        for pileupcolumn in iter:

            pos = pileupcolumn.pos

            if pos < start or end < pos:

                continue

            pileupcolumn.set_min_base_quality(0)

            haplo_variant_dict[pos] = {}

            for pileupread in pileupcolumn.pileups:

                if pos < pileupread.alignment.reference_start or pileupread.alignment.reference_end < pos:

                    continue

                # Makes sure that the read has been used in the grouping

                if pileupread.alignment.is_secondary:

                    read_type = "s"
                    start_pos = pileupread.alignment.reference_start

                else:

                    read_type = "p"
                    start_pos = "-"

                read_name = pileupread.alignment.query_name

                if read_type == "s" and read_name not in read_dict["p"]:

                    continue

                if read_name not in read_dict[read_type]:

                    continue

                elif start_pos not in read_dict[read_type][read_name]:

                    continue

                # Makes sure that the read has been grouped with a haplotype

                hap_name = read_dict[read_type][pileupread.alignment.query_name][start_pos][3]

                if hap_name in ["removed"]:

                    continue

                if pileupread.alignment.is_secondary:

                    # Gets the reference sequence for the secondary read

                    if read_name not in secondary_dict[read_type]:

                        secondary_dict[read_type][read_name] = {}

                    if start_pos not in secondary_dict[read_type][read_name]:

                        secondary_dict[read_type][read_name][start_pos] = get_reference_sequence(read_type, read_name, start_pos, read_dict, True)

                    if pos-pileupread.alignment.reference_start >= len(secondary_dict[read_type][read_name][start_pos]):

                        continue

                    alt_base = secondary_dict[read_type][read_name][start_pos][pos-start_pos]

                else:

                    if pileupread.is_del:

                        alt_base = "D"

                    elif pileupread.indel > 0:

                        # Gets the reference sequence for the primary read

                        if read_name not in secondary_dict[read_type]:

                            secondary_dict[read_type][read_name] = {}

                        if start_pos not in secondary_dict[read_type][read_name]:

                            secondary_dict[read_type][read_name][start_pos] = get_reference_sequence(read_type, read_name, start_pos, read_dict, True)

                        if pos-pileupread.alignment.reference_start >= len(secondary_dict[read_type][read_name][start_pos]):

                            continue

                        alt_base = secondary_dict[read_type][read_name][start_pos][pos-pileupread.alignment.reference_start]

                    else:

                        alt_base = pileupread.alignment.query_sequence[pileupread.query_position]

                if len(haploblock_dict) == 0:

                    add_base_to_haplo_variant_dict(pos, alt_base, hap_name, haplo_variant_dict)

                elif hap_name == "unphasable":

                    for hap in haploblock_dict:

                        if hap_unphasale_range[hap][0] <= pos and pos <= hap_unphasale_range[hap][1]:

                            add_base_to_haplo_variant_dict(pos, alt_base, hap, haplo_variant_dict)

                elif hap_name in ["ungrouped", "unchosen"]:

                    read_variant_dict = read_dict[read_type][read_name][start_pos][2]

                    fit_list = make_fit_list(read_variant_dict, read_dict[read_type][read_name][start_pos][0].reference_name, read_dict[read_type][read_name][start_pos][4], haploblock_dict, 0, 1)

                    for i in range(len(fit_list)):

                        add_base_to_haplo_variant_dict(pos, alt_base, fit_list[i][0], haplo_variant_dict)

                else:

                    add_base_to_haplo_variant_dict(pos, alt_base, hap_name, haplo_variant_dict)

            for hap in haplo_variant_dict[pos]:

                if hap not in haplo_asm:

                    haplo_asm[hap] = {}
                    haplo_asm[hap]["asm"] = ""
                    haplo_asm[hap]["count"] = []

                base = max(haplo_variant_dict[pos][hap], key = haplo_variant_dict[pos][hap].get)

                if base == "D":

                    continue

                if len(base) > 1:

                    for _ in range(len(base)):

                        haplo_asm[hap]["count"].append(haplo_variant_dict[pos][hap])

                elif len(base) == 1:

                    haplo_asm[hap]["count"].append(haplo_variant_dict[pos][hap])

                else: # FIXME: test

                    print("If this prints something is wrong in creare_haplo_variant_dict.")

                haplo_asm[hap]["asm"] += base

    if args.output_raw_assemblies:

        asm_outfile = f'{args.outdir}/{args.prefix}_raw_asm_all.fasta'

        if verbosity > 2:

            print(f"Writing all raw haplotype assemblies to: {asm_outfile}", file=sys.stderr)

        outfile = open(asm_outfile, "w")

        for hap in haplo_asm:

            outfile.write(">" + hap + "\n" + haplo_asm[hap]["asm"] + "\n")

        outfile.close()

        sam_file = f'{args.outdir}/{args.prefix}_raw_asm_all.sam'
        sorted_sam_file = f'{args.outdir}/{args.prefix}_raw_asm_all_sorted.bam'

        subprocess.run(["minimap2", "-a", "-x", "asm5", "--cs", "-r2k", "-o", sam_file, args.reference, asm_outfile])
        subprocess.run(["samtools", "sort", "-o", sorted_sam_file, sam_file])
        subprocess.run(["samtools", "index", sorted_sam_file])

    return haplo_asm

def find_best_Levenshtein_distance(haplo_asm, left_hap, right_hap):

    max_overlap = min([len(haplo_asm[left_hap]["asm"]), len(haplo_asm[right_hap]["asm"]), args.maximal_block_overlap])

    print(max_overlap)

    approximate_best_overlap = None
    approximate_best_distance = max_overlap

    for overlap in range(min_overlap, max_overlap, overlap_jump):

        left_flank = haplo_asm[left_hap]["asm"][len(haplo_asm[left_hap]["asm"])-overlap:]
        right_flank = haplo_asm[right_hap]["asm"][:overlap]

        distance = Levenshtein.distance(left_flank, right_flank)

        if verbosity > 4:

            print(overlap, distance, Levenshtein.ratio(left_flank, right_flank), file=sys.stderr)

        if distance < approximate_best_distance and Levenshtein.ratio(left_flank, right_flank) > min_overlap_identity: # FIXME: maybe this should just be ratio-based

            approximate_best_overlap = overlap
            approximate_best_distance = distance

    best_overlap = None
    best_distance = max_overlap

    for overlap in range(approximate_best_overlap-approximate_best_distance, approximate_best_overlap+approximate_best_distance+1):

        left_flank = haplo_asm[left_hap]["asm"][len(haplo_asm[left_hap]["asm"])-overlap:]
        right_flank = haplo_asm[right_hap]["asm"][:overlap]

        distance = Levenshtein.distance(left_flank, right_flank)

        if verbosity > 4:

            print(overlap, Levenshtein.distance(left_flank, right_flank), Levenshtein.ratio(left_flank, right_flank), file=sys.stderr)

        if distance < best_distance:

            best_overlap = overlap
            best_distance = distance

    return best_overlap, best_distance

def make_haplotype_consensus_sequence(haplo_dict): # Calls create_haplo_variant_dict

    print_subprocess_title("Making haplotypes consensus sequences.")

    haplo_asm = create_haplo_variant_dict()

    name_list = []
    seq_list = []
    seen_haplotypes = []

    for block1, block2, merge_status in merge_list:

        if merge_status == "within":

            continue

        left_hap = haplo_dict[block1]
        right_hap = haplo_dict[block2]

        seen_haplotypes.append(left_hap)
        seen_haplotypes.append(right_hap)

        if verbosity > 2:

            print(f'Evaluating optimal overlap between {left_hap} and {right_hap}.', file=sys.stderr)

        i = 0

        left_found = False
        right_found = False

        left_count = haplo_asm[left_hap]["count"]
        right_count = haplo_asm[right_hap]["count"]

        while i < 40000:

            if sum(left_count[-(i+1)].values()) >= minimum_count:

                if not left_found:

                    left_i = i

                    left_found = True

            else:

                left_found = False

            if sum(right_count[i].values()) >= minimum_count:

                if not right_found:

                    right_i = i

                    right_found = True

            else:

                right_found = False

            i += 1

        haplo_asm[left_hap]["asm"] = haplo_asm[left_hap]["asm"][:len(haplo_asm[left_hap]["asm"])-left_i]
        haplo_asm[left_hap]["count"] = haplo_asm[left_hap]["count"][:len(haplo_asm[left_hap]["count"])-left_i]

        haplo_asm[right_hap]["asm"] = haplo_asm[right_hap]["asm"][right_i:]
        haplo_asm[right_hap]["count"] = haplo_asm[right_hap]["count"][right_i:]

        best_overlap, best_distance = find_best_Levenshtein_distance(haplo_asm, left_hap, right_hap)

        left_flank = haplo_asm[left_hap]["asm"][len(haplo_asm[left_hap]["asm"])-best_overlap:]
        left_count = haplo_asm[left_hap]["count"][len(haplo_asm[left_hap]["count"])-best_overlap:]
        right_flank = haplo_asm[right_hap]["asm"][:best_overlap]
        right_count = haplo_asm[right_hap]["count"][:best_overlap]

        flank = 10
        best_asm = (None, None, 1, None)

        for overlap in range(best_overlap-best_distance, best_overlap+best_distance+1):

            left_flank = haplo_asm[left_hap]["asm"][len(haplo_asm[left_hap]["asm"])-overlap:]
            left_count = haplo_asm[left_hap]["count"][len(haplo_asm[left_hap]["count"])-overlap:]
            right_flank = haplo_asm[right_hap]["asm"][:overlap]
            right_count = haplo_asm[right_hap]["count"][:overlap]

            right_overlap = overlap

            total_overlap = 0

            i = 0
            j = 0
            running_overlap = 0
            merged_asm = ""
            mismatch_count = 0
            mismatch_list = []

            while i < overlap:

                if left_flank[i] == right_flank[j]:

                    merged_asm += left_flank[i]

                    i += 1
                    j += 1
                    running_overlap += 1

                    total_overlap += 1

                else:

                    mismatch_list.append(["match", i, j, running_overlap])
                    running_overlap = 0

                    found = False
                    indel = 1

                    if i+indel+flank > overlap:

                        for x in range(overlap-i):

                            if left_flank[i+x] != right_flank[j+x]:

                                mismatch_count += 1

                        i += x + 1
                        j += x + 1

                        total_overlap += x

                        found = True

                        break

                    while indel < 10000:
                    # while True:

                        if i+indel+flank == overlap:

                            for x in range(indel+flank):

                                if left_flank[i+x] != right_flank[j+x]:

                                    mismatch_count += 1

                            i += indel+flank
                            j += indel+flank
                            total_overlap += indel+flank

                            found = True

                            break

                        mismatch_indel = 0

                        while mismatch_indel < 10:

                            if left_flank[i+mismatch_indel:i+mismatch_indel+flank] == right_flank[j+mismatch_indel+indel:j+mismatch_indel+indel+flank]:

                                merged_asm += right_flank[j:j+indel]

                                i += mismatch_indel
                                j += indel + mismatch_indel
                                mismatch_count += indel + mismatch_indel
                                total_overlap += indel + mismatch_indel

                                found = True

                                break

                            elif left_flank[i+mismatch_indel+indel:i+mismatch_indel+indel+flank] == right_flank[j+mismatch_indel:j+mismatch_indel+flank]:

                                merged_asm += left_flank[i:i+indel]

                                i += indel + mismatch_indel
                                j += mismatch_indel
                                mismatch_count += indel + mismatch_indel
                                total_overlap += indel + mismatch_indel

                                found = True

                                break

                            elif left_flank[i+indel:i+indel+flank] == right_flank[j+indel:j+indel+flank]:

                                if left_count[i][max(left_count[i], key = left_count[i].get)] > right_count[j][max(right_count[j], key = right_count[j].get)]:

                                    merged_asm += left_flank[i:i+indel]

                                else:

                                    merged_asm += right_flank[j:j+indel]

                                for x in range(indel):

                                    if left_flank[i+x] != right_flank[j+x]:

                                        mismatch_count += 1

                                i += indel
                                j += indel
                                total_overlap += indel

                                found = True

                                break

                            mismatch_indel += 1

                        if found:

                            if right_overlap-j != overlap-i:

                                right_overlap = j+overlap-i
                                right_flank = haplo_asm[right_hap]["asm"][:right_overlap]
                                right_count = haplo_asm[right_hap]["count"][:right_overlap]

                            break

                        indel += 1

                    if not found:

                        break

            if i == overlap:

                if best_asm[2] > mismatch_count/total_overlap:

                    best_asm = (overlap, right_overlap, mismatch_count/total_overlap, merged_asm)

        left_flank = haplo_asm[left_hap]["asm"][len(haplo_asm[left_hap]["asm"])-best_asm[0]:]
        right_flank = haplo_asm[right_hap]["asm"][:best_asm[1]]

        if verbosity > 2:

            print(f'The optimal overlap between {left_hap} and {right_hap} is {best_asm[0]} with a levenshtein distance of {Levenshtein.distance(left_flank, right_flank)}.', file=sys.stderr)

        if verbosity > 4:
            
            print(best_asm[0], best_asm[1], Levenshtein.distance(left_flank, right_flank), Levenshtein.ratio(left_flank, right_flank), file=sys.stderr)
            print(Levenshtein.editops(left_flank, right_flank), file=sys.stderr)

        name_list.append(left_hap.split("_")[0])
        seq_list.append(haplo_asm[left_hap]["asm"][:len(haplo_asm[left_hap]["asm"])-best_asm[0]] + best_asm[3] + haplo_asm[right_hap]["asm"][best_asm[1]:])

    for hap in haplo_asm:

        if hap not in seen_haplotypes and hap not in ["ungrouped", "unchosen", "removed"]:

            seen_haplotypes.append(hap)

            name_list.append(hap)
            seq_list.append(haplo_asm[hap]["asm"])

    return name_list, seq_list

def write_bam_file(bam_file):

    if verbosity > 2:

        print(f'\nWriting filtered alignment outfile.\n', file=sys.stderr)

    outfile = pysam.AlignmentFile(bam_file, "wb", template=samfile)

    for read_type in read_dict:
        for read in read_dict[read_type]:
            for start_pos in read_dict[read_type][read]:

                if read_type == "s" and read in read_dict["p"] and read_dict["p"][read]["-"][3] == "removed":

                    read_dict[read_type][read][start_pos][0].flag &= ~pysam.FSECONDARY
                    read_dict[read_type][read][start_pos][0].mapping_quality = 60

                if args.merge_overlap and ("left" in read_dict[read_type][read][start_pos][3] or "right" in read_dict[read_type][read][start_pos][3]):

                    hap_name = read_dict[read_type][read][start_pos][3].split("_")[0]

                else:

                    hap_name = read_dict[read_type][read][start_pos][3]

                if not show_removed and hap_name == "removed":

                    continue

                read_alignment = read_dict[read_type][read][start_pos][0]

                read_alignment.set_tag("HP", hap_name, replace=True)

                outfile.write(read_alignment)

    for read in supplementary_read_dict:
        for start_pos in supplementary_read_dict[read]:

            if read not in supplementary_read_dict: # FIXME: Dunno why this is necessary

                continue

            hap_name = supplementary_read_dict[read][start_pos][3]

            read_alignment = supplementary_read_dict[read][start_pos][0]

            read_alignment.set_tag("HP", hap_name, replace=True)

            outfile.write(read_alignment)

    for read in unmatched_secondary_dict:
        for start_pos in unmatched_secondary_dict[read]:

            if read not in unmatched_secondary_dict: # FIXME: Dunno why this is necessary

                continue

            hap_name = unmatched_secondary_dict[read][start_pos][3]

            read_alignment = unmatched_secondary_dict[read][start_pos][0]

            read_alignment.set_tag("HP", hap_name, replace=True)

            outfile.write(read_alignment)

    outfile.close()

def produce_outfiles(name_list, seq_list):

    print_subprocess_title("Producing outfiles.")

    alignment_outfile = f'{args.outdir}/{args.prefix}_NanoCAH_filtered_alignment.bam'
    sorted_alignment_outfile = f'{args.outdir}/{args.prefix}_NanoCAH_filtered_alignment_sorted.bam'

    if verbosity > 2:

        print(f"Writing filtered alignment file to: {alignment_outfile}", file=sys.stderr)

    write_bam_file(alignment_outfile)

    subprocess.run(["samtools", "sort", "-o", sorted_alignment_outfile, alignment_outfile])
    subprocess.run(["samtools", "index", sorted_alignment_outfile])

    asm_outfile = f'{args.outdir}/{args.prefix}_asm_all.fasta'

    if verbosity > 2:

        print(f"Writing all haplotype assemblies to: {asm_outfile}", file=sys.stderr)

    outfile = open(asm_outfile, "w")

    for i in range(len(seq_list)):

        outfile.write(">" + name_list[i] + "\n" + seq_list[i] + "\n")

    outfile.close()

    sam_file = f'{args.outdir}/{args.prefix}_asm_all.sam'
    sorted_sam_file = f'{args.outdir}/{args.prefix}_asm_all_sorted.bam'

    subprocess.run(["minimap2", "-a", "-x", "asm5", "--cs", "-r2k", "-o", sam_file, args.reference, asm_outfile])
    subprocess.run(["samtools", "sort", "-o", sorted_sam_file, sam_file])
    subprocess.run(["samtools", "index", sorted_sam_file])

    for i in range(len(name_list)):

        hap = name_list[i]

        asm_outfile = f'{args.outdir}/{args.prefix}_asm_{hap}.fasta'

        if verbosity > 2:

            print(f"Writing {hap} assembly to: {asm_outfile}", file=sys.stderr)

        outfile = open(asm_outfile, "w")

        outfile.write(">" + name_list[i] + "\n" + seq_list[i] + "\n")

        outfile.close()

        sam_file = f'{args.outdir}/{args.prefix}_asm_{hap}.sam'
        sorted_sam_file = f'{args.outdir}/{args.prefix}_asm_{hap}_sorted.bam'

        subprocess.run(["minimap2", "-a", "-x", "asm5", "--cs", "-r2k", "-o", sam_file, args.reference, asm_outfile])
        subprocess.run(["samtools", "sort", "-o", sorted_sam_file, sam_file])
        subprocess.run(["samtools", "index", sorted_sam_file])

    if "Hap1" in name_list and "Hap2" in name_list:

        subprocess.run(["svim-asm", "diploid", f'{args.outdir}/{args.prefix}_SVIM', f'{args.outdir}/{args.prefix}_asm_Hap1_sorted.bam', f'{args.outdir}/{args.prefix}_asm_Hap2_sorted.bam', args.reference])

        subprocess.run(["mv", f'{args.outdir}/{args.prefix}_SVIM/variants.vcf', f'{args.outdir}/{args.prefix}_SV.vcf'])

        subprocess.run(["rm", "-r", f'{args.outdir}/{args.prefix}_SVIM'])
    
    else:

        if verbosity > 2:

            print(f"It was not possible to run structural variant calling with SVIM-ASM, because there is not at least two haplotypes.", file=sys.stderr)


def make_test_output(alignment_outfile):

    if verbosity > 2:

        print(f"Writing filtered alignment file to: {alignment_outfile}", file=sys.stderr)

    write_bam_file(f'{alignment_outfile}.temp')

    subprocess.run(["samtools", "sort", "-o", alignment_outfile, f'{alignment_outfile}.temp'])
    subprocess.run(["samtools", "index", alignment_outfile])
    subprocess.run(["rm", f'{alignment_outfile}.temp'])

def test_read_dict():

    for read_type in read_dict:
        for read in read_dict[read_type]:
            for start_pos in read_dict[read_type][read]:

                if read_dict[read_type][read][start_pos][3] != "unphasable":

                    continue

                print("\n", file=sys.stderr)
                print(read_type, read, start_pos, len(read_dict[read_type][read][start_pos][1]), len(read_dict[read_type][read][start_pos][2]), read_dict[read_type][read][start_pos][3], file=sys.stderr)
                print(make_fit_list(read_dict[read_type][read][start_pos][2], read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], read_dict[read_type][read][start_pos][0].reference_name, haploblock_dict, 0, 1), file=sys.stderr)
                print(make_fit_list(read_dict[read_type][read][start_pos][2], read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], read_dict[read_type][read][start_pos][0].reference_name, haploblock_dict, maxdiff, 1), file=sys.stderr)

def test_haploblock_dict(haploblock_dict):

    for hap in haploblock_dict:

        hap_variant_dict = haploblock_dict[hap][0]
        chrom = haploblock_dict[hap][2]

        fit_list = make_fit_list(hap_variant_dict, chrom, haploblock_dict[hap][3], haploblock_dict, maxdiff, 1)

        print(hap, len(hap_variant_dict), len(haploblock_dict[hap][1]), len(set([f'{read_type}_{read}_{start_pos}' for read_type, read, start_pos in haploblock_dict[hap][1]])), fit_list, file=sys.stderr)

        for block, _, _, mismatch_list in fit_list:

            if block == hap:

                continue

            else:

                print(hap, block, [(pos, base1, sum([pos in read_dict[read_type][read][start_pos][2] and read_dict[read_type][read][start_pos][2][pos] == base1 for read_type, read, start_pos in haploblock_dict[block][1]]), base2, sum([pos in read_dict[read_type][read][start_pos][2] and read_dict[read_type][read][start_pos][2][pos] == base2 for read_type, read, start_pos in haploblock_dict[hap][1]]), [(base, len(variant_dict[chrom][pos][base])) for base in variant_dict[chrom][pos]]) for pos, base1, base2 in mismatch_list], file=sys.stderr)

    new_dict = {}

    for read_type in read_dict:
        for read in read_dict[read_type]:
            for start_pos in read_dict[read_type][read]:

                if read_dict[read_type][read][start_pos][3] not in new_dict:

                    new_dict[read_dict[read_type][read][start_pos][3]] = 0

                new_dict[read_dict[read_type][read][start_pos][3]] += 1

    print(new_dict, "\n", file=sys.stderr)

def test_not_unique_list(not_unique_list, haploblock_dict):

    for read_type, read, start_pos in not_unique_list:

        print(read_type, read, start_pos, len(read_dict[read_type][read][start_pos][1]), len(read_dict[read_type][read][start_pos][2]), read_dict[read_type][read][start_pos][3], file=sys.stderr)
        print(make_fit_list(read_dict[read_type][read][start_pos][2], read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, maxdiff, 1), file=sys.stderr)

        if len(make_fit_list(read_dict[read_type][read][start_pos][2], read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, maxdiff, 1)) == 0:

            print(make_fit_list(read_dict[read_type][read][start_pos][2], read_dict[read_type][read][start_pos][0].reference_name, read_dict[read_type][read][start_pos][4], haploblock_dict, maxdiff, 0), file=sys.stderr)
    
    print("\n", file=sys.stderr)

def test_variant_dict():

    for chrom in variant_dict:
        for pos in variant_dict[chrom]:

            new_dict = {}

            for base in variant_dict[chrom][pos]:

                new_dict[base] = {}

                for read_type, read, start_pos in variant_dict[chrom][pos][base]:

                    if read_dict[read_type][read][start_pos][3] not in new_dict[base]:

                        new_dict[base][read_dict[read_type][read][start_pos][3]] = 0

                    new_dict[base][read_dict[read_type][read][start_pos][3]] += 1

            print(chrom, pos, [(base, len(variant_dict[chrom][pos][base])) for base in variant_dict[chrom][pos]], new_dict, sum([len(variant_dict[chrom][pos][base]) if base != "deleted" else 0 for base in variant_dict[chrom][pos]]), depth_dict[chrom][pos], sum([len(variant_dict[chrom][pos][base]) if base != "deleted" else 0 for base in variant_dict[chrom][pos]])/depth_dict[chrom][pos], file=sys.stderr)
    
    print("\n", file=sys.stderr)


####################################################################################################
# CODE                                                                                             #
####################################################################################################

# TODO: Depth_dict not used for anything yet
# TODO: Note large deletions.
# TODO: Make input variable print

variant_dict, read_dict, depth_dict, supplementary_read_dict, unmatched_secondary_dict = create_dicts(samfile)

remove_bad_variants(variant_dict)

haploblock_dict, not_unique_list, hap_count = add_reads_to_haplotype_dict()

haploblock_dict, not_unique_list = perform_no_variance_merging(haploblock_dict, not_unique_list)

haploblock_dict, not_unique_list = perform_low_variance_merging(haploblock_dict, not_unique_list)
# test_haploblock_dict(haploblock_dict)
# test_variant_dict()
hap_dict, merge_list, not_unique_list, unchosen_list = reevaluate_secondary_reads(not_unique_list, haploblock_dict)

# exit()

# TODO: Make haploblock beds

haplo_dict, not_unique_list = create_haplotypes(haploblock_dict, not_unique_list, merge_list)

name_list, seq_list = make_haplotype_consensus_sequence(haplo_dict)

produce_outfiles(name_list, seq_list)
