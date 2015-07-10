#!/usr/bin/env python3

import re
from Bio import SeqIO

# This script provides several functions for the detection of sensitive
# crispr motifs.

def fasta_processor(filename):

    """
    Given a fasta generates a unique string with it to be processed by the
    next functions
    """

    lines = [line.rstrip('\n') for line in open(filename)]

    # remove > line
    del lines[0]

    # Join sequence
    sequence = "".join(lines)

    return sequence


def insertion_processor(inputfilename, ins_type, outputfilename):
    """
    Processes the file with insertions and types and generates a plain text file with the positions of insertions
    for the insertion type selected
    """

    filehandle = open(inputfilename, "r")
    fo = open(outputfilename, "w")

    count = 0
    for line in filehandle:
        words = line.split()
        try:
            if words[0] == ins_type:
                fo.write(words[1]+'\n')
                count += 1
        except:
            pass
    print(count)

    filehandle.close()
    fo.close()


def genename_associator(genbank_file, genenames_file):
    """
    Given a genbank file and a list of genenames, it generates a file with only the MPN nomenclature names.
    """

    # Generate a list with the genenames
    fi = open(genenames_file, "r")

    classic_genenames = []

    for line in fi:
        line = line.split("\t")
        genename = line[0]
        if not genename.startswith("MPN"):
            classic_genenames.append(genename)

    # Generate the dictionary of results:
    genenames_d = {}

    for record in SeqIO.parse(genbank_file, "genbank"):
        for f in record.features:
            if f.type == "CDS" and "gene" in f.qualifiers:
                gene = f.qualifiers["gene"][0]
                if gene in classic_genenames:
                    try:
                        genenames_d[gene] = f.qualifiers["locus_tag"][0]
                    except:
                        print("error parsing the gb file")

    # Generate the new file:
    fo = open("../ref_data/mpn_genenames_coord.txt", "w")
    fu = open(genenames_file, "r")

    for line in fu:
        line = line.split("\t")
        genename = line[0]
        start = line[1]
        end = line[2]

        if genename.startswith("MPN"):
            fo.write(genename+"\t"+"\t".join(line))
        else:
            try:
                fo.write(genenames_d[genename]+"\t"+genename+"\t"+start+"\t"+end)
            except:
                fo.write(genename+"\t"+genename+"\t"+start+"\t"+end)

    fi.close()
    fo.close()
    fu.close()


def coding_coordinates(inputfilename, outputfilename, negative = None):
    """
    Processes the file with genes and coordinates
    """

    filehandle = open(inputfilename, "r")
    fo = open(outputfilename, "w")

    count = 0

    for line in filehandle:
        words = line.split()
        try:
            if words[0].startswith(">"):
                count +=1
                genename = words[1][6:-1]
                coordinates = words[-1][10:-1]
                coordinates = coordinates.replace("location=", "")
                try:
                    coordinates = coordinates.replace("complement(", "")
                    coordinates = coordinates.replace(")","")
                except:
                    pass

                coordinates = coordinates.replace("..", "\t")
                coordinates = coordinates.split("\t")
                start = coordinates[0]
                end   = coordinates[1]

                fo.write(genename+'\t'+start+'\t'+end+'\n')
        except:
            pass


    filehandle.close()
    fo.close()
    print(count)

def reverse_complement(seq):

    """
    Given a sequence it returns its reverse complementary
    """

    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

    for k,v in alt_map.items():
        seq = seq.replace(k,v)

    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)

    for k,v in alt_map.items():
        bases = bases.replace(v,k)

    return bases


def motif_finder(sequence, motif):

    """
    Given a sequence and a motif, it returns a dictionary= {position: sequence}
    """

    offset = 0
    listindex = []

    i = sequence.find(motif, offset)
    while i >= 0:
        listindex.append(i)
        i = sequence.find(motif, i + 1)

    return listindex


def subseq_extractor(sequence, indices, length_motif, negative_strand = None):

    """
    Given a list of indexes it finds them in the original sequence and returns a dictionary
    with the start of the motif (defined with the length of the motif defined and the subsequence
    found

    Structure of seq:

    5'---------------------PAM3'
                 ...7654321

    1 = start
    20 = end
    """

    list_start_seq = []

    for index in indices:
        new_result = []
        start = index-2
        end =  index-length_motif

        if end < 0:
            subsequence = sequence[end:] + sequence[0:start+1]
        else:
            subsequence = sequence[end:start+1]

        if negative_strand:
            new_start = len(sequence) - start
            new_result.append(new_start)
        else:
            new_result.append(start + 1)

        new_result.append(subsequence)

        list_start_seq += [new_result,]

    list_start_seq.sort(key = lambda x: x[0])

    return list_start_seq

def remove_redundant(genome, list_start_seq):
    """
    For all the 20 mer given returns a list of list having those that appears only one time in the
    genome
    """

    redundant_list = []

    m = 0

    circ_genome = genome[-20:] + genome[0:21]
    rev_comp_genome = reverse_complement(genome)

    while m < len(list_start_seq):
        current_seq = list_start_seq[m][1]

        # Direct genome
        counts = motif_finder(genome, current_seq)

        # Circular
        circ_counts = motif_finder(circ_genome, current_seq)

        # Reversible
        rev_comp_counts = motif_finder(rev_comp_genome, current_seq)


        if len(counts) == 1 and len(circ_counts) == 0 and len(rev_comp_counts) == 0:
            unique_list.append(list_start_seq[m])
        elif len(counts) == 0 and len(circ_counts) == 1 and len(rev_comp_counts) == 0:
            unique_list.append(list_start_seq[m])
        elif len(counts) == 0 and len(circ_counts) == 0 and len(rev_comp_counts) == 1:
            unique_list.append(list_start_seq[m])
        else:
            pass

        m += 1

    return unique_list


def take_redundant(genome, list_start_seq):
    """
    For all the 20 mer given returns a list of list having only those that appears more than one time in the
    genome
    """

    unique_list = []

    m = 0

    circ_genome = genome[-24:] + genome[0:26]
    rev_comp_genome = reverse_complement(genome)

    num_repeats = 0

    while m < len(list_start_seq):
        Acurrent_seq = list_start_seq[m][1]+"AGG"
        Ccurrent_seq = list_start_seq[m][1]+"CGG"
        Gcurrent_seq = list_start_seq[m][1]+"GGG"
        Tcurrent_seq = list_start_seq[m][1]+"TGG"

        # Direct genome
        counts = len(motif_finder(genome, Acurrent_seq)) + len(motif_finder(genome, Ccurrent_seq)) + len(motif_finder(genome, Gcurrent_seq)) + len(motif_finder(genome, Tcurrent_seq))

        # Circular
        circ_counts = len(motif_finder(circ_genome, Acurrent_seq)) + len(motif_finder(circ_genome, Ccurrent_seq)) + len(motif_finder(circ_genome, Gcurrent_seq)) + len(motif_finder(circ_genome, Tcurrent_seq))

        # Reversible
        rev_counts = len(motif_finder(rev_comp_genome, Acurrent_seq)) + len(motif_finder(rev_comp_genome, Ccurrent_seq)) + len(motif_finder(rev_comp_genome, Gcurrent_seq)) + len(motif_finder(rev_comp_genome, Tcurrent_seq))


        counts += circ_counts + rev_counts

        if counts > 1:
            result = list_start_seq[m]
            result[1] = result[1]+"-"+str(counts)
            if result not in unique_list:
                unique_list.append(result)
        else:
            pass

        m += 1

    return unique_list



def interruption(gene_coordinates, results, motif_length = 20, negative_strand = None, index=0):

    """
    Given the file of gene coordinates and results it generates a list of list having all the positions and the possible
    intersections
    """
    # We start defining a list, having the starts and ends of the genes:
    genes = []

    filehandle = open(gene_coordinates, "r")

    for gene_result in filehandle:
        new_gene = []
        gene_result = gene_result.split()

        genename    = gene_result[index]
        start_gene  = int(gene_result[2])
        end_gene    = int(gene_result[3])

        new_gene = [genename, start_gene, end_gene]
        genes += [new_gene,]

    # We want to add a genename or nothing to each one of the results, in case a oligo overlaps with more than one result
    # we need to consider the result both; to do this we need to iterate all the list of genes for each result:

    output = []

    for pos_bef_PAM, subsequence in results:
        new_result = [pos_bef_PAM, subsequence]

        # Different data settings depending on the strand we work:
        if not negative_strand:
            start_motif = pos_bef_PAM - 19
            end_motif = int(pos_bef_PAM)
        else:
            start_motif = int(pos_bef_PAM)
            end_motif = pos_bef_PAM + 19

        # Iterate the gene results
        for genename, start_gene, end_gene in genes:

            # Find intersections
            if (start_motif >= start_gene and start_motif <= end_gene) or (end_motif >= start_gene and end_motif <= end_gene):
                new_result.append(genename)
            else:
                pass

        # Considering circularity. If the list of genes is ordered, the last coordinates arriving here
        # are the ones of the last gene.

        if start_motif < 0:
            new_start_motif = start_motif + 816394
            if (new_start_motif >= 815526 and new_start_motif <= 816338):
                new_result.append(genename)
            else:
                pass

        # If we arrive here without any intersection (len(result) = 2) , add the results without genename
        if len(new_result) == 2:
            genename = "-"
            new_result.append(genename)

        # Add the new result:
        output += [new_result,]

    return output


def insertion_detector_counter(ins_positions_41, ins_positions_7, unique_inter_results, negative_strand):
    """
    Given the list of insertion positions and the list of unique results it assigns the number of insertions
    InsTh41 and InsTh7 before and after the 17 position and the total
    """

    # Process the list of 41 and 7 insertion elements
    filehandle_41 = open(ins_positions_41, "r")
    filehandle_7 = open(ins_positions_7, "r")

    def generate_set(filehandle):
        a_set = set()
        for line in filehandle:
            line = line.split()
            a_set.add(int(line[0]))
        return a_set

    s_41 = generate_set(filehandle_41)
    s_7  = generate_set(filehandle_7)

    filehandle_41.close()
    filehandle_7.close()

    # Process the input list in order to transform all the genenames (in case a result having more than one match)
    # in only one string

    processed_list = []

    for result in unique_inter_results:
        new = []
        position = result[0]
        subsequence = result[1]
        gene = ", ".join(result[2:])
        new = [position, subsequence, gene]

        processed_list += [new,]

    # Start the counting process
    results = []

    for pos_bef_PAM, subsequence, genename in processed_list:
        new_result = []
        bef17_41, aft17_41, bef17_7, aft17_7, total41, total7, total = (0,)*7

        # Different data settings depending on the strand we work:
        if not negative_strand:
            start_motif = pos_bef_PAM - 19
            end_motif = int(pos_bef_PAM) + 2
        else:
            start_motif = int(pos_bef_PAM) - 2
            end_motif = pos_bef_PAM + 19

        # Count the number of insertions
        for position in range(start_motif , start_motif + 17):
            if position in s_41:
                bef17_41 += 1

            if position in s_7:
                bef17_7 += 1

        for position in range(end_motif-2 , end_motif + 1):
            if position in s_41:
                aft17_41 += 1

            if position in s_7:
                aft17_7 += 1

        total41 = bef17_41 + aft17_41
        total7  = bef17_7  + aft17_7
        total = total41 + total7

        new_result = [pos_bef_PAM, subsequence, genename, bef17_41, aft17_41, total41, bef17_7, aft17_7, total7, total]
        results += [new_result,]

    return results


def file_generator(list_of_lists, fileoutname):
    """
    It generates plain text files containing the list o lists content
    """

    fo = open(fileoutname, "w")

    fo.write("pos_next_PAM\tsubsequence\tgenename\t41_bef17\t41_aft17\ttotal41\t7_bef17\t7_aft17\ttotal7\ttotal\n")

    for result in list_of_lists:
        result = map(str, result)
        fo.write("\t".join(result)+"\n")

    fo.close()

def order_total_file(total_file):
    """
    Given the total file it orders it by the first column (position in the genome)
    """

    fi = open(total_file, "r")

    lines = []
    for line in fi:
        line = line.strip("\n")
        line = line.split("\t")
        line[0] = int(line[0])
        lines.append(line)

    lines.sort()

    file_generator(lines, "../results/total_file_processed.txt")

    fi.close()
