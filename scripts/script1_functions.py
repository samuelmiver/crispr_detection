#!/usr/bin/env python3

import re

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
    For all the 20 mer given returns a list of list having only those that appears only one time in the
    genome
    """

    unique_list = []

    m = 0

    while m < len(list_start_seq):
        current_seq = list_start_seq[m][1]

        # Direct genome
        counts = motif_finder(genome, current_seq)

        # Circular
        circ_genome = genome[-20:] + genome[0:21]
        circ_counts = motif_finder(circ_genome, current_seq)

        # Reversible
        rev_comp_genome = reverse_complement(genome)
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


def interruption(gene_coordinates, results, motif_length = 20):

    """
    Given the file of gene coordinates and results it generates a list of list having all the positions and the possible
    intersections
    """

    filehandle = open(gene_coordinates, "r")

    output = []

    for pos_bef_PAM, subsequence in results:
        new_result = []

        start_motif = pos_bef_PAM - 19
        end_motif = int(pos_bef_PAM)

        for line in filehandle:
            line = line.split()
            genename = line[0]
            start_gene = int(line[1])
            end_gene = int(line[2])

            # Find intersections
            if (start_motif >= start_gene and start_motif <= end_gene) or (end_motif >= start_gene and end_motif <= end_gene):
                new_result = [pos_bef_PAM, subsequence, genename]
                output += [new_result,]
            else:
                pass

        # Considering circularity. If the list of genes is ordered, the last coordinates arriving here
        # are the ones of the last gene.

        if start_motif < 0:
            new_start_motif = start_motif + 50
            if (new_start_motif >= start_gene and new_start_motif <= end_gene) or (end_motif >= start_gene and end_motif <= end_gene):
                new_result = [pos_bef_PAM, subsequence, genename]
                output += [new_result,]
            else:
                pass

        # If we arrive here without any intersection, add the results without genename
        if new_result == []:
            genename = "-"
            new_result = [pos_bef_PAM, subsequence, genename]
            output += [new_result,]

    return output


def insertion_detector_counter(ins_positions_41, ins_positions_7, unique_inter_results):
    """
    Given the list of insertion positions and the list of unique results it assigns the number of insertions
    InsTh41 and InsTh7 before and after the 17 position and the total
    """

    filehandle_41 = open(ins_positions_41, "r")
    filehandle_7 = open(ins_positions_7, "r")

    def generate_list(filehandle):
        l = []
        for line in filehandle:
            line = line.split()
            l.append(int(line))
        return l

    l_41 = generate_list(filehandle_41)
    l_7  = generate_list(filehandle_7)

    # Start the counting process

    results = []

    for result in unique_inter_results:
        new_result = []
        for pos_bef_PAM, subsequence, genename in results:
            bef17_41 = 0
            aft17_41 = 0
            bef17_7  = 0
            aft17_7  = 0

            for ins_position in l_41:
                if ins_position <= pos_bef_PAM + 2 and ins_position >= pos_bef_PAM - 16:
                    bef17_41 += 1
                if ins_position < pos_bef_PAM - 21 and ins_position < pos_bef_PAM - 16:
                    aft17_41 += 1

            for ins_position in l_7:
                if ins_position <= pos_bef_PAM + 2 and ins_position >= pos_bef_PAM - 16:
                    bef17_7 += 1
                if ins_position < pos_bef_PAM - 21 and ins_position < pos_bef_PAM - 16:
                    aft17_7 += 1

            total41 = bef17_41 + aft17_41
            total7  = bef17_7  + aft17_7
            total = total41 + total7

            new_result = [pos_bef_PAM, subsequence, genename, bef17_41, aft17_41, total41, bef17_7, aft17_7, total7, total]


def file_generator(list_of_lists, fileoutname):
    """
    It generates plain text files containing the list o lists content
    """

    fo = open(fileoutname, "w")

    fo.write("pos_bef_PAM\tsubsequence\t41_bef17\t41_aft\ttotal41\t7_bef17\t7_aft17\ttotal7\ttotal\n")

    for start, sequence, genename in list_of_lists:

        fo.write(str(start)+"\t"+sequence+'\t'+genename+"\n")

    fo.close()

