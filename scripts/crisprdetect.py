#!/usr/bin/env python3

import script1_functions as func


def general_pipeline(genome, motif, length_motif, outputfile, neg_strand = None):
    """
    Generate a file containing all the 20mer sequences and its position:
    """

    # Obtain the indices for the 22G
    indices = func.motif_finder(genome, motif)

    # Return a list of lists having the start of the motif and the sequence
    raw_results = func.subseq_extractor(genome, indices, length_motif, negative_strand = neg_strand)

    # Remove all those sequence that appears more than once in the genome
    unique_results = func.remove_redundant(genome, raw_results)

    # Find sections interrumping genes:
    gene_results = func.interruption("../ref_data/gene_coordinates.txt", unique_results)

    # Generate the output file
    func.file_generator(gene_results, outputfile)


if __name__ == "__main__":

    # Define the genome
    # complete_sequence = func.fasta_processor(filename = "../ref_data/mpn_genome_NC_000912.fasta")
    positive_strand_sequence = func.fasta_processor(filename = "../ref_data/prueba.fasta")
    negative_strand_sequence = func.reverse_complement(positive_strand_sequence)


    # Run the pipelines
    general_pipeline(positive_strand_sequence, "GG", 21, "../results/data1_positive_unique.txt")
    general_pipeline(negative_strand_sequence, "GG", 21, "../results/data1_negative_unique.txt", neg_strand = True)


