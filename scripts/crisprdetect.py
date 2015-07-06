#!/usr/bin/env python3

import script1_functions as func


def general_pipeline(genome, motif, length_motif, outputfile, neg_strand = None):
    """
    Generate a file containing all the 20mer sequences and its position:
    """

    # Informing about the strand
    if not neg_strand:
        print("Starting with the positive strand")
    else:
        print("Now the negative")


    # Obtain the indices for the 22G
    indices = func.motif_finder(genome, motif)


    # Return a list of lists having the start of the motif and the sequence
    print("Generating raw results...")
    raw_results = func.subseq_extractor(genome, indices, length_motif, negative_strand = neg_strand)
    print(str(len(raw_results))+" possible oligos found")

    # Remove all those sequence that appears more than once in the genome
    print("Ensuring they are unique results...")
    unique_results = func.remove_redundant(genome, raw_results)
    print(str(len(unique_results))+" unique oligos found")

    # Find if they are interrumping genes:
    print("Adding information about the genes they can be interrumping...")
    gene_results = func.interruption("../ref_data/gene_coordinates.txt", unique_results, negative_strand = neg_strand)
    print(str(len(gene_results)) + " = Have you lose results?")

    # Find the number of insertion 41 and 7 in the motif:
    # Define InsTh files:
    file41 = "../ref_data/insTh41_positions.txt"
    file7 = "../ref_data/insTh7_positions.txt"

    print("Counting the insertion sites overlapping")
    final_results = func.insertion_detector_counter(file41, file7, gene_results)
    print(str(len(final_results)) + " = Have you lose results?")

    # Generate the output file
    func.file_generator(final_results, outputfile)


if __name__ == "__main__":

    # Define the genome
    # complete_sequence = func.fasta_processor(filename = "../ref_data/mpn_genome_NC_000912.fasta")
    positive_strand_sequence = func.fasta_processor(filename = "../ref_data/mpn_genome_NC_000912.fasta")
    negative_strand_sequence = func.reverse_complement(positive_strand_sequence)

    # Run the pipelines
    general_pipeline(positive_strand_sequence, "GG", 21, "../results/data1_positive_unique.txt")
    general_pipeline(negative_strand_sequence, "GG", 21, "../results/data1_negative_unique.txt", neg_strand = True)


