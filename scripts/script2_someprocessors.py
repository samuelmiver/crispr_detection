#!/usr/bin/env python3

import script1_functions as func


# # Process insertion files:
# func.insertion_processor("../ref_data/NC_000912transposoninsertions.gbk", "InsTh41", "../ref_data/insTh41_positions.txt")
# func.insertion_processor("../ref_data/NC_000912transposoninsertions.gbk", "InsTh7", "../ref_data/insTh7_positions.txt")

# # Process coding  fasta file:
# func.coding_coordinates("../ref_data/coding_mpnM129.txt", "../ref_data/gene_coordinates.txt")

# Order the total file (after concatenation!)
func.order_total_file("../results/total.txt")

# # Convert the gene names file to other having the same coordinates but with the MPN nomenclature:
# func.genename_associator("../ref_data/gene_mpn_relation.gbk", "../ref_data/gene_coordinates.txt")
