# -*- coding: utf-8 -*-
import csv
import argparse
import sys
import logging
from src.oligo_extractor import OligoExtractor
from Bio.Seq import Seq

DATA_DIR_OLIGO = "/Users/ngocht/data/oligos/"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Run the ASO thermodynamics pipeline to retrieve the ddG landscape for "
                    "selective oligos by the Ensemble gene ID of interest"
    )

    parser.add_argument(dest='gene_id',
                        type=str,
                        help="Gene Name to be searched for candidate oligos")
    parser.add_argument(dest='k',
                        type=int,
                        help='Length of oligo'
                        )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Outputs progress information on the console."
    )
    parser.add_argument(
        "--species", "-s",
        type=str,
        default="human",
        help="mouse or human (default)"
    )
    parser.add_argument(
        "--ensembl-release", "-e",
        type=int,
        default=111,
        help="Ensemble release, default=111"
    )

    args = parser.parse_args()
    logging.basicConfig(
        # force=True,
        level=logging.INFO,
        stream=sys.stdout,
        format='### INFO - %(asctime)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logging.info("%s starting up" % sys.argv[0])

    if args.species == "mouse":
        scaffold_path = None
        bowtie_index = f"GRCm38_{args.ensembl_release}"
    elif args.species == "human":
        scaffold_path = f"/Users/ngocht/Library/Caches/pyensembl/GRCh38/ensembl{args.ensembl_release}/Homo_sapiens.GRCh38.{args.ensembl_release}.chr_patch_hapl_scaff.gtf"
        bowtie_index = "GRCh38"
    else:
        raise ValueError("Only mouse or human species implemented.")

    logging.info(args)
    oligo_obj = OligoExtractor(args.gene_id, args.ensembl_release, args.species, args.k, None, scaffold_path)
    oligo_obj.get_candidate_oligos_by_gene()
    oligo_obj.run_bowtie()

    # TODO: make optional
    with open(f"{DATA_DIR_OLIGO}{bowtie_index}_{args.gene_id}_filtered_{args.k}mers_test.csv", "w") as filteredkmerfile:
        writer = csv.writer(filteredkmerfile)
        writer.writerows([[str(Seq(x).reverse_complement())] for x in oligo_obj.filtered_kmers])


    # Sampling kmers from tissues, assumes existence of KMC kmer counts with k length


