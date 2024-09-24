import csv
import subprocess
import shlex
import os
from Bio.SeqUtils import gc_fraction
from pyensembl import EnsemblRelease, Genome
from utils.sequence_analysis import get_chromosomal_positions_per_transcript
from utils.sequence_analysis import get_exon_id, gc_content, longest_at_run, longest_t_run

import logging
import time
import configparser 
import pandas as pd
from Bio.Seq import Seq
import polars as pl


# Create a configparser object
config = configparser.ConfigParser()

# Read the configuration file
config.read('config.ini')



class OligoExtractor:
    """
    A class to extract and analyze oligonucleotide (k-mer) sequences from a specified gene, 
    using data from the Ensembl database and aligning them with Bowtie2.

    This class provides functionalities to:
    - Extract k-mer sequences from a specified gene.
    - Align k-mers using Bowtie2 and analyze alignment results to find viable kmers for ASO Design.
    - Compute result k-mers along with their Intrinsic and Extrinsic feautres.

    Attributes:
        gene_id (str): The Ensembl gene ID for the target gene.
        e_release (int): The Ensembl release version to use.
        g_assembly (str): The genome assembly version (e.g., '38' for GRCh38).
        species (str): The species of interest, either "mouse" or "human".
        k (int): The length of k-mers to extract.
        bowtie_index (str): The path to the Bowtie2 index file.
        gc_bounds (tuple, optional): A tuple specifying the lower and upper GC content bounds for filtering k-mers.
        scaffold_path (str, optional): Path to a GTF file for scaffold annotations.
        filtered_kmers (list): A list to store k-mers that pass all filters.
        bowtie_infile (str): The path to the input FASTA file for Bowtie2.
        ensembl_obj (Genome): An instance of EnsemblRelease for querying gene and transcript data.
        ensembl_obj_scaffolds (Genome, optional): An optional Genome object for scaffold annotations.
        gene (Gene): The Gene object representing the target gene.
        transcript_lookup (dict): A dictionary mapping transcript IDs to gene IDs.
    
    Methods:
        _kmers(s: str) -> set:
            Generate k-mers from the input sequence and optionally filter them by GC content.

        _runCommand(command: str) -> int:
            Execute a shell command and return the exit code.
            
        get_gene_transcript_mapping(save_to_file: str = None) -> dict:
            Create a mapping of transcript IDs to gene information, and optionally save to a file.

        run_bowtie():
            Execute Bowtie2 alignment for the k-mers and filter the aligned k-mers.

        get_candidate_oligos_by_gene():
            Extract k-mers from the target gene, map them to chromosomal positions, and save the candidate oligos to a CSV file.

        get_kmer_occurrences() -> dict:
            Calculate the number of occurrences of each k-mer across the whole Genome.

        store_kmer_results(cofoldOutFile: str):
            Combine the k-mer analysis results and save the final results to a CSV file.
    """

    def __init__(self, gene_id, e_release, g_assembly, species, k, bowtie_index, gc_bounds= None, scaffold_path=None):
        self.gene_id = gene_id
        self.k = k
        self.g_assembly =  g_assembly
        self.filtered_kmers = []
        self.gc_bounds=gc_bounds
        self.bowtie_index = bowtie_index
        self.bowtie_infile = f"{config['DEFAULT']['DataDir']}/bowtie2Home/{self.gene_id}_{self.k}mers.fa"

        if species == "mouse":
            self.species = "mus_musculus"
            # mouse doesn't have scaffold so far...
        elif species == "human":
            self.species = "homo_sapiens"
        else:
            raise ValueError("Only mouse or human species implemented.")
        
        self.ensembl_obj = EnsemblRelease(release=e_release, species=self.species)
        self.ensembl_obj.download()
        self.ensembl_obj.index()
        self.scaffold_path = scaffold_path


        if scaffold_path:
            self.ensembl_obj_scaffolds = Genome(
                reference_name=f'GRCh{g_assembly}',
                annotation_name='scaffolds',
                gtf_path_or_url=scaffold_path,
            )
            self.ensembl_obj_scaffolds.download()
            self.ensembl_obj_scaffolds.index()
        else:
            self.ensembl_obj_scaffolds = None
        

        self.gene = self.ensembl_obj.gene_by_id(gene_id=gene_id)
        
        logging.info(f"Gene name: {self.gene.gene_name}")
        logging.info(f"Build transcript gene references")
        self.transcript_lookup = self.get_gene_transcript_mapping(save_to_file=f"transcript_gene_mapping_GRC{self.species[0]}{g_assembly}.csv")

    def _kmers(self, s):
        """
        Generate k-mers from the input sequence and filter them based on GC bounds if specified.

        Parameters:
            s (str): The input DNA sequence from which k-mers are generated.

        Returns:
            set: A set of tuples, where each tuple contains a k-mer and its starting position in the sequence.
        """
        kmers_list = [(s[i:i + self.k], i+1) for i in range(len(s) - self.k + 1)]
        if self.gc_bounds:
            kmers_list = [seq for seq in kmers_list if self.gc_bounds[0] <= gc_fraction(seq[0]) <= self.gc_bounds[1]]
        kmers_set = set(kmers_list)

        return kmers_set

    def get_gene_transcript_mapping(self, save_to_file=None):
        """
        Create a mapping of transcript IDs to gene information.

        Optionally, save the mapping to a CSV file, including details about each transcript's exons.

        Parameters:
            save_to_file (str, optional): If provided, the path to the file where the mapping will be saved.

        Returns:
            dict: A dictionary mapping transcript IDs to gene information.
        """
        # TODO: might be extended to exon mapping
        transcript_lookup = dict()
        transcripts = self.ensembl_obj.transcripts()
        if self.scaffold_path:
            transcripts.extend(self.ensembl_obj_scaffolds.transcripts())

        if save_to_file:
            file = open(f"{config['DEFAULT']['DataDir']}/{save_to_file}", "w")
        for t in transcripts:
            if t.transcript_id not in transcript_lookup.keys():
                transcript_lookup[t.transcript_id] = t.gene_id
                if save_to_file:
                    for e in t.exons:
                        file.write(
                            f"{t.transcript_id},{t.gene_id},{t.contig}:{t.start},{t.contig}:{t.end},{e.exon_id},{e.start},{e.end},{t.gene_name}\n")
        if save_to_file:
            file.close()
        return transcript_lookup

    def _runCommand(self, command):
        """ Execute command while immediately printing its stdout to outFile and stderr to the terminal (or logging)."""
        return_code = subprocess.call(shlex.split(command))
        return return_code

    def run_bowtie(self, gene_only = False):
        """
        Execute Bowtie2 alignment for the k-mers and filter the aligned k-mers.

        The method runs Bowtie2 with the specified parameters, processes the alignment output.
        """

                
        # Run RNAcofold
        logging.info("Running Bowtie2")
        
        start = time.time()

        if gene_only:
            outFile = os.path.splitext(self.bowtie_infile)[0] + '_local' + ".sam"
            command = f'bowtie2 -x {config["DEFAULT"]["DataDir"]}/bowtie2Home/{self.bowtie_index}_{self.gene_id}_only -U {self.bowtie_infile} -S {outFile} {config["DEFAULT"]["BowtieLocalArgs"]}'
        else:
            outFile = os.path.splitext(self.bowtie_infile)[0] + ".sam"
            command = f'bowtie2 -x {config["DEFAULT"]["DataDir"]}/bowtie2Home/{self.bowtie_index} -U {self.bowtie_infile} -S {outFile} {config["DEFAULT"]["BowtieArgs"]}'

        logging.info("Command: {}".format(command))
        return_code = self._runCommand(command)
        logging.info("Return Code: {}".format(return_code))
        
        end = time.time() - start

        logging.info(f"Bowtie Processing time: {end}")
        
        return outFile


    def get_viable_kmers(self, inFile):
        """
        Filter the aligned k-mers.

        The method processes the alignment output given by bowtie2.
        """
        columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "ALIGN SCORE", "XS", "XN", "XM", "XO", "XG", "EDIT DIST REF", "MISMATCH POS", "YT"]

        align_file = pl.read_csv(inFile ,separator='\t', has_header=False, new_columns=columns)
        
        res = (align_file
               .group_by('SEQ')
               .agg(
                    pl.col('RNAME')
                    .str.split(".")  
                    .list.first().alias('transcript_id')
                    .replace(self.transcript_lookup)
                    .alias('genes'),
                pl.col('QNAME').first().alias('seq_id'),  
                ).with_columns(
                    pl.col('genes').list.set_difference([self.gene_id])
                ).filter(pl.col('genes').list.len() == 0)
                .select([pl.exclude('genes')])  
            )
        
        self.filtered_kmers = res.select(["seq_id", "SEQ"]).to_numpy().tolist()
                            
        logging.info(f"Viable  {self.k}mers candidates after Bowtie: {len(self.filtered_kmers)}")

    def get_candidate_oligos_by_gene(self):
        """
        Extract and process candidate oligos (k-mers) from a specified gene and save results to files.

        This method identifies candidate oligos of length `k` from the gene transcripts associated with the gene id. 
        It calculates their absolute chromosomal positions, associates them with the relevant transcripts and exons, and then 
        saves the results to a CSV file. Additionally, it creates a Bowtie input file for sequence alignment.

        The method performs the following steps:
        1. Retrieves the transcripts for the specified gene.
        2. Extracts k-mers from each transcript and computes their chromosomal positions, transcript IDs, and exon IDs.
        3. Aggregates k-mer data into a DataFrame, grouping by sequence and chromosomal position.
        4. Saves the k-mer data to a CSV file with custom indexing `SXXXXXX`.
        5. Creates a Bowtie input file with the k-mer sequences formatted for alignment.

        Parameters:
            None: The method uses instance variables `self.gene_id`, `self.k`, `self.ensembl_obj`, `self.ensembl_obj_scaffolds`, 
                and `self.bowtie_infile`.

        Returns:
            None: The method saves results to a CSV file and writes to a Bowtie input file.

        Notes:
            - The CSV file is saved in a directory specified by the configuration file, with a filename based on the gene ID 
            and k-mer length.
            - The Bowtie input file is created with k-mers in FASTA format, suitable for sequence alignment tools.
        """
        logging.info(f"Extract {self.k}mers from gene {self.gene_id}")

        transcripts = self.gene.transcripts
        candidate_oligos = set()
        for t in transcripts:
            # rev_comp_t = Seq(t.sequence).reverse_complement()
            # TODO: make GC content bounds a parameter
            kmers_set = self._kmers(t.sequence)
            
            kmers_set = {(tup[0], 
                          get_chromosomal_positions_per_transcript(t.transcript_id, tup[1], self.ensembl_obj, self.k, self.ensembl_obj_scaffolds), 
                          t.transcript_id,
                          get_exon_id(tup[1], t)) for tup in kmers_set}
            
            candidate_oligos.update(kmers_set)
            
        logging.info(f"{len(candidate_oligos)} candidate {self.k}mers found")
        
        
        columns = ['seq', 'chromosomal_position', 'transcripts', 'exons']
        candidate_oligos = pd.DataFrame(columns=columns, data=candidate_oligos)
        
        candidate_oligos = candidate_oligos.groupby(['seq', 'chromosomal_position']).agg({
            'exons' : lambda x: list(x),
            'transcripts' : lambda x: list(x),
        }).reset_index()
        
        custom_index = [f'S{str(i).zfill(6)}' for i in range(1, len(candidate_oligos) + 1)]
        candidate_oligos.index = custom_index
        
        candidate_oligos.to_csv(f'{config["DEFAULT"]["DataDir"]}/oligos/{self.gene_id}_{self.k}mer_candidates.csv')
        
        os.makedirs(f'{config["DEFAULT"]["DataDir"]}/bowtie2Home', exist_ok=True)
        
        with open(self.bowtie_infile, "w") as tmp_bowtie_in:
            candidate_oligos.apply(lambda x: tmp_bowtie_in.write(">" + str(x.name) + "\n" + x['seq'] + "\n"), axis = 1)

        
    
       
    def get_kmer_occurrences(self, inFile):
        """
        Calculate the number of unique chromosomal position occurences in the Gene for each k-mer based on local Bowtie2 alignment data.

        This method processes the local Bowtie2 alignment results to determine the number of unique absolute chromosomal positions 
        associated with each k-mer. It reads a SAM file generated by Bowtie2, groups the data by k-mer sequence, 
        and calculates the number of unique chromosomal positions for each group.

        The method performs the following steps:
        1. Reads the SAM file into a DataFrame.
        2. Groups the data by k-mer sequence and calculates the number of unique chromosomal positions for each group.
        3. Stores the results in a dictionary where keys are k-mer sequences and values are the number of unique positions.
        
        The output dictionary is also stored as occurrence_dictionary in the object.

        Returns:
            dict: A dictionary where keys are k-mer sequences and values are the number of unique chromosomal positions.
                The dictionary provides the count of unique positions for each k-mer based on the alignment results.

        """
        logging.info(f"Calculating number of kmer occurrences")

        
        def calculate_occurrences(group):
            # Extract relevant columns
            positions = group.apply(lambda row: get_chromosomal_positions_per_transcript(row['RNAME'], row['POS'], self.ensembl_obj, int(config["DEFAULT"]["MinConsecutiveMatches"]), self.ensembl_obj_scaffolds), axis=1)
            unique_positions = positions.unique()
            return unique_positions
    
        columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "ALIGN SCORE", "XS", "XN", "XM", "XO", "XG", "EDIT DIST REF", "MISMATCH POS", "YT"]

        sam_out = pd.read_csv(inFile, sep="\t", header=None, names = columns)
        # Group by SEQ
        sam_out_agg = sam_out.groupby('QNAME').apply(calculate_occurrences).reset_index(name='UniquePositions')
        print(sam_out_agg.head())
        # Rename columns for clarity (optional)
        sam_out_agg.columns = ['QNAME', 'UniquePositions']
        
        # Convert to dictionary
        self.occurrence_dictionary = sam_out_agg.set_index('QNAME').to_dict()['UniquePositions']
        
        return self.occurrence_dictionary
    
    def store_kmer_results(self, cofoldOutFile): 
        """
        Generate a CSV file with detailed results for each k-mer, including various properties and metrics.

        This method processes k-mer candidates and RNAcofold output to compile a comprehensive results file. It merges 
        information about each k-mer's sequence, its reverse complement, GC content, longest AT-run, longest T-run, 
        chromosomal positions, associated transcripts and exons in the Gene, multiplicity, and RNAcofold binding energy. The results 
        are saved to a CSV file.

        Parameters:
            cofoldOutFile (str): The path to the RNAcofold output file in CSV format, containing binding energy information 
                                for the k-mers.

        Returns:
            None: The results are saved directly to a CSV file in the specified directory. No value is returned by the method.
        """
        logging.info(f"Completing final results")

        cofold_out = pd.read_csv(cofoldOutFile)
        cofold_out.set_index('seq_id', inplace=True)

        oligo_candidates = pd.read_csv(f'{config["DEFAULT"]["DataDir"]}/oligos/{self.gene_id}_{self.k}mer_candidates.csv', index_col=0)

        os.makedirs(f"{config['DEFAULT']['DataDir']}/oligos", exist_ok=True)
        
        # result csv column names
        columns = ['seq_num',  
                   'oligo_reverse_comp', 
                   'oligo_gc_content',
                   'oligo_longest_at_run',
                   'oligo_longest_t_run',
                   'target', 
                   'absolute_loc', 
                   'ordered_transcripts', 
                   'ordered_exons',
                   'prone_multiplicity', 
                   'dG_binding']
        
        kmer_indices = [x[0] for x in self.filtered_kmers]
        res_temp = []
        
        for idx in kmer_indices:
            can = oligo_candidates.loc[idx]
            res_temp.append((idx,                                             # seq_num
                             can['seq'],                                      # oligo_reverse_comp
                             gc_content(can['seq']),                          # oligo_gc_content
                             longest_at_run(can['seq']),                      # oligo_longest_at_run
                             longest_t_run(can['seq']),                       # oligo_longest_t_run
                             str(Seq(can['seq']).reverse_complement()),       # target
                             can['chromosomal_position'],                     # absolute_loc
                             can['transcripts'],                              # ordered_transcripts
                             can['exons'],                                    # ordered_exons
                             len(self.occurrence_dictionary.get(idx, 0)),          # prone_multiplicity
                             cofold_out.loc[idx]['dG_binding'])               # dG_binding
                            )
            
        kmer_results = pd.DataFrame(res_temp, columns=columns)
        
        kmer_results.set_index('seq_num', inplace=True)
        kmer_results.to_csv(f'{config["DEFAULT"]["DataDir"]}/oligos/{self.gene_id}_{self.k}mer_results.csv')