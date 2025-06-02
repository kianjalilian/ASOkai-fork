import os
import gzip
import logging
from typing import Optional, List, Tuple, Dict, Any, Union
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from src.utils.genome_utils import Genome
import gget
import urllib.request, urllib.parse
from src.utils.time_utils import timed, format_duration

# Assuming Gene and Transcript are part of genome_utils or a similar accessible module
# Adjust the import based on your actual structure for Gene and Transcript classes
try:
    from src.utils.genome_utils import Genome, Gene
except ImportError:
    logging.warning("genome_utils.Gene or not found, using placeholder. Functionality might be limited.")
    # Placeholder if Gene class is not directly importable for type hinting
    class Gene: pass 


class GenomeDownloader:
    def __init__(self, species: str, e_release: int, genome_dir: str):
        """
        Initializes the GenomeDownloader.

        Args:
            species (str): Species name. e.g., 'homo_sapiens'.
            e_release (int): The Ensembl release version.
            genome_dir (str): Path to the directory where the genome files will be saved.
        """
        self.species = species
        self.e_release = e_release
        self.genome_dir = genome_dir
        os.makedirs(self.genome_dir, exist_ok=True)

    def download(self) -> Tuple[str, str, str, Optional[str]]:
        """
        Downloads the genome files for the specified species.
        
        Returns:
            Tuple[str, str, str, Optional[str]]: Paths to the GTF, cDNA, genome FASTA files, and an optional scaffold GTF path.
        """
        gtf_url, cdna_url, genome_url = tuple(
            gget.ref(self.species, which=["gtf", "cdna", "dna"], release=self.e_release, ftp=True, verbose=False)
        )

        gtf_name = os.path.basename(urllib.parse.urlparse(gtf_url).path)
        cdna_name = os.path.basename(urllib.parse.urlparse(cdna_url).path)
        genome_name = os.path.basename(urllib.parse.urlparse(genome_url).path)

        scaffold_gtf_path = None
        if self.species == 'homo_sapiens':     
            scaffold_gtf_url = gtf_url.replace('.gtf.gz', '.chr_patch_hapl_scaff.gtf.gz')
            
            scaffold_gtf_name = os.path.basename(urllib.parse.urlparse(scaffold_gtf_url).path)
            scaffold_gtf_path = os.path.join(self.genome_dir, scaffold_gtf_name)
            if not os.path.exists(scaffold_gtf_path):
                logging.info(f"Downloading scaffold GTF file to '{scaffold_gtf_path}'")
                urllib.request.urlretrieve(scaffold_gtf_url, scaffold_gtf_path)
            else:
                logging.info(f"Scaffold GTF file already exists at '{scaffold_gtf_path}'")
                
        gtf_path = os.path.join(self.genome_dir, gtf_name)
        cdna_path = os.path.join(self.genome_dir, cdna_name)
        genome_path = os.path.join(self.genome_dir, genome_name)

        if not os.path.exists(gtf_path):
            logging.info(f"Downloading GTF file to '{gtf_path}'")
            urllib.request.urlretrieve(gtf_url, gtf_path)
        else:
            logging.info(f"GTF file already exists at '{gtf_path}'")
            
        if not os.path.exists(cdna_path):
            logging.info(f"Downloading cDNA file to '{cdna_path}'")
            urllib.request.urlretrieve(cdna_url, cdna_path)
        else:
            logging.info(f"cDNA file already exists at '{cdna_path}'")
            
        if not os.path.exists(genome_path):
            logging.info(f"Downloading genome FASTA file to '{genome_path}'")
            urllib.request.urlretrieve(genome_url, genome_path)
        else:
            logging.info(f"Genome FASTA file already exists at '{genome_path}'")

        return gtf_path, cdna_path, genome_path, scaffold_gtf_path


class GenomeDataManager:
    def __init__(self,
                 gene_id: str,
                 species: str,
                 e_release: int,
                 genome_assembly: int,
                 genome_dir: str,
                 tsl_config_str: Optional[str] = None,
                 force_overwrite: bool = False
                 ):
        """
        Manages downloading, processing, and providing paths to essential genome data files.
        """
        self.gene_id = gene_id
        self.species = species
        self.e_release = e_release
        self.genome_assembly = genome_assembly
        self.genome_dir = genome_dir
        self.force_overwrite = force_overwrite

        logging.info(f"Initializing GenomeDataManager for gene: {self.gene_id}, species: {self.species}, release: {self.e_release}")

        self.tsl_is_active, self.tsl_list_to_keep = self._convert_tsl_list(tsl_config_str) if tsl_config_str else (False, None)
        if self.tsl_is_active:
            logging.info(f"TSL filtering is active. Levels to keep: {self.tsl_list_to_keep}")
        else:
            logging.info("TSL filtering is not active or all levels are included.")

        # 1. Download base files using GenomeDownloader
        downloader = GenomeDownloader(
            species=self.species,
            e_release=self.e_release,
            genome_dir=self.genome_dir,
        )
        # These paths are to the raw downloaded files as returned by GenomeDownloader
        self.raw_gtf_path, self.raw_cdna_path, self.raw_genome_path, self.raw_scaffold_gtf_path = downloader.download()

        # 2. Initialize Genome objects
        # Construct reference name for genome-utils
        # Example: GRCm38 for mus_musculus, GRCh38 for homo_sapiens
        if self.species.lower() == "mus_musculus":
            self.species = "Mus_musculus"
            reference_prefix = "GRCm"
        elif self.species.lower() == "homo_sapiens":
            reference_prefix = "GRCh"
        else:
            # Fallback or error for other species, assuming genome-utils needs this convention
            reference_prefix = f"GRC{self.species[0].upper()}" 
            logging.warning(f"Species '{self.species}' is not 'mus_musculus' or 'homo_sapiens'. Using prefix '{reference_prefix}'. Ensure this matches genome-utils expectations.")
        
        main_reference_name = f'{reference_prefix}{self.genome_assembly}'
        self.genome_file_prefix = f'{self.species}.{reference_prefix}{self.genome_assembly}'
        self.genome = Genome(
            reference_name=main_reference_name,
            annotation_version=str(self.e_release), # genome-utils might expect string
            gtf_path=self.raw_gtf_path,
            transcript_fasta_paths=self.raw_cdna_path,
            primary_assembly_path=self.raw_genome_path,
            tsl_to_keep=self.tsl_list_to_keep if self.tsl_is_active else None,
            protein_coding_only=True,
        )
        # Index the genome object after initialization (will apply TSL filtering from GTF)
        self.genome.index()
        logging.info(f"Main genome ('{main_reference_name}') indexed. Found {len(self.genome.genes)} genes.")

        self.genome_scaffolds: Optional[Genome] = None
        if self.raw_scaffold_gtf_path:
            # Assuming scaffold reference name might be different or need a suffix
            scaffold_reference_name = f'{main_reference_name}_scaffolds' 
            self.genome_scaffolds = Genome(
                reference_name=scaffold_reference_name,
                gtf_path=self.raw_scaffold_gtf_path,
            )
            self.genome_scaffolds.index(overwrite=False)
            logging.info(f"Scaffold genome ('{scaffold_reference_name}') indexed.")

        self.target_gene: Gene = self.genome.gene_by_id(self.gene_id)
        if not self.target_gene:
            # Try scaffolds if not in primary assembly, though typically genes are on primary
            if self.genome_scaffolds:
                logging.warning(f"Target gene {self.gene_id} not found in primary assembly, checking scaffolds...")
                self.target_gene = self.genome_scaffolds.gene_by_id(self.gene_id)
            if not self.target_gene: # Still not found
                raise ValueError(f"Target gene {self.gene_id} not found in the main or scaffold genome assemblies.")
        logging.info(f"Target gene '{self.target_gene.gene_name if hasattr(self.target_gene, 'gene_name') else self.gene_id}' loaded.")

        # --- File path attributes for prepared files will be populated by _prepare_xxx methods ---
        self.processed_cdna_path: str = self.raw_cdna_path
        self.cdna_excluding_target_path: Optional[str] = None
        self.genes_pre_mrna_fasta_path: Optional[str] = None
        self.target_gene_transcripts_fasta_path: Optional[str] = None
        
        # --- Run preparation methods ---
        self._prepare_tsl_filtered_cdna()       # Sets self.processed_cdna_path
        self._prepare_pre_mrna_fastas(genes_to_exclude=[self.gene_id], force=self.force_overwrite) 
        self._prepare_cdna_fasta_excluding_target()   # Sets self.cdna_excluding_target_path
        self._prepare_target_gene_transcripts_fasta()         # For specific transcript operations

        logging.info("GenomeDataManager initialization complete.")

    def _convert_tsl_list(self, tsl_str: str) -> Tuple[bool, Optional[List[Optional[int]]]]:
        """
        Convert a string of transcript support levels into a list of integers or None values.
        
        Args:
            tsl_str (str): Comma-separated string of transcript support levels (e.g., "tsl1,tsl2,tslNA")
            
        Returns:
            Tuple[bool, Optional[List[Optional[int]]]]: A tuple containing:
                - bool: True if specific TSLs were provided, False if all TSLs are included
                - Optional[List[Optional[int]]]: List of TSL values (1-5 or None for NA)
        """
        tsl_tokens = [token.strip() for token in tsl_str.split(',') if token.strip()]
        converted_tsls = []

        for token in tsl_tokens:
            if token.lower().startswith("tsl"):
                value_part = token[3:]
                if value_part.lower() == "na":
                    converted_tsls.append(None)
                else:
                    try:
                        converted_tsls.append(int(value_part))
                    except ValueError:
                        logging.warning("Invalid transcript support level '%s'; using None.", token)
                        converted_tsls.append(None)
            else:
                try:
                    converted_tsls.append(int(token))
                except ValueError:
                    logging.warning("Invalid transcript support level '%s'; using None.", token)
                    converted_tsls.append(None)

        all_tsls = [1, 2, 3, 4, 5, None]
        if converted_tsls == all_tsls:
            return False, None
        return True, converted_tsls
    
    
    def _get_gene_transcript_mapping(self) -> Dict[str, str]:
        """
        Create a mapping of transcript IDs to gene information.
        Uses the main genome and scaffold genome objects if available.

        Returns:
            A dictionary mapping transcript IDs to gene information
        """
        all_transcripts = []
        if self.genome:
            all_transcripts.extend(self.genome.transcripts) # Use .transcripts property
        if self.genome_scaffolds:
            all_transcripts.extend(self.genome_scaffolds.transcripts) # Use .transcripts property
        
        mapping = {}
        for t in all_transcripts:
            if hasattr(t, 'transcript_id') and hasattr(t, 'gene_id'):
                 mapping[t.transcript_id] = t.gene_id
            else:
                logging.debug(f"Transcript object missing transcript_id or gene_id: {t}")
        return mapping

    def _prepare_tsl_filtered_cdna(self):
        """
        If TSL filtering is active, this method now primarily ensures a TSL-filtered
        cDNA FASTA file is written out if needed, using the already TSL-filtered
        transcripts from the Genome object.
        It updates self.processed_cdna_path to point to this (potentially new) file.
        The Genome object itself is ALREADY filtered by TSL upon initialization if tsl_is_active.
        """
        if not self.tsl_is_active:
            self.processed_cdna_path = self.raw_cdna_path
            return

        
        base_raw_cdna_name = os.path.basename(self.raw_cdna_path)
        tsl_suffix = f".tsl{'_'.join(str(tsl) if tsl is not None else 'NA' for tsl in self.tsl_list_to_keep)}"
        
        if base_raw_cdna_name.endswith(".fa.gz"):
            tsl_filtered_fasta_name = base_raw_cdna_name.replace(".fa.gz", f"{tsl_suffix}.fa.gz")
        elif base_raw_cdna_name.endswith(".fasta.gz"):
            tsl_filtered_fasta_name = base_raw_cdna_name.replace(".fasta.gz", f"{tsl_suffix}.fasta.gz")
        elif base_raw_cdna_name.endswith(".fa"):
            tsl_filtered_fasta_name = base_raw_cdna_name.replace(".fa", f"{tsl_suffix}.fa")
        elif base_raw_cdna_name.endswith(".fasta"):
            tsl_filtered_fasta_name = base_raw_cdna_name.replace(".fasta", f"{tsl_suffix}.fasta")
        else:
            tsl_filtered_fasta_name = f"{base_raw_cdna_name}{tsl_suffix}.fasta"
            logging.warning(f"Raw cDNA path '{self.raw_cdna_path}' has unknown FASTA extension. Outputting TSL-filtered as '{tsl_filtered_fasta_name}'")

        potential_path = os.path.join(self.genome_dir, tsl_filtered_fasta_name)

        if os.path.exists(potential_path) and os.path.getsize(potential_path) > 0:
            logging.info(f"TSL-filtered cDNA FASTA already exists: {potential_path}")
            self.processed_cdna_path = potential_path
            # Update Genome's transcript fasta path to the filtered one and re-parse sequences from it.
            # This is important if the existing file should be the source of sequences.
            if self.genome:
                self.genome.transcript_fasta_paths = [self.processed_cdna_path]
                logging.info(f"Genome object will use existing TSL-filtered FASTA: {self.processed_cdna_path} for transcript sequences.")
                self.genome._parse_transcript_fasta() # Re-parse sequences from this specific file
            return
        
        logging.info("Preparing TSL-filtered cDNA FASTA output file.")

        count_written = 0
        try:
            open_func = gzip.open if potential_path.endswith(".gz") else open
            mode = 'wt' if potential_path.endswith(".gz") else 'w'
            with open_func(potential_path, mode) as f_out:
                # Iterate TSL-filtered transcripts from the already filtered Genome object
                # The Genome object's transcripts are already filtered by TSL during its initialization.
                for transcript in self.genome.transcripts: # These are already TSL and biotype filtered
                    # We still need to ensure the transcript sequence is loaded from the *raw* cDNA
                    # if it hasn't been, before writing. The genome.index() called earlier handles this.
                    if transcript.sequence:
                        f_out.write(f">{transcript.transcript_id}\n{transcript.sequence}\n")
                        count_written += 1
                    else:
                        # This case implies the transcript was in the GTF (and passed TSL/biotype filters)
                        # but its sequence was not found in the raw_cdna_path.
                        logging.debug(f"Transcript {transcript.transcript_id} (gene {transcript.gene_id}) passed GTF filters but has no sequence from {self.raw_cdna_path}. Not written to TSL-filtered FASTA.")
                        
            
            if count_written > 0:
                self.processed_cdna_path = potential_path
                logging.info(f"TSL-filtered cDNA FASTA created with {count_written} transcripts: {self.processed_cdna_path}")
                # Ensure Genome object uses this newly created TSL-filtered FASTA for its sequences
                if self.genome:
                    self.genome.transcript_fasta_paths = [self.processed_cdna_path]
                    logging.info(f"Genome object updated to use newly created TSL-filtered FASTA: {self.processed_cdna_path}")
                    self.genome._parse_transcript_fasta() # Re-parse sequences
            else:
                logging.warning(f"TSL filtering active, but no transcripts were written to {potential_path}. "
                                "This might be due to no transcripts passing filters or sequences not found in raw cDNA. "
                                f"Processed cDNA path remains {self.raw_cdna_path}.")
                self.processed_cdna_path = self.raw_cdna_path
                
        except Exception as e:
            logging.error(f"Error during TSL filtering of cDNA: {e}. Processed cDNA path set to raw path: {self.raw_cdna_path}.")
            self.processed_cdna_path = self.raw_cdna_path
            if self.genome:
                self.genome.transcript_fasta_paths = [self.raw_cdna_path] # Ensure it points to raw
                logging.info("Reverting Genome object to use raw cDNA path for transcript sequences due to error.")
                self.genome._parse_transcript_fasta()


    def _prepare_pre_mrna_fastas(self, genes_to_exclude: Optional[List[str]] = None, force: bool = False):
        """
        Prepares a pre-mRNA FASTA file for genes, potentially excluding specified ones.
        The path to the generated FASTA file is stored in self.genes_pre_mrna_fasta_path.

        Args:
            genes_to_exclude: A list of gene IDs to exclude from the FASTA file.
                              If None or empty, all genes will be included.
            force: If True, overwrite the FASTA file if it already exists.
                   If False and the file exists, generation is skipped.
        """
        logging.info("Preparing pre-mRNA FASTA file...")

        base_name = f"{self.genome_file_prefix}.premrna"
        log_message_exclusion_details = ""

        if genes_to_exclude:
            if len(genes_to_exclude) == 1:
                excluded_gene_id = genes_to_exclude[0]
                suffix = f".all_except_{excluded_gene_id}.fa"
                log_message_exclusion_details = f" (excluding gene '{excluded_gene_id}')"
            else:
                suffix = f".all_except_{len(genes_to_exclude)}_genes.fa"
                log_message_exclusion_details = f" (excluding {len(genes_to_exclude)} genes)"
        else:  # No genes to exclude
            suffix = ".all.fa"
            log_message_exclusion_details = " (all genes)"
        
        genes_fasta_name = f"{base_name}{suffix}"
        # Assuming self.genome_dir is the correct directory for this genome-wide FASTA.
        target_fasta_path = os.path.join(self.genome_dir, genes_fasta_name)

        logging.info(f"Target pre-mRNA FASTA path: {target_fasta_path}{log_message_exclusion_details}")

        if os.path.exists(target_fasta_path) and not force:
            logging.info(
                f"Pre-mRNA FASTA file already exists at {target_fasta_path} and force is False. Skipping generation."
            )
            self.genes_pre_mrna_fasta_path = target_fasta_path  # Ensure attribute is set
            return

        try:
            os.makedirs(self.genome_dir, exist_ok=True)

            records_batch = []
            batch_size = 1000  # Adjust batch size as needed
            sequences_written_count = 0
            
            # Open in 'w' mode to create/overwrite.
            # The 'force' logic above determines if we reach this point when file exists.
            with open(target_fasta_path, 'w') as output_handle:
                for seq_record in self.genome.yield_premrna_seqrecords(exclude_genes=genes_to_exclude):
                    records_batch.append(seq_record)
                    sequences_written_count += 1
                    if len(records_batch) >= batch_size:
                        SeqIO.write(records_batch, output_handle, "fasta")
                        records_batch = []
                
                if records_batch:  # Write any remaining records in the last batch
                    SeqIO.write(records_batch, output_handle, "fasta")
            
            if sequences_written_count > 0:
                logging.info(
                    f"Successfully generated pre-mRNA FASTA with {sequences_written_count} sequences: {target_fasta_path}"
                )
            else:
                logging.info(
                    f"Generated pre-mRNA FASTA (0 sequences, possibly all relevant genes were excluded or no genes found): {target_fasta_path}"
                )
            self.genes_pre_mrna_fasta_path = target_fasta_path

        except FileNotFoundError as fnf_error:
            # This typically means the primary assembly FASTA was not found by the Genome class
            logging.error(f"Error preparing pre-mRNA FASTA: {fnf_error}. Ensure primary assembly path is correct.")
            self.genes_pre_mrna_fasta_path = None
        except Exception as e:
            logging.error(
                f"An unexpected error occurred while preparing pre-mRNA FASTA for {target_fasta_path}: {e}"
            )
            self.genes_pre_mrna_fasta_path = None
            # Clean up potentially incomplete file
            if os.path.exists(target_fasta_path):
                try:
                    os.remove(target_fasta_path)
                    logging.info(f"Removed potentially incomplete file: {target_fasta_path}")
                except OSError as oe:
                    logging.error(f"Error removing incomplete file {target_fasta_path}: {oe}")

    def _prepare_cdna_fasta_excluding_target(self):
        """
        Creates a cDNA FASTA file that includes all transcripts *except* those from the target gene.
        Uses the (potentially TSL-filtered) self.processed_cdna_path as the reference for transcript sequences.
        """
        
        target_gene_transcript_ids = {t.transcript_id for t in self.target_gene.transcripts if hasattr(t, 'transcript_id')}
        
        # Simplify name generation and ensure .gz is preserved if input was gzipped
        input_is_gzipped = self.processed_cdna_path.endswith(".gz")
        name_part = os.path.basename(self.processed_cdna_path)
        if input_is_gzipped:
            name_part = name_part[:-3] # Remove .gz
        
        if name_part.endswith(".fa"):
            name_part = name_part[:-3]
            ext = ".fa.gz" if input_is_gzipped else ".fa"
        elif name_part.endswith(".fasta"):
            name_part = name_part[:-6]
            ext = ".fasta.gz" if input_is_gzipped else ".fasta"
        else: # Fallback, just append
            ext = ".fa.gz" if input_is_gzipped else ".fa"

        cdna_excluding_target_fasta_name = f"{name_part}_no_{self.gene_id}{ext}"
        self.cdna_excluding_target_path = os.path.join(self.genome_dir, cdna_excluding_target_fasta_name)

        # Check if file already exists
        if os.path.exists(self.cdna_excluding_target_path) and os.path.getsize(self.cdna_excluding_target_path) > 0:
            logging.info(f"cDNA FASTA excluding target gene already exists: {self.cdna_excluding_target_path}")
            return
        
        logging.info(f"Preparing cDNA FASTA excluding target gene: {self.gene_id}, using input: {self.processed_cdna_path}")

        count_written = 0
        try:
            open_func_out = gzip.open if self.cdna_excluding_target_path.endswith(".gz") else open
            mode_out = 'wt' if self.cdna_excluding_target_path.endswith(".gz") else 'w'

            with open_func_out(self.cdna_excluding_target_path, mode_out) as f_out:
                # Iterate all transcripts from the main genome object.
                # These sequences come directly from the Genome object, which should have handled
                # the TSL filtering if self.processed_cdna_path was updated and genome re-indexed.
                for transcript in self.genome.transcripts:
                    if transcript.transcript_id in target_gene_transcript_ids:
                        continue 
                    
                    # The TSL check here is removed as self.genome.transcripts should already be TSL-filtered.
                    if transcript and hasattr(transcript, 'sequence') and transcript.sequence:
                        f_out.write(f">{transcript.transcript_id}\n{transcript.sequence}\n")
                        count_written += 1
            
            if count_written > 0:
                logging.info(f"Created cDNA FASTA excluding target gene '{self.gene_id}' with {count_written} transcripts: {self.cdna_excluding_target_path}")
            else:
                logging.warning(f"No transcripts written to cDNA FASTA excluding target gene '{self.gene_id}'. File might be empty: {self.cdna_excluding_target_path}")

        except Exception as e:
            logging.error(f"Error creating cDNA FASTA excluding target gene: {e}")
            self.cdna_excluding_target_path = None


    def _prepare_target_gene_transcripts_fasta(self):
        """
        Prepares a FASTA file containing all transcripts for the target gene.
        Uses transcripts from self.target_gene, which are already TSL and biotype filtered.
        The path is stored in self.target_gene_transcripts_fasta_path.
        """
        if not self.target_gene or not hasattr(self.target_gene, 'transcripts') or not self.target_gene.transcripts:
            logging.warning(f"Target gene {self.gene_id} has no transcripts. Skipping target gene transcript FASTA creation.")
            self.target_gene_transcripts_fasta_path = None
            return

        gene_id_safe = self.gene_id.replace(":", "_").replace(" ", "_")
        # Use the processed_cdna_path as a base for naming to reflect TSL filtering if applied
        base_cnda_name_part = os.path.basename(self.processed_cdna_path)
        if base_cnda_name_part.endswith(".fa.gz"):
            base_cnda_name_part = base_cnda_name_part[:-len(".fa.gz")]
        elif base_cnda_name_part.endswith(".fasta.gz"):
            base_cnda_name_part = base_cnda_name_part[:-len(".fasta.gz")]
        elif base_cnda_name_part.endswith(".fa"):
            base_cnda_name_part = base_cnda_name_part[:-len(".fa")]
        elif base_cnda_name_part.endswith(".fasta"):
            base_cnda_name_part = base_cnda_name_part[:-len(".fasta")]
        
        # Construct a meaningful name
        # Example: Mus_musculus.GRCm39.cdna.tsl1_tslNA.ENSMUSG00000020275_transcripts.fa
        output_filename = f"{base_cnda_name_part}.{gene_id_safe}_transcripts.fa"
        if self.processed_cdna_path.endswith(".gz"): # Preserve compression if original was gzipped
            output_filename += ".gz"
        
        self.target_gene_transcripts_fasta_path = os.path.join(self.genome_dir, output_filename)

        if os.path.exists(self.target_gene_transcripts_fasta_path) and not self.force_overwrite:
            logging.info(f"Target gene transcript FASTA already exists: {self.target_gene_transcripts_fasta_path}. Skipping.")
            return

        logging.info(f"Preparing target gene transcript FASTA: {self.target_gene_transcripts_fasta_path}")
        count_written = 0
        try:
            open_func = gzip.open if output_filename.endswith(".gz") else open
            mode = 'wt' if output_filename.endswith(".gz") else 'w'
            
            with open_func(self.target_gene_transcripts_fasta_path, mode) as f_out:
                for transcript in self.target_gene.transcripts:
                    if transcript.sequence:
                        # Create a SeqRecord for proper FASTA formatting if desired, or simple write
                        # record = SeqRecord(Seq(transcript.sequence), id=transcript.transcript_id, description=f"gene_id={self.gene_id}")
                        # SeqIO.write(record, f_out, "fasta-2line") # fasta-2line for no wrapping
                        f_out.write(f">{transcript.transcript_id} gene_id={self.gene_id}\n{transcript.sequence}\n")
                        count_written += 1
                    else:
                        logging.debug(f"Transcript {transcript.transcript_id} for target gene {self.gene_id} has no sequence. Not written.")
            
            if count_written > 0:
                logging.info(f"Target gene transcript FASTA created with {count_written} transcripts: {self.target_gene_transcripts_fasta_path}")
            else:
                logging.warning(f"No transcripts written for target gene {self.gene_id} to {self.target_gene_transcripts_fasta_path}. File might be empty.")

        except Exception as e:
            logging.error(f"Error creating target gene transcript FASTA: {e}")
            self.target_gene_transcripts_fasta_path = None # Indicate failure
            if os.path.exists(output_filename):
                try:
                    os.remove(output_filename)
                except OSError as oe:
                    logging.error(f"Error removing incomplete file {output_filename}: {oe}")

    # --- Getter methods ---
    def get_target_gene_object(self) -> Gene:
        return self.target_gene

    def get_main_genome_object(self) -> Genome:
        return self.genome

    def get_scaffold_genome_object(self) -> Optional[Genome]:
        return self.genome_scaffolds

    def get_processed_cdna_path(self) -> str:
        """Path to cDNA FASTA, potentially TSL-filtered."""
        return self.processed_cdna_path

    def get_cdna_excluding_target_path(self) -> Optional[str]:
        """Path to cDNA FASTA excluding target gene's transcripts."""
        return self.cdna_excluding_target_path

    def get_genes_pre_mrna_fasta_path(self) -> Optional[str]:
        """Path to pre-mRNA FASTA for all genes excluding the target."""
        return self.genes_pre_mrna_fasta_path

    def get_transcript_gene_mapping(self) -> Dict[str, str]:
        """
        Returns a dictionary mapping transcript IDs to gene IDs.
        """
        return self._get_gene_transcript_mapping()
    
    def get_raw_download_paths(self) -> Dict[str, Optional[str]]:
        return {
            "gtf": self.raw_gtf_path,
            "cdna": self.raw_cdna_path,
            "genome": self.raw_genome_path,
            "scaffold_gtf": self.raw_scaffold_gtf_path
        }

    def get_target_gene_transcripts_fasta_path(self) -> Optional[str]:
        """Path to FASTA file containing transcripts of the target gene."""
        return self.target_gene_transcripts_fasta_path




