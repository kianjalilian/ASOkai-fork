import os
import subprocess
import tempfile
import gzip
import logging
import shutil
from typing import Dict, List, Set, Tuple, Optional, Iterable, Any
import concurrent.futures
import threading

from src.utils.time_utils import ProgressTracker

try:
    import polars as pl
except ImportError:
    logging.warning("Polars library not found. calculate_per_gene_counts_matrix will not be available.")
    pl = None

# Configure basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class KmerData:
    """Manages and prepares k-mer sets derived from potential secondary sites."""
    def __init__(self, potential_secondary_sites: Dict[str, List[Tuple[str, float]]], k: int):
        self.k = k
        self.aso_ids: List[str] = list(potential_secondary_sites.keys())
        self.potential_kmers_by_aso: Dict[str, List[str]] = {
            aso_id: [site[0] for site in sites]
            for aso_id, sites in potential_secondary_sites.items()
        }

        self.all_unique_potential_kmers: Set[str] = set()
        for sites_list in self.potential_kmers_by_aso.values():
            for kmer_seq in sites_list:
                if len(kmer_seq) != self.k:
                    raise ValueError(
                        f"K-mer '{kmer_seq}' has length {len(kmer_seq)}, "
                        f"but k is set to {self.k}. All potential secondary "
                        "site k-mers must have the same length k."
                    )
                self.all_unique_potential_kmers.add(kmer_seq)
        
        if not self.all_unique_potential_kmers:
            logging.warning("KmerData: No potential k-mers found from potential_secondary_sites input.")
        logging.info(f"KmerData initialized with k={k}. Found {len(self.all_unique_potential_kmers)} unique potential k-mers for {len(self.aso_ids)} ASOs.")

class CommandRunner:
    """Handles execution of shell commands."""
    @staticmethod
    def run(command: List[str], work_dir: Optional[str] = None) -> Tuple[str, str]:
        """Runs a shell command and returns its stdout and stderr."""
        current_work_dir = work_dir or os.getcwd()
        logging.debug(f"Running command: {' '.join(command)} in dir: {current_work_dir}")
        try:
            process = subprocess.run(
                command,
                cwd=current_work_dir,
                capture_output=True,
                text=True,
                check=True
            )
            return process.stdout.strip(), process.stderr.strip()
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed with exit code {e.returncode}: {' '.join(e.cmd)}")
            if e.stdout:
                logging.error(f"Stdout: {e.stdout.strip()}")
            if e.stderr:
                logging.error(f"Stderr: {e.stderr.strip()}")
            raise RuntimeError(f"Command failed: {' '.join(e.cmd)} RC: {e.returncode} Stderr: {e.stderr.strip() if e.stderr else 'N/A'}") from e
        except FileNotFoundError:
            logging.error(f"Executable not found for command: {command[0]}")
            raise

class KMC:
    """Wrapper for KMC (K-mer Counter) executable."""
    def __init__(self, kmc_path: str, default_threads: int = 1, default_memory_gb: int = 4):
        if not shutil.which(kmc_path):
            raise FileNotFoundError(f"KMC executable not found at: {kmc_path}")
        self.kmc_path = kmc_path
        self.default_threads = default_threads
        self.default_memory_gb = default_memory_gb
        self.runner = CommandRunner()

    def build_database(self,
                       input_path: str,
                       db_prefix_path: str,
                       k_value: int,
                       temp_dir_path: str,
                       input_type: str = "fm", 
                       min_count: int = 1,
                       threads: Optional[int] = None,
                       memory_gb: Optional[int] = None,
                       additional_args: Optional[List[str]] = None) -> None:
        """
        Builds a KMC database.
        Args:
            input_path: Path to input file (FASTA/FASTQ).
            db_prefix_path: Prefix for the output KMC database files.
            k_value: K-mer length.
            temp_dir_path: Path to the temporary directory for KMC.
            input_type: input in FASTA format (-fa), FASTQ format (-fq), multi FASTA (-fm) or BAM (-fbam) or KMC (-fkmc); default: multi FASTA.
            min_count: Minimum count for k-mers to be stored (-ci).
            threads: Number of threads to use. Defaults to self.default_threads.
            memory_gb: Memory in GB to use. Defaults to self.default_memory_gb.
            additional_args: List of any other KMC arguments. (Should be in kmc arguments format, e.g. ['-ci10'])
        """
        resolved_threads = threads if threads is not None else self.default_threads
        resolved_memory_gb = memory_gb if memory_gb is not None else self.default_memory_gb

        cmd = [
            self.kmc_path,
            f"-k{k_value}",
            f"-{input_type}", # -fm (auto), -fa (fasta), -fq (fastq)
            f"-t{resolved_threads}",
            f"-m{resolved_memory_gb}",
            f"-ci{min_count}",
        ]
        if additional_args:
            cmd.extend(additional_args)
        cmd.extend([input_path, db_prefix_path, temp_dir_path])
        
        logging.info(f"Building KMC database from {input_path} with k={k_value}...")
        self.runner.run(cmd, work_dir=temp_dir_path)
        logging.info(f"KMC database built at {db_prefix_path}")

class KMCTools:
    """Wrapper for KMC_Tools executable."""
    def __init__(self, kmc_tools_path: str):
        if not shutil.which(kmc_tools_path):
            raise FileNotFoundError(f"KMC_Tools executable not found at: {kmc_tools_path}")
        self.kmc_tools_path = kmc_tools_path
        self.runner = CommandRunner()

    def simple_operation(self,
                         left_db_prefix: str,
                         right_db_prefix: str,
                         operation: str, # e.g., 'intersect', 'union', 'kmers_subtract'
                         output_db_prefix: str,
                         work_dir: Optional[str] = None,
                         additional_args: Optional[List[str]] = None) -> None:
        """
        Performs a simple operation (like intersect, union) between two KMC databases.
        Args:
            left_db_prefix: Prefix of the first (left) input KMC database.
            right_db_prefix: Prefix of the second (right) input KMC database.
            operation: The operation to perform (e.g., 'intersect').
            output_db_prefix: Prefix for the output KMC database.
            work_dir: Working directory for the command.
            additional_args: List of any other kmc_tools simple arguments.
        """
        cmd = [self.kmc_tools_path, 'simple', left_db_prefix, right_db_prefix, operation, output_db_prefix]
        if additional_args:
            cmd.extend(additional_args)
        self.runner.run(cmd, work_dir=work_dir)

    def transform_database(self,
                           input_db_prefix: str,
                           operation: str, # e.g., 'dump'
                           output_path: Optional[str] = None, # Required for some operations like dump
                           work_dir: Optional[str] = None,
                           additional_args: Optional[List[str]] = None) -> None:
        """
        Transforms a KMC database (e.g., dumps k-mers and counts).
        Args:
            input_db_prefix: Prefix of the input KMC database.
            operation: The transformation operation (e.g., 'dump').
            output_path: Path for the output file (e.g., for dump operation).
            work_dir: Working directory for the command.
            additional_args: List of any other kmc_tools transform arguments.
        """
        cmd = [self.kmc_tools_path, 'transform', input_db_prefix, operation]
        if output_path: # Some operations like dump require an output path directly in command
             cmd.extend(['-s', output_path]) # Assume -s for sorted dump, adjust if needed
        if additional_args:
            cmd.extend(additional_args)
        
        self.runner.run(cmd, work_dir=work_dir)

class KmerDbQuerier:
    """Encapsulates logic for querying a KMC database for counts of specific k-mers."""
    def __init__(self, kmc_executor: KMC, kmc_tools_executor: KMCTools, k: int):
        self.kmc_executor = kmc_executor
        self.kmc_tools_executor = kmc_tools_executor
        self.k = k

    def get_counts(self, main_kmc_db_prefix: str, target_kmers: Set[str], temp_dir: str) -> Dict[str, int]:
        """
        Queries a KMC database for counts of specific target k-mers.
        Args:
            main_kmc_db_prefix: Prefix of the KMC database to query.
            target_kmers: A set of k-mer sequences to query.
            temp_dir: Temporary directory to write intermediate files.
        Returns:
            A dictionary mapping k-mer sequences to their counts.
        """
        if not target_kmers:
            return {}

        target_kmers_fasta_file = os.path.join(temp_dir, "query_target_kmers.fasta")
        with open(target_kmers_fasta_file, 'w') as f:
            for i, kmer in enumerate(target_kmers):
                f.write(f">query_kmer_{i}\n{kmer}\n")

        target_kmers_db_prefix = os.path.join(temp_dir, "query_target_kmers_db")
        self.kmc_executor.build_database(
            input_path=target_kmers_fasta_file,
            db_prefix_path=target_kmers_db_prefix,
            k_value=self.k,
            temp_dir_path=temp_dir,
            input_type="fm",
            min_count=1,
            threads=1, 
            memory_gb=2 
        )

        intersect_db_prefix = os.path.join(temp_dir, "query_intersect_db")
        self.kmc_tools_executor.simple_operation(
            left_db_prefix=main_kmc_db_prefix,
            right_db_prefix=target_kmers_db_prefix,
            operation='intersect',
            output_db_prefix=intersect_db_prefix,
            work_dir=temp_dir,
            additional_args=["-ocleft"]
        )

        output_counts_file = os.path.join(temp_dir, "query_kmer_counts.txt")
        self.kmc_tools_executor.transform_database(
            input_db_prefix=intersect_db_prefix,
            operation='dump',
            output_path=output_counts_file,
            work_dir=temp_dir
        )
        
        kmer_counts: Dict[str, int] = {kmer: 0 for kmer in target_kmers}
        if os.path.exists(output_counts_file):
            with open(output_counts_file, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) == 2:
                        kmer, count_str = parts[0], parts[1]
                        if kmer in kmer_counts:
                           kmer_counts[kmer] = int(count_str)
        return kmer_counts

class KmerCounter:
    """
    Orchestrates k-mer counting workflows using KMC tools, KmerData, and KmerDbQuerier.
    """

    def __init__(self,
                 pre_mrna_fasta_path: str,
                 potential_secondary_sites: Dict[str, List[Tuple[str, float]]],
                 k: int,
                 kmc_path: str = "kmc",
                 kmc_tools_path: str = "kmc_tools",
                 kmc_db_threads: int = 4,
                 kmc_db_memory_gb: int = 12,
                 gene_processing_workers: int = os.cpu_count() or 1,
                 temp_dir_base: Optional[str] = None,
                 kmc_min_count: int = 1,
                 total_genes_for_matrix: Optional[int] = None):
        """
        Initializes the KmerCounter.

        Args:
            pre_mrna_fasta_path: Path to the gzipped pre-mRNA FASTA file.
            potential_secondary_sites: Dict mapping ASO/target IDs to a list of
                                       tuples (kmer_sequence, ddg_value).
            k: The k-mer length.
            kmc_path: Path to the KMC executable.
            kmc_tools_path: Path to the KMC_Tools executable.
            kmc_db_threads: Number of threads for KMC database creation.
            kmc_db_memory_gb: Memory (in GB) for KMC database creation.
            gene_processing_workers: Number of worker threads for per-gene processing.
            temp_dir_base: Optional base directory for temporary files. If None,
                           the system's default temporary directory is used.
            kmc_min_count: Minimum count for k-mers to be stored by KMC (-ci).
            total_genes_for_matrix: Optional total number of genes for matrix calculation progress tracking.
        """
        if not os.path.exists(pre_mrna_fasta_path):
            raise FileNotFoundError(f"Pre-mRNA FASTA file not found: {pre_mrna_fasta_path}")
        
        self.kmc_executor = KMC(kmc_path, default_threads=kmc_db_threads, default_memory_gb=kmc_db_memory_gb)
        self.kmc_tools_executor = KMCTools(kmc_tools_path)
        
        self.kmer_data = KmerData(potential_secondary_sites, k)
        self.db_querier = KmerDbQuerier(self.kmc_executor, self.kmc_tools_executor, k)

        self.pre_mrna_fasta_path = pre_mrna_fasta_path
        self.k = k
        self.kmc_db_threads = kmc_db_threads
        self.kmc_db_memory_gb = kmc_db_memory_gb
        self.gene_processing_workers = gene_processing_workers
        self.temp_dir_base = temp_dir_base
        self.kmc_min_count = kmc_min_count
        self.total_genes_for_matrix = total_genes_for_matrix

        logging.info(f"KmerCounter initialized with k={k}. Using KmerData with {len(self.kmer_data.all_unique_potential_kmers)} unique potential k-mers for {len(self.kmer_data.aso_ids)} ASOs.")

    def calculate_aggregate_counts(self) -> Dict[str, int]:
        """
        Calculates aggregate off-target counts for each ASO across the entire pre-mRNA FASTA file.
        Returns:
            A dictionary mapping ASO IDs to their total k-mer counts.
        """
        logging.info("Starting aggregate k-mer counting for all ASOs.")
        if not self.kmer_data.all_unique_potential_kmers:
            logging.warning("No unique k-mers to count (from KmerData). Returning empty results.")
            return {aso_id: 0 for aso_id in self.kmer_data.aso_ids}

        main_temp_dir = tempfile.mkdtemp(dir=self.temp_dir_base)
        logging.info(f"Created temporary directory for aggregate counts: {main_temp_dir}")

        try:
            kmc_db_prefix = os.path.join(main_temp_dir, "global_premrna_kmc_db")
            
            logging.info(f"Building global KMC database from {self.pre_mrna_fasta_path}...")
            self.kmc_executor.build_database(
                input_path=self.pre_mrna_fasta_path,
                db_prefix_path=kmc_db_prefix,
                k_value=self.k,
                temp_dir_path=main_temp_dir,
                input_type="fm", 
                min_count=self.kmc_min_count,
                threads=self.kmc_db_threads, 
                memory_gb=self.kmc_db_memory_gb 
            )
            logging.info(f"Global KMC database built at {kmc_db_prefix}")

            logging.info(f"Querying KMC database for {len(self.kmer_data.all_unique_potential_kmers)} unique k-mers...")
            all_kmer_counts = self.db_querier.get_counts(kmc_db_prefix, self.kmer_data.all_unique_potential_kmers, main_temp_dir)
            logging.info("Finished querying KMC database.")

            aso_aggregate_counts: Dict[str, int] = {}
            for aso_id, kmer_list in self.kmer_data.potential_kmers_by_aso.items():
                total_count = 0
                for kmer_seq in kmer_list:
                    total_count += all_kmer_counts.get(kmer_seq, 0)
                aso_aggregate_counts[aso_id] = total_count
            
            logging.info(f"Aggregate counts calculated for {len(aso_aggregate_counts)} ASOs.")
            return aso_aggregate_counts

        finally:
            logging.info(f"Cleaning up temporary directory: {main_temp_dir}")
            shutil.rmtree(main_temp_dir, ignore_errors=True)

    def _parse_fasta(self, file_path: str) -> Iterable[Tuple[str, str, str]]:
        """
        Parses a gzipped FASTA file, yielding (raw_header, gene_id, sequence).
        Gene ID is extracted assuming format like >ENSG...|... or >ENSG...
        """
        open_func = gzip.open if file_path.endswith(".gz") else open
        try:
            with open_func(file_path, 'rt') as f:
                header = None
                sequence_parts = []
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if header and sequence_parts:
                            gene_id_part = header.split('|')[0].lstrip('>')
                            yield header, gene_id_part, "".join(sequence_parts)
                        header = line
                        sequence_parts = []
                    else:
                        if header: #Only collect sequence if a header has been seen
                            sequence_parts.append(line)
                if header and sequence_parts: # Yield the last sequence
                    gene_id_part = header.split('|')[0].lstrip('>')
                    yield header, gene_id_part, "".join(sequence_parts)
        except Exception as e:
            logging.error(f"Error parsing FASTA file {file_path}: {e}")
            raise

    def _process_gene_for_matrix(self, gene_id: str, gene_sequence: str, gene_header: str) -> Tuple[str, Dict[str, int]]:
        """
        Processes a single gene: builds KMC DB, queries k-mers, and returns counts per ASO for this gene.
        """
        thread_id = threading.get_ident() # For unique temp dir names if needed, though gene_id should be unique
        gene_temp_dir = tempfile.mkdtemp(dir=self.temp_dir_base, prefix=f"kmercounter_gene_{gene_id.replace(':','_')}_{thread_id}_")
        # logging.info(f"[Gene: {gene_id}] Processing in temp dir: {gene_temp_dir}")

        try:
            temp_gene_fasta = os.path.join(gene_temp_dir, f"{gene_id.replace(':','_')}.fa")
            with open(temp_gene_fasta, 'w') as f:
                f.write(f"{gene_header}\n{gene_sequence}\n")

            gene_kmc_db_prefix = os.path.join(gene_temp_dir, f"{gene_id.replace(':','_')}_kmc_db")
            
            # KMC for single gene: less threads/memory by default, customize if needed
            # Using -t1 as parallelism is managed by ThreadPoolExecutor
            self.kmc_executor.build_database(
                input_path=temp_gene_fasta,
                db_prefix_path=gene_kmc_db_prefix,
                k_value=self.k,
                temp_dir_path=gene_temp_dir,
                input_type="fa", # Single gene FASTA
                min_count=self.kmc_min_count,
                threads=1, # Single thread as parallelized by ThreadPoolExecutor
                memory_gb=2  # Modest memory for single gene KMC
            )

            gene_kmer_counts = self.db_querier.get_counts(gene_kmc_db_prefix, self.kmer_data.all_unique_potential_kmers, gene_temp_dir)

            aso_counts_for_gene: Dict[str, int] = {}
            for aso_id, kmer_list in self.kmer_data.potential_kmers_by_aso.items():
                total_count = 0
                for kmer_seq in kmer_list:
                    total_count += gene_kmer_counts.get(kmer_seq, 0)
                aso_counts_for_gene[aso_id] = total_count
            
            # logging.info(f"[Gene: {gene_id}] Processed. Found hits for {sum(1 for x in aso_counts_for_gene.values() if x > 0)} ASOs.")
            return gene_id, aso_counts_for_gene

        except Exception as e:
            logging.error(f"[Gene: {gene_id}] Error during processing: {e}")
            # Return empty counts for this gene in case of error to not break the whole matrix
            return gene_id, {aso_id: 0 for aso_id in self.kmer_data.aso_ids}
        finally:
            logging.debug(f"[Gene: {gene_id}] Cleaning up temp dir: {gene_temp_dir}")
            shutil.rmtree(gene_temp_dir, ignore_errors=True)

    def calculate_per_gene_counts_matrix(self) -> Optional[Any]: # polars.DataFrame
        """
        Calculates per-gene off-target counts, producing an ASO x Gene matrix.
        Returns:
            A polars.DataFrame with ASOs as rows, Genes as columns, and k-mer counts as values.
            Returns None if polars is not available.
        """
        if pl is None:
            logging.error("Polars library is not installed. Cannot create per-gene counts matrix.")
            return None
        
        logging.info("Starting per-gene k-mer counting for ASO x Gene matrix.")
        if not self.kmer_data.all_unique_potential_kmers:
            logging.warning("No unique k-mers to count (from KmerData). Returning empty matrix.")
            return pl.DataFrame(schema={aso_id: pl.Int64 for aso_id in self.kmer_data.aso_ids}).transpose(include_header=True, header_column_name='ASO_ID', column_names=[])

        if not self.total_genes_for_matrix or self.total_genes_for_matrix == 0:
            logging.warning("total_genes_for_matrix not provided or is zero. Progress tracking for gene matrix calculation will be limited or disabled. Returning empty matrix.")
            return pl.DataFrame(schema={aso_id: pl.Int64 for aso_id in self.kmer_data.aso_ids}).transpose(include_header=True, header_column_name='ASO_ID', column_names=[])

        results_per_gene: Dict[str, Dict[str, int]] = {} 

        progress_tracker = ProgressTracker(total_items=self.total_genes_for_matrix, description="Processing genes for KMC matrix", update_interval=200)
        actual_genes_submitted_count = 0

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.gene_processing_workers) as executor:
            futures = []
            for raw_header, gene_id, gene_sequence in self._parse_fasta(self.pre_mrna_fasta_path):
                if not gene_id or not gene_sequence: 
                    logging.warning(f"Skipping invalid entry from FASTA: Header='{raw_header}'")
                    continue
                futures.append(executor.submit(self._process_gene_for_matrix, gene_id, gene_sequence, raw_header))
                actual_genes_submitted_count += 1
            
            if actual_genes_submitted_count == 0:
                logging.warning("No valid genes found in FASTA to process after parsing. Returning empty matrix.")
                return pl.DataFrame(schema={aso_id: pl.Int64 for aso_id in self.kmer_data.aso_ids}).transpose(include_header=True, header_column_name='ASO_ID', column_names=[])

            if abs(actual_genes_submitted_count - self.total_genes_for_matrix) > 0.05 * self.total_genes_for_matrix: 
                 logging.warning(f"Number of genes submitted for processing ({actual_genes_submitted_count}) "
                                 f"differs from the initially provided total_genes_for_matrix ({self.total_genes_for_matrix}). "
                                 "Progress ETA might be inaccurate.")

            logging.info(f"Submitted {len(futures)} genes for processing. Waiting for completion...")

            for future in concurrent.futures.as_completed(futures):
                try:
                    gene_id_processed, aso_counts_for_gene = future.result()
                    results_per_gene[gene_id_processed] = aso_counts_for_gene
                except Exception as e:
                    logging.error(f"Error processing a gene future: {e}")
                finally:
                    progress_tracker.update(1) 
        
        logging.info(f"Finished parallel gene processing. Successfully obtained results for {len(results_per_gene)} genes.")

        if not results_per_gene:
            logging.warning("No gene processing results obtained after parallel execution. Returning an empty matrix.")
            return pl.DataFrame(schema={aso_id: pl.Int64 for aso_id in self.kmer_data.aso_ids}).transpose(include_header=True, header_column_name='ASO_ID', column_names=[])

        data_for_df: Dict[str, List[Any]] = {'ASO_ID': self.kmer_data.aso_ids}
        all_processed_gene_ids = sorted(list(results_per_gene.keys()))

        for gene_id_col in all_processed_gene_ids:
            data_for_df[gene_id_col] = [results_per_gene[gene_id_col].get(aso_id, 0) for aso_id in self.kmer_data.aso_ids]
        
        df = pl.DataFrame(data_for_df)

        logging.info(f"Successfully created ASO x Gene matrix with shape: {df.shape}")
        return df
