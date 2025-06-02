import os
import gzip
import shlex
import subprocess
import logging
from typing import Optional, List, Tuple, Dict, Any, Union
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from src.utils.genome_utils import Genome, TargetSite, Site
import time
import gget
import urllib.request, urllib.parse
from src.utils.time_utils import timed, format_duration


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



def create_job_config_summary(job_dir: str, config: dict) -> None:
    """
    Create a summary file containing all configuration parameters used for the job.
    
    Parameters:
        job_dir (str): Path to the job directory
        config (dict): Dictionary containing all configuration parameters
    """
    config_file = os.path.join(job_dir, 'job_config.txt')
    
    with open(config_file, 'w') as f:
        f.write("Job Configuration Summary\n")
        f.write("=======================\n\n")
        f.write(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        for key, value in config.items():
            f.write(f"{key}: {value}\n")
            
    logging.info(f"Job configuration summary written to {config_file}")


