import subprocess
import shlex
import configparser 
import logging


# Create a configparser object
config = configparser.ConfigParser()

# Read the configuration file
config.read('config.ini')

def build_bowtie_index(e_release, g_assembly, species, bowtie_index):
            
    logging.info("Running Bowtie2 index build")

    if species == 'human':
        command = f'bowtie2-build {config["DEFAULT"]["DataDir"]}/pyensembl/GRCh{g_assembly}/ensembl{e_release}/Homo_sapiens.GRCh{g_assembly}.cdna.all.fa.gz {config["DEFAULT"]["DataDir"]}/bowtie2Home/{bowtie_index}'
    elif species == 'mouse':
        command = f'bowtie2-build {config["DEFAULT"]["DataDir"]}/pyensembl/GRCm{g_assembly}/ensembl{e_release}/Homo_sapiens.GRCm{g_assembly}.cdna.all.fa.gz {config["DEFAULT"]["DataDir"]}/bowtie2Home/{bowtie_index}'
    
    logging.info("Command: {}".format(command))

        
    return_code = subprocess.call(shlex.split(command))
    logging.info("Return Code: {}".format(return_code))

    return return_code
