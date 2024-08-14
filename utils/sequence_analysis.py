import logging
import multiprocessing
from tqdm import tqdm
import pandas as pd
import configparser
import os
import subprocess
import shlex


# Create a configparser object
config = configparser.ConfigParser()

# Read the configuration file
config.read('config.ini')

def get_chromosomal_positions_per_transcript(transcript, position_in_transcript, ensembl_obj, ensembl_obj_scaffolds = None):
    transcript_id = transcript.split(".")[0]

    try:
        transcript = ensembl_obj.transcript_by_id(transcript_id=transcript_id)
    except Exception as e:
        try:
            transcript = ensembl_obj_scaffolds.transcript_by_id(transcript_id=transcript_id)
        except Exception as e:
            logging.warning(e)
            return

    start_pos, end_pos = calculate_chromosomal_positions(
        transcript.exon_intervals,
        position_in_transcript,
        transcript.strand
    )
    
    key = f'{transcript.contig}:{start_pos}-{end_pos}:{transcript.strand}'
    return key


def calculate_chromosomal_positions(exon_intervals, pos, strand):
    
    accumulated = 0
    
    exon_intervals = sorted(exon_intervals, key=lambda x: x[0], reverse=True)
    
    if strand == '-':
        
        for exon in exon_intervals:
            if (accumulated + (exon[1] - exon[0] + 1)) > pos:
                start_pos = (exon[1] - (pos - accumulated) + 1 - 16)
                return start_pos, start_pos + 16
            
            accumulated += (exon[1] - exon[0] + 1)
            
    elif strand == '+':
        
        for exon in reversed(exon_intervals):
            if (accumulated + (exon[1] - exon[0] + 1)) > pos:
                start_pos = (exon[0] + (pos - accumulated) - 1)
                return start_pos, start_pos + 16
            
            accumulated += (exon[1] - exon[0] + 1)




def calculate_occurances(ensembl_obj, candidates, ensembl_obj_scaffolds = None):
    """
    Calculate the number of distinct positions a sequence occurs at.
    """
    distinct_positions = set()


                
    
    def rcheck(row):

        key = get_chromosomal_positions_per_transcript(row[2], row[3], ensembl_obj, ensembl_obj_scaffolds)
        
        if key not in distinct_positions:
            distinct_positions.add(key)

    candidates.apply(lambda row: rcheck(row), axis=1)
    
    return len(distinct_positions)

def worker(seq, ensembl_obj, candidates, ensembl_obj_scaffolds = None):
    """
    Worker function to calculate occurrences for a given sequence.
    """
    return seq, calculate_occurances(ensembl_obj, candidates, ensembl_obj_scaffolds)

def create_occurrence_dict(unique_seqs, ensembl_obj, sam_out, ensembl_obj_scaffolds = None):
    """
    Create a dictionary of sequence occurrences using multiprocessing.
    """
    num_cpus = os.cpu_count()
    with multiprocessing.Pool(processes=num_cpus) as pool:
        pbar = tqdm(total=len(unique_seqs), desc="Processing Sequences", position=0, leave=True, mininterval=10)

        results = []
        async_results = [pool.apply_async(worker, (seq, ensembl_obj, sam_out[sam_out[9] == seq], ensembl_obj_scaffolds), 
                                          callback=lambda _: pbar.update()) for seq in unique_seqs]

        for r in async_results:
            results.append(r.get())

        pbar.close()

    return dict(results)

def get_kmer_occurances(sam_out, ensembl_obj, ensembl_obj_scaffolds = None):
    """
    Main function to get k-mer occurrences and add them to the SAM output. 
    Returns a Dictionary with kmers as key and the number of occurences as value
    """
    unique_seqs = sam_out[9].unique()
    print(f"Unique sequences count: {len(unique_seqs)}")
    occurrence_dict = create_occurrence_dict(unique_seqs, ensembl_obj, sam_out, ensembl_obj_scaffolds)
    
    return occurrence_dict

def getRNAcofoldEnergy(rnaCofoldInFile):
    rcfOutFileName = os.path.splitext(rnaCofoldInFile)[0] + ".rnacofoldout"
    outFile = os.path.splitext(rnaCofoldInFile)[0] + "_cofold_out.csv"

    
    # Run RNAcofold
    logging.info("Running RNAcofold")
    command = f'RNAcofold -p0  --output-format=D --jobs=0 --noPS --noconv {rnaCofoldInFile}'
    logging.info("Command: {}".format(command))
    
    with open(outFile, 'w') as rcfOutFile:
        process = subprocess.Popen(shlex.split(command), stdout=rcfOutFile, stderr=subprocess.PIPE)
        while True:
            output = process.stderr.readline().decode()
            if output == '' and process.poll() is not None:
                break
            if output:
                logging.info(output.strip())
        rc = process.poll()

    return outFile

def get_exon_id(pos_in_transcript, transcript):
    
    
    accumulated = 0
    
    exons = sorted(transcript.exons, key=lambda x: x.start, reverse=True)
    
    if transcript.strand == '-':
        
        for exon in exons:
            
            if (accumulated + (exon.end - exon.start + 1)) > pos_in_transcript:
                return exon.exon_id
            
            accumulated += (exon.end - exon.start + 1)
            
    elif transcript.strand == '+':
        
        for exon in reversed(exons):
            if (accumulated + (exon.end - exon.start + 1)) > pos_in_transcript:
                return exon.exon_id

            
            accumulated += (exon.end - exon.start + 1)

