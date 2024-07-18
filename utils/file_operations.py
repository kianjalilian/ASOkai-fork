import os
from ftplib import FTP

def collect_scaffold(path, genome_assembly, ensembl_release):
    """ Download the specified human scaffold from ENSEMBL if it is not already present in
        current directory. In either case, return the file path. """
        
    filepath = path + f"/pyensembl/GRCh{genome_assembly}/ensembl{ensembl_release}/"
    filename = f"Homo_sapiens.GRCh{genome_assembly}.{ensembl_release}.chr_patch_hapl_scaff.gtf.gz"

    if not os.path.exists(filepath+filename):  # Don't re-download.
        ftp = FTP('ftp.ensembl.org')
        ftp.login()
        ftp.cwd(f'pub/release-{ensembl_release}/gtf/homo_sapiens')
        
        os.makedirs(filepath, exist_ok=True)

        with open(filepath + filename, 'wb') as fp:
            
            try:
                ftp.retrbinary("RETR " + filename, fp.write)
            except:
                print('# Could not collect Scaffold')
                return None
            
            print('# Downloaded', filename, 'Scaffold')
    else:
        print('# Using', filename, 'Scaffold')
    return filepath + filename
