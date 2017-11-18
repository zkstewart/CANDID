# Import external packages
import argparse, os, shutil
# Import classes from included script folder
from domfind import hmm_dl

cdd_prefixes = ('cd', 'COG', 'KOG', 'LOAD', 'MTH', 'pfam', 'PHA', 'PRK', 'PTZ', 'sd', 'smart', 'TIGR', 'PLN', 'CHL')    # We use this so we can tell if the user has been using the output directories for anything other than the program. It won't be entirely foolproof, but it should prevent anything major.

#### USER INPUT SECTION
usage = """Usage: <output directory> <molecule type> [-options]
----
%(prog)s will find novel globular domain regions using all-against-all PSI-BLAST.
This script, during the individual HMM building step, will automatically continue previously cancelled runs. Individual HMM building is very time consuming.
Note: This script should be easily hackable to add your own msa files into the cdd_extraction directory to create .hmm models,
or to add your own .hmm models to the cdd_individual_models directory to add these to the final CDD.hmm file.
"""

# Required
p = argparse.ArgumentParser(description=usage)
p.add_argument("outdir", type = str, help="Specify the name of the directory to download and process the CDD into a .hmm readable by HMMER3. This directory MUST be dedicated solely to this program; if this program encounters an error, certain directories and their contents will be removed.")
# Optional cmds
p.add_argument("-hmmer3dir", "--hmmer3dir", dest="hmmer3dir", type = str, default = '',
                  help="Specify the directory where HMMER3 executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-clean", "--clean", dest="clean", choices = ['y', 'n'], default = 'n',
                  help="Specify 'y' if you wish for all intermediate files to be deleted once the final CDD.hmm file is created.")
p.add_argument("-threads", "--threads", dest="threads", type = int, default = 1,
                  help="Optionally specify the number of worker threads when performing the individual HMM building step.")
p.add_argument("-superfamily", "--superfamily", dest="superfamily", type = str, default = 'n',
                  help="Optionally specify the name of the SUPERFAMILY HMM database file if downloaded and placed within the parent output directory for incorporation into larger database file.")
p.add_argument("-cath", "--cath", dest="cath", type = str, default = 'n',
                  help="Optionally specify the name of the CATH HMM database file if downloaded and placed within the parent output directory for incorporation into larger database file.")
# Future proofing
p.add_argument("-url", "--url", dest="url", type = str, default = 'ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/fasta.tar.gz',
                  help="If the CDD fasta.tar.gz file changes location in the future, specify its new location using this option. Otherwise, leave this blank.")
#p.add_argument("-pfamurl", "--pfamurl", dest="pfamurl", type = str, default = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.fasta.gz',
#                  help="If the Pfam-A.fasta.gz file changes location in the future, specify its new location using this option. Otherwise, leave this blank.")

args = p.parse_args()
cdd_filename = os.path.basename(args.url)
#if args.superfamily != 'n':
#        superfamily_filename = os.path.basename(args.superfamily)       

# Create output directory if needed
if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

#### CORE PROCESSES

### CHECK IF RUNNING THIS IS NECESSARY
if os.path.isfile(os.path.join(args.outdir, 'CDD.hmm')):
        print('The CDD.hmm file is already present in this directory. Specify a new output directory, or delete this file if it is outdated.')
        quit()

### DOWNLOAD CDD
# Download
if not os.path.isfile(os.path.join(args.outdir, cdd_filename)):
        hmm_dl.fasta_dl(args, cdd_filename)

# Untar
if not os.path.isdir(os.path.join(args.outdir, 'cdd_extraction')):
        os.mkdir(os.path.join(args.outdir, 'cdd_extraction'))
        try:
                hmm_dl.untar(args, os.path.join(args.outdir, cdd_filename), 'cdd_extraction')
        except Exception, e:
                #shutil.rmtree(os.path.join(args.outdir, 'cdd_extraction'))      # This could be dangerous if this directory were being used for something other than this program...
                print('Untarring file ' + os.path.join(args.outdir, cdd_filename) + ' failed. Check the error log below')
                print(str(e))
                # Delete the temporary file(s) intelligently
                tmpdir_contents = os.listdir(os.path.join(args.outdir, 'cdd_extraction'))
                safe = 'y'
                for file in tmpdir_contents:
                        if not file.startswith(cdd_prefixes) and not file.endswith('.FASTA'):
                                print('I think I can detect files in the ' + os.path.join(args.outdir, 'cdd_extraction') + ' directory that should not exist (i.e., ' + file + ')')
                                print('
                        if file != os.path.join(args.outdir, cdd_filename)
                
if args.superfamily != 'n' and not os.path.isdir(os.path.join(args.outdir, 'superfamily_extraction')):
        os.mkdir(os.path.join(args.outdir, 'superfamily_extraction'))
        hmm_dl.untar(args, args.superfamily, 'superfamily_extraction')

if args.cath != 'n' and not os.path.isdir(os.path.join(args.outdir, 'cath_extraction')):
        os.mkdir(os.path.join(args.outdir, 'cath_extraction'))
        hmm_dl.untar(args, args.cath, 'cath_extraction')

### BUILD HMM MODEL
# Individual HMM construction from CDD MSAs
if not os.path.isdir(os.path.join(args.outdir, 'cdd_individual_models')):
        os.mkdir(os.path.join(args.outdir, 'cdd_individual_models'))

#if not os.path.isdir(os.path.join(args.outdir, 'cdd_db')):
#        hmm_dl.hmmbuild(args, cdd_filename)
hmm_dl.hmmbuild(args, cdd_filename)     # We can afford to spend a bit of time validating that the files are not 0kb (which is all we do if this step has been previously completed)

# Check that all HMMs are present
buildresult = hmm_dl.filenum_check(args)
if buildresult == 'Error':
        print('Not all MSAs present in the ' + os.path.join(args.outdir, 'cdd_extraction') + ' directory appear to have built HMMs successfully. Recommend that you re-run the program to fix this and identify the problem file(s).')
        quit()
else:
        print('Validated that HMM building appears to have worked.')

# Concatenate individual CDD HMMs
#if not os.path.isfile(os.path.join(args.outdir, 'CDD.hmm')):
if not os.path.isdir(os.path.join(args.outdir, 'cdd_db')):
        os.mkdir(os.path.join(args.outdir, 'cdd_db'))
        hmm_dl.concat_hmms(args)

# Concatenate any additional databases specified (SUPERFAMILY, CATH) into the .hmm file (oh boy is this going to be large!)
# Convert SUPERFAMILY to 3.1
if args.superfamily != 'n':
        if args.superfamily.endswith('_3.1'):
                donothing = 0
        elif args.superfamily.endswith('.gz'):
                if not os.path.isfile(os.path.join(args.outdir, 'superfamily_extraction', args.superfamily[0:-3] + '_3.1')):
                        hmm_dl.convert(args, 'superfamily_extraction')
        else:
                if not os.path.isfile(os.path.join(args.outdir, 'superfamily_extraction', args.superfamily + '_3.1')):
                        hmm_dl.convert(args, 'superfamily_extraction')

# Combine files
if args.superfamily != '\n' or args.cath != '\n':
        hmm_dl.concat_additional(args)


# All done!
print('Finished formatting a .hmm file representing the CDD (and SUPERFAMILY/CATH if specified). Running the domain_finder.py script will automatically \'hmmpress\' this file for further use.')
