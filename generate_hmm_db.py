#! python3
# generate_hmm_db.py
# This program will generate a comprehensive HMM database consisting of most
# domain models currently within sequence databases. It was developed as an
# assistant for identifying novel domains via the CANDID program 
# (https://github.com/zkstewart/CANDID) but has more generalised uses as well.

# Import external packages
import argparse, os, shutil

# Define functions for later use
# Argument validation
def validate_args(args):
        import os
        # Validate program execution is successful
        program_execution_check(os.path.join(args.hmmer3dir, 'hmmbuild -h'))
        # Validate that SUPERFAMILY and/or CATH files are locatable if specified
        if args.superfamily != None:
                if not os.path.isfile(args.superfamily):
                        print('"' + args.superfamily + '" was not able to be located.')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        if args.cath != None:
                if not os.path.isfile(args.cath):
                        print('"' + args.cath + '" was not able to be located.')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        # Validate that integer arguments are sensible
        if args.threads < 1:
                print('Threads argument must be an integer greater than 0.')
                print('Fix your input and try again.')
                quit()

def program_execution_check(cmd):
        import subprocess
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        if cmderr.decode("utf-8") != '':
                print('Failed to execute program "' + cmd + '". Is this executable in the location specified/discoverable in your PATH, or does the executable even exist? I won\'t be able to run properly if I can\'t execute this program.')
                print('---')
                print('stderr is below for debugging purposes.')
                print(cmderr.decode("utf-8"))
                print('Program closing now.')
                quit()

# HMM generation-related functions
def fasta_dl(args, filename):
        # Function adapted from https://stackoverflow.com/questions/32288113/python-3-how-to-create-a-text-progress-bar-for-downloading-files
        import urllib.request, os, traceback
        print('Downloading CDD ' + filename + ' file...')
        url = args.url
        u = urllib.request.urlopen(url)
        meta = u.info()
        metaInfo = str(meta).split()
        print ("Download size: " + str(round(int(metaInfo[3])/1000000, 1)) + " megabytes")
        fileTotalbytes=int(metaInfo[3])

        data_blocks = []
        total=0
        # Put rest of function in try:except clause to delete the file if the download fails
        try:
                # Open output file
                out_name = os.path.join(args.outdir, filename)
                f = open(out_name, "wb")
                # Download loop
                while True:
                        block = u.read(1024)
                        data_blocks.append(block)
                        total += len(block)
                        hash = ((60*total)//fileTotalbytes)
                        print("[{}{}] {}%".format('#' * hash, ' ' * (60-hash), int(total/fileTotalbytes*100)), end="\r")

                        # Save to file as we go to reduce memory usage
                        if len(data_blocks) == 10000:   # Saves every 10mb
                                f.write(b''.join(data_blocks))
                                data_blocks = []

                        # Break out of dl loop when complete
                        if not len(block):
                                break
                
                # Final save
                f.write(b''.join(data_blocks))
                f.close()
                u.close()
                print('')
        except:
                var = traceback.format_exc()
                print("Unexpected error:")
                print(var)
                f.close()
                u.close()
                os.remove(out_name)
                quit()

def pfam_dl(args, pfam_filename):
        import urllib.request, os, traceback
        print('Downloading PFAM ' + filename + ' file...')
        url = args.pfamurl
        u = urllib.request.urlopen(url)
        meta = u.info()
        #print(str(meta).split())
        metaInfo = str(meta).split()
        #print(len(metaInfo))
        print ("Download size: " + str(round(int(metaInfo[3])/1000000, 1)) + " megabytes")
        fileTotalbytes=int(metaInfo[3])

        data_blocks = []
        total=0
        # Put rest of function in try:except clause to delete the file if the download fails
        try:
                # Open output file
                out_name = os.path.join(args.outdir, filename)
                f = open(out_name, "wb")
                # Download loop
                while True:
                        block = u.read(1024)
                        data_blocks.append(block)
                        total += len(block)
                        hash = ((60*total)//fileTotalbytes)
                        print("[{}{}] {}%".format('#' * hash, ' ' * (60-hash), int(total/fileTotalbytes*100)), end="\r")

                        # Save to file as we go to reduce memory usage
                        if len(data_blocks) == 10000:   # Saves every 10mb
                                f.write(b''.join(data_blocks))
                                data_blocks = []

                        # Break out of dl loop when complete
                        if not len(block):
                                break
                
                # Final save
                f.write(b''.join(data_blocks))
                f.close()
                u.close()
        except:
                var = traceback.format_exc()
                print("Unexpected error:")
                print(var)
                f.close()
                u.close()
                os.remove(out_name)
                quit()
        
def untar(args, filename, extractdir):
        import os, tarfile, gzip
        print('Extracting ' + filename + '...')
        # Get file details
        file = filename
        extract_dir = os.path.join(args.outdir, extractdir)
        # Jump into output directory
        if file.endswith(".tar.gz"):
                tar = tarfile.open(file, "r:gz")
                tar.extractall(extract_dir)
                tar.close()
        elif file.endswith(".tar"):
                tar = tarfile.open(file, "r:")
                tar.extractall(extract_dir)
                tar.close()
        elif file.endswith(".tar.bz2"):
                tar = tarfile.open(file, "r:bz2")
                tar.extractall(extract_dir)
                tar.close()
        elif file.endswith(".gz"):
                with gzip.open(file, "rb") as tar, open(os.path.join(extract_dir, filename[0:-3]), 'wb') as outfile:
                        for line in tar:
                                outfile.write(line)
        else:
                print('I don\'t recognise the file format for ' + filename)
                print('This script assume it is compressed with .tar, .tar.gz, or .tar.bz2 extension')
                print('To proceed further, I recommend you uncompress this file yourself into a directory titled')
                print(extract_dir)
                print('and then re-run this script, as it will automatically detect this directory.')

def hmmbuild(args, filename):
        # Define thread function for later use
        def build_loop(loop_args, startNum, endNum):
                import subprocess, os
                for i in range(startNum, endNum):
                # Loop through all individual msas to produce hmms
                        file_prefix = '.'.join(loop_args[2][i].split('.')[0:-1])
                        # Continue from previous run if cancelled
                        hmmout_name = os.path.join(loop_args[1], file_prefix + '.hmm')
                        if os.path.isfile(hmmout_name):
                                if os.path.getsize(hmmout_name) != 0:
                                        continue
                        cmd = loop_args[3] + ' ' + hmmout_name + ' ' + os.path.join(loop_args[0], loop_args[2][i])
                        run_hmmbuild = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                        hmmout, hmmerr = run_hmmbuild.communicate()
                        # Error handler (inside loop)
                        if hmmerr.decode("utf-8") != '':
                                raise Exception('hmmbuild error text below\n' + str(hmmerr.decode("utf-8")) + '\nMake sure that you define the -hmmer3dir argument if this directory is not in your PATH. The program crashed when processing file ' + file_prefix)
                        
        # Continue defining the thread function-using function
        import os, threading
        print('Building individual HMMs from ' + filename + ' MSAs...')
        # Get file details
        extract_dir = os.path.join(args.outdir, 'cdd_extraction')
        ind_hmm_dir = os.path.join(args.outdir, 'cdd_individual_models')
        extracted_msas = os.listdir(extract_dir)
        cmd_prefix = os.path.join(args.hmmer3dir, 'hmmbuild')
        dirsize = len(extracted_msas)
        # Set up threading requirements
        threads = args.threads
        chunk_size = int(dirsize / threads)
        loop_args = [extract_dir, ind_hmm_dir, extracted_msas, cmd_prefix, dirsize]
        processing_threads = []
        # Begin the loop
        for i in range(threads):
                start = chunk_size * i
                if i+1 != threads:
                        end = chunk_size * (i+1)
                else:
                        end = dirsize
                build = threading.Thread(target=build_loop, args=(loop_args, start, end))
                processing_threads.append(build)
                build.start()
                print('Initiated thread num ' + str(i+1) + ' for individual hmm building...')

        # Wait for all threads to end.
        for process_thread in processing_threads:
                process_thread.join()
        print('Individual HMM building completed')

def filenum_check(args):
        import os
        # Get the directory of the cdd extraction & count how many files are present
        extract_dir = os.path.join(args.outdir, 'cdd_extraction')
        extracted_msas = os.listdir(extract_dir)
        extr_count = len(extracted_msas)
        # Get the directory of the cdd individual models & count how many files are present
        ind_hmm_dir = os.path.join(args.outdir, 'cdd_individual_models')
        hmms = os.listdir(ind_hmm_dir)
        hmm_count = len(hmms)
        # Return error message if the two values do not correspond
        if extr_count != hmm_count:
                return 'Error'
        else:
                return 'Success'
        
def concat_hmms(args):
        import os
        print('Concatenating individual HMMs into a single .hmm...')
        # Get file details
        ind_hmm_dir = os.path.join(args.outdir, 'cdd_individual_models')
        ind_hmms = os.listdir(ind_hmm_dir)
        outhmm_name = os.path.join(args.outdir, 'cdd_db', 'CDD.hmm')
        outhmm_file = open(outhmm_name, 'w')
        # Concatenate files
        for ind_hmm in ind_hmms:
                outhmm_file.write(open(os.path.join(ind_hmm_dir, ind_hmm), 'r').read())
        outhmm_file.close()

def convert(args, extractdir):
        import subprocess, os
        print('Converting SUPERFAMILY database to HMMER 3.1 version...')
        # Get file details
        if args.superfamily.endswith('.gz'):
                supfam_name = args.superfamily[0:-3]                # If, for whatever reason, SUPERFAMILY either is compressed with something other than just gz, or if it has already been decompressed and the user is providing its uncompressed name as an argument, we need to handle its file name
        else:
                supfam_name = args.superfamily                          # I assume that, if the file name doesn't end with .gz, it was already decompressed. This may cause issues if SUPERFAMILY is distributed with another type of compression in the future.
        file_in = os.path.join(args.outdir, extractdir, supfam_name)
        out_name = file_in + '_3.1'
        # Format command
        cmd_prefix = os.path.join(args.hmmer3dir, 'hmmconvert')
        cmd = os.path.join(args.hmmer3dir, 'hmmconvert') + ' ' + file_in + ' > ' + out_name
        # Run
        run_hmmconvert = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        hmmout, hmmerr = run_hmmconvert.communicate()
        # Error handler (inside loop)
        if hmmerr.decode("utf-8") != '':
                raise Exception('hmmconvert error text below\n' + str(hmmerr.decode("utf-8")) + '\nMake sure that you define the -hmmer3dir argument if this directory is not in your PATH')

def concat_additional(args):
        import os
        print('Concatenating additional HMM databases into a single .hmm...')
        # Open CDD hmm for append
        cddhmm_name = os.path.join(args.outdir, 'cdd_db', 'CDD.hmm')
        cddhmm_file = open(cddhmm_name, 'a')
        # Concatenate SUPERFAMILY file if arguments were provided
        if args.superfamily != 'n':
                print('Concatenating SUPERFAMILY...')
                if args.superfamily.endswith('.gz'):
                        supfam_name = args.superfamily[0:-3]
                else:
                        supfam_name = args.superfamily
                if supfam_name.endswith('_3.1'):
                        donothing = 0
                else:
                        supfam_name += '_3.1'
                supfam_hmm_name = os.path.join(args.outdir, 'superfamily_extraction', supfam_name)
                supfam_hmm_file = open(supfam_hmm_name, 'r')
                for line in supfam_hmm_file:
                        cddhmm_file.write(line)
                supfam_hmm_file.close()
        # Concatenate CATH file if arguments were provided
        if args.cath != 'n':
                print('Concatenating CATH...')
                if args.cath.endswith('.gz'):
                        cath_name = args.cath[0:-3]
                else:
                        cath_name = args.cath
                cath_hmm_name = os.path.join(args.outdir, 'cath_extraction', cath_name)
                cath_hmm_file = open(cath_hmm_name, 'r')
                for line in cath_hmm_file:
                        cddhmm_file.write(line)
                cath_hmm_file.close()
        cddhmm_file.close()
                
        # Leave a note saying that the CDD file also includes SUPERFAMILY/CATH
        note_name = os.path.join(args.outdir, 'cdd_db', 'note.txt')
        note_file = open(note_name, 'w')
        if args.superfamily != 'n':
                note_file.write('This .hmm file contains SUPERFAMILY models as well\n')
        if args.superfamily != 'n':
                note_file.write('This .hmm file contains CATH models as well\n')
        note_file.close()

# Import classes from included script folder
from domfind import hmm_dl

cdd_prefixes = ('cd', 'COG', 'KOG', 'LOAD', 'MTH', 'pfam', 'PHA', 'PRK', 'PTZ', 'sd', 'smart', 'TIGR', 'PLN', 'CHL')    # We use this so we can tell if the user has been using the output directories for anything other than the program. It won't be entirely foolproof, but it should prevent anything major.

#### USER INPUT SECTION
usage = """Usage: <output directory> [-options]
----
%(prog)s is designed to easily facilitate the conversion of the CDD into a format
compatible with HMMER 3.1+. %(prog)s will attempt to automatically continue previously
cancelled runs at the nearest completed step. This script is easily hackable to add your
own msa files into the cdd_extraction directory to create .hmm models, or to add your own
.hmm models to the cdd_individual_models directory to add these to the final CDD.hmm file.
"""

# This directory MUST be dedicated solely to this program; if this program encounters an error, certain directories and their contents will be removed

# Required
p = argparse.ArgumentParser(description=usage, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("outdir", type = str, help="""Specify the name of the directory to download and process the
               CDD into a .hmm readable by HMMER3. If you have already downloaded this file, this argument
               should specify the directory it is contained within, noting that this directory should be
               dedicated solely to this program's operations""")
# Optional cmds
p.add_argument("-h", dest="hmmer3dir", type = str, default = '',
                  help="Specify the directory where HMMER3 executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-t",  dest="threads", type = int, default = 1,
                  help="Optionally specify the number of worker threads when performing the individual HMM building step.")
p.add_argument("-s",  dest="superfamily", type = str, default = None,
                  help="""Optionally specify the location of the SUPERFAMILY HMM file if you
                  want to incorporate these models into the final HMM database file.""")
p.add_argument("-c", dest="cath", type = str, default = None,
                  help="""Optionally specify the name of the CATH HMM file if you want to
                   want to incorporate these models into the final HMM database file.""")
p.add_argument("-r", dest="clean", action = "store_true", default = False,
                  help="Provide this tag if you wish for all intermediate files to be removed once the final CDD.hmm file is created.")
# Future proofing
p.add_argument("-u", dest="url", type = str, default = 'ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/fasta.tar.gz',
                  help="""If the CDD fasta.tar.gz url location changes in the future, specify the new url using this option.
                  Otherwise, do not specify this argument.""")

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
                        if not file.startswith(cdd_prefixes) and not file.lower().endswith('.fasta'):
                                print('I think I can detect files in the ' + os.path.join(args.outdir, 'cdd_extraction') + ' directory that should not exist (i.e., ' + file + ')')
                                print('Are you using this directory for anything other than this script? Move this file and any similar ones out of it to resume the program.')
                                quit()
                
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
