#! python3
# generate_hmm_db.py
# This program will generate a comprehensive HMM database consisting of most
# domain models currently within sequence databases. It was developed as an
# assistant for identifying novel domains via the CANDID program 
# (https://github.com/zkstewart/CANDID) but has more generalised uses as well.

# Import external packages
import argparse, os

# Define functions for later use
# Argument validation
def validate_args(args):
        import os
        # Get the full path of outdir argument
        args.outdir = os.path.abspath(args.outdir)
        # Validate program execution is successful
        program_execution_check(os.path.join(args.hmmer3dir, 'hmmbuild -h'))
        # Validate that SUPERFAMILY and/or CATH files are locatable if specified
        if args.superfamily != None:
                if not os.path.isfile(args.superfamily):
                        print('"' + args.superfamily + '" was not able to be located.')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
                args.superfamily = os.path.abspath(args.superfamily)    # Use the full path to make later options easier
        if args.cath != None:
                if not os.path.isfile(args.cath):
                        print('"' + args.cath + '" was not able to be located.')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
                args.cath = os.path.abspath(args.cath)
        # Validate that integer arguments are sensible
        if args.threads < 1:
                print('Threads argument must be an integer greater than 0.')
                print('Fix your input and try again.')
                quit()
        return args

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
def file_dl(url, outdir, filename):
        'Function adapted from https://stackoverflow.com/questions/32288113/python-3-how-to-create-a-text-progress-bar-for-downloading-files'
        # Set up
        import urllib.request, os, traceback
        print('# Downloading "' + filename + '"...')
        # Extract file information
        u = urllib.request.urlopen(url)
        meta = u.info()
        metaInfo = str(meta).split()
        fileTotalbytes=int(metaInfo[3])
        print ("# Download size: " + str(round(int(metaInfo[3])/1000000, 1)) + " megabytes")
        # Main download loop [try:except clause is used to delete the file if the download fails]
        data_blocks = []
        total=0
        try:
                # Open output file
                out_name = os.path.join(outdir, filename)
                f = open(out_name, "wb")
                # Download loop
                while True:
                        block = u.read(1024)
                        data_blocks.append(block)
                        total += len(block)
                        hash = ((60*total)//fileTotalbytes)
                        print("# [{}{}] {}%".format('#' * hash, ' ' * (60-hash), int(total/fileTotalbytes*100)), end="\r")

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
                print("# Unexpected error:")
                print(var)
                f.close()
                u.close()
                os.remove(out_name)
                quit()

def untar(outdir, file, skipUntar):
        # Set up
        import os, shutil, tarfile, gzip
        # Main function
        if file.endswith(".tar.gz"):
                if skipUntar == False:
                        tar = tarfile.open(file, "r:gz")
                        tar.extractall(outdir)
                        tar.close()
                return os.path.join(outdir, os.path.basename(file[:-7]))
        elif file.endswith(".tar"):
                if skipUntar == False:
                        tar = tarfile.open(file, "r:")
                        tar.extractall(outdir)
                        tar.close()
                return os.path.join(outdir, os.path.basename(file[:-4]))
        elif file.endswith(".tar.bz2"):
                if skipUntar == False:
                        tar = tarfile.open(file, "r:bz2")
                        tar.extractall(outdir)
                        tar.close()
                return os.path.join(outdir, os.path.basename(file[:-8]))
        elif file.endswith(".gz"):
                if skipUntar == False:
                        with gzip.open(file, "rb") as tar, open(os.path.join(outdir, os.path.basename(file[0:-3])), 'wb') as fileOut:
                                for line in tar:
                                        fileOut.write(line)
                return os.path.join(outdir, os.path.basename(file[0:-3]))
        else:
                if skipUntar == False:
                        print('## "' + file + '" does not have an extension I recognise for decompression; I assume this has already occurred.')
                        if not os.path.isfile(os.path.join(outdir, os.path.basename(file))):
                                print('## I will make a copy of this file to the specified outdir "' + outdir + '".')
                                shutil.copy(file, outdir)
                        else:
                                print('## A file with name "' + os.path.basename(file) + '" exists at the specified outdir "' + outdir + '".')
                                print('## I won\'t overwrite this file; if this isn\'t the most up-to-date version of this file, you\'ll need to delete it and re-run this program.')
                return os.path.join(outdir, os.path.basename(file))

def cluster_hmms(msaDir, hmmer3dir):
        # Set up
        import os, subprocess
        # Build HMMs from MSA directory
        msaFastas = os.listdir(msaDir)
        hmms = []
        for msa in msaFastas:
                outputFileName = msa.rsplit('.', maxsplit=1)[0] + '.hmm'
                hmms.append(outputFileName)
                # Format cmd
                cmd = os.path.join(hmmer3dir, 'hmmbuild') + ' "' + os.path.join(msaDir, outputFileName) + '" "' + os.path.join(msaDir, msa) + '"'
                # Run hmmbuild
                run_hmmbuild = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                hmmout, hmmerr = run_hmmbuild.communicate()
                if hmmerr.decode("utf-8") != '':
                        raise Exception('hmmbuild error text below\n' + str(hmmerr.decode("utf-8")) + '\nMake sure that you define the -h3dir argument if this directory is not in your PATH')

def hmmbuild_dir(msaDir, hmmOutDir, hmmer3dir, threads):
        # Set up
        import os, threading, math, subprocess
        # Define functions integral to this one
        def run_hmmbuild_dir(fileNameList, msaDir, outputDir, start, end):
                for msa in fileNameList[start:end]:
                        outputFileName = os.path.join(outputDir, msa.rsplit('.', maxsplit=1)[0] + '.hmm')
                        # Enable resumption of a previously failed run
                        if os.path.isfile(outputFileName):
                                if os.path.getsize(outputFileName) != 0:
                                        continue
                        # Format cmd
                        cmd = os.path.join(hmmer3dir, 'hmmbuild') + ' "' + os.path.join(outputFileName) + '" "' + os.path.join(msaDir, msa) + '"'
                        # Run hmmbuild
                        run_hmmbuild = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                        hmmout, hmmerr = run_hmmbuild.communicate()
                        if hmmerr.decode("utf-8") != '':
                                raise Exception('hmmbuild error text below\n' + str(hmmerr.decode("utf-8")) + '\nMake sure that you define the hmmer3dir argument if this directory is not in your PATH.\nProgram crashed when processing "' + outputFileName + '".')
        # Get our directory's file contents exclusive of marker file
        msaFileNames = []
        for file in os.listdir(msaDir):
                if not file.endswith('.complete'):
                        msaFileNames.append(file)
        # Validate that file inputs are correct
        for file in msaFileNames:
                if os.path.isdir(file) == True:
                        print('hmmbuild_dir: A directory exists within "' + msaDir + '"; only files should exist here.')
                        print('The culprit is "' + file + '". Move or delete this directory and try again.')
                        quit()
        # Set up threading requirements
        rawNum = len(msaFileNames) / int(threads)               # In cases where threads > len(msaFileNames), rawNum will be less than 1. numRoundedUp will equal the number of threads, and so we'll end up rounding these to 1. Yay!
        numRoundedUp = round((rawNum % 1) * threads, 0)         # By taking the decimal place and multiplying it by the num of threads, we can figure out how many threads need to be rounded up to process every MSA
        breakPoints = [0]                                       # Seed this list with 0 at the start so we can loop through it more easily
        ongoingCount = 0
        for i in range(threads):
                if i+1 <= numRoundedUp:                                         # i.e., if two threads are being rounded up, we'll round up the first two loops of this
                        breakPoints.append(math.ceil(rawNum) + ongoingCount)    # Round up the rawNum, and also add our ongoingCount which corresponds to the number of MSA already put into a chunk
                        ongoingCount += math.ceil(rawNum)
                else:
                        breakPoints.append(math.floor(rawNum) + ongoingCount)
                        ongoingCount += math.floor(rawNum)
                if ongoingCount >= len(msaFileNames):                           # Without this check, if we have more threads than MSA files, we can end up with "extra" numbers in the list (e.g., [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]).
                        break  
        # Begin the loop
        processing_threads = []
        ongoingCount = 0    # This will keep track of what MSA number we are on
        for i in range(len(breakPoints) - 1):                                   # This lets us dictate how many threads we actually run since the user may have specified more threads than MSAs
                start = breakPoints[i]
                end = breakPoints[i+1]
                build = threading.Thread(target=run_hmmbuild_dir, args=(msaFileNames, msaDir, hmmOutDir, start, end))
                processing_threads.append(build)
                build.start()

        # Wait for all threads to end.
        for process_thread in processing_threads:
                process_thread.join()

def filenum_check(directory1, directory2):
        # Set up
        import os
        # Count the number of files in each directory exclusive of '.complete' files [these are our marker files that the process completed successfully and shouldn't count]
        dir1Files = os.listdir(directory1)
        dir1Count = 0
        for i in range(len(dir1Files)):
                if not dir1Files[i].endswith('.complete'):
                        dir1Count += 1
        dir2Files = os.listdir(directory2)
        dir2Count = 0
        for i in range(len(dir2Files)):
                if not dir2Files[i].endswith('.complete'):
                        dir2Count += 1
        # Return True if the two values were identical, else False
        if dir1Count == dir2Count:
                return True
        else:
                return False

def concat_hmms(hmmFilesDirList, concatOutDir, concatFileName):
        # Set up
        import os
        # Ensure hmmFilesDirList is properly formatted
        if type(hmmFilesDirList) == str:
                hmmFilesDirList = [hmmFilesDirList]
        elif type(hmmFilesDirList) != list:
                print('concat_hmms: hmmFilesDirList type is not recognisable. It should be a list, but instead it is ' + str(type(hmmFilesDirList)) + '.')
                print('Fix your input and try again.')
                quit()
        if hmmFilesDirList == []:
                print('concat_hmms: hmmFilesDirList is empty. I don\'t know what to do in this situation since it shouldn\'t happen.')
                print('Code your call to this function properly to skip it.')
                quit()
        # Get file details from each directory
        hmmFiles = []
        for value in hmmFilesDirList:
                if os.path.isdir(value) == True:
                        hmmDirFiles = os.listdir(value)
                        for hmm in hmmDirFiles:
                                hmmFiles.append(os.path.join(value, hmm))
                else:
                        hmmFiles.append(value)
        concatHmmName = os.path.join(concatOutDir, concatFileName)
        # Concatenate HMMs
        with open(concatHmmName, 'w') as fileOut:
                for hmm in hmmFiles:
                        fileOut.write(open(hmm, 'r').read())

def convert_hmm_db(hmmDbFileIn, hmmDbFileOut):
        # Set up
        import subprocess, os
        # Format command
        cmd = os.path.join(args.hmmer3dir, 'hmmconvert') + ' ' + hmmDbFileIn + ' > ' + hmmDbFileOut
        # Run
        run_hmmconvert = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        hmmout, hmmerr = run_hmmconvert.communicate()
        # Error handler (inside loop)
        if hmmerr.decode("utf-8") != '':
                raise Exception('hmmconvert error text below\n' + str(hmmerr.decode("utf-8")) + '\nMake sure that you define the hmmer3dir argument if this directory is not in your PATH')

def hmmpress(hmmDbFile, hmmer3dir):
        # Set up
        import subprocess
        # Format command
        cmd = os.path.join(hmmer3dir, 'hmmpress') + ' -f "' + hmmDbFile + '"'
        # Run
        run_hmmpress = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        hmmout, hmmerr = run_hmmpress.communicate()
        if hmmerr.decode("utf-8") != '':
                raise Exception('hmmpress error text below' + str(hmmerr.decode("utf-8")))

def create_blank_file(fileName):
        fileOut = open(fileName, 'w')
        fileOut.close()

# Define tuple which may be useful later to tell if the user has been using the output directories for anything other than the program; it won't be entirely foolproof, but it should prevent major mistakes from occurring
cddPrefixes = ('cd', 'COG', 'KOG', 'LOAD', 'MTH', 'pfam', 'PHA', 'PRK', 'PTZ', 'sd', 'smart', 'TIGR', 'PLN', 'CHL')

#### USER INPUT SECTION
usage = """
%(prog)s is designed to easily facilitate the conversion of the CDD and 
SUPERFAMILY/CATH into a format compatible with HMMER 3.1+. %(prog)s will
attempt to automatically continue previously cancelled runs at the nearest completed
step. This script is easily hackable to add or remove msa files in the cdd_extraction
directory to alter the output .hmm file.
"""

# This directory MUST be dedicated solely to this program; if this program encounters an error, certain directories and their contents will be removed

# Required
p = argparse.ArgumentParser(description=usage, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("outdir", type = str, help="""Specify the name of the directory to download and process the
               CDD into a .hmm readable by HMMER3. If you have already downloaded this file, this argument
               should specify the directory it is contained within, noting that this directory should be
               dedicated solely to this program's operations""")
# Optional cmds
p.add_argument("-hmm", "-hmmer3dir", dest="hmmer3dir", type = str, default = '',
                  help="Specify the directory where HMMER3 executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-t",  "-threads", dest="threads", type = int, default = 1,
                  help="Optionally specify the number of worker threads when performing the individual HMM building step.")
p.add_argument("-s",  "-superfamily", dest="superfamily", type = str, default = None,
                  help="""Optionally specify the location of the SUPERFAMILY HMM file if you
                  want to incorporate these models into the final HMM database file.""")
p.add_argument("-c", "-cath", dest="cath", type = str, default = None,
                  help="""Optionally specify the name of the CATH HMM file if you want to
                   want to incorporate these models into the final HMM database file.""")
# Future proofing
p.add_argument("-u", "-url", dest="url", type = str, default = 'ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/fasta.tar.gz',
                  help="""If the CDD fasta.tar.gz url location changes in the future, specify the new url using this option.
                  Otherwise, do not specify this argument.""")

args = p.parse_args()
args = validate_args(args)

### PROGRAM SET UP

# Create output directory if needed
if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

# Figure out our file names and directories
cddFileName = os.path.basename(args.url)
cddExtractDir = os.path.join(args.outdir, 'cdd_extraction')
cddHmmDir = os.path.join(args.outdir, 'cdd_individual_models')
spfamExtractDir = os.path.join(args.outdir, 'superfamily_extraction')
cathExtractDir = os.path.join(args.outdir, 'cath_extraction')

### CDD HANDLING

# Download
if not os.path.isfile(os.path.join(args.outdir, cddFileName)):
        print('# Downloading CDD files')
        file_dl(args.url, args.outdir, cddFileName)

# Untar CDD
if not os.path.isdir(cddExtractDir):
        os.mkdir(cddExtractDir)
if not os.path.isfile(os.path.join(cddExtractDir, 'cdd_untar.complete')):               # This should let us know if untarring was previously performed successfully while still being "hackable" by removing unwanted models from the directory
        try:
                print('# Decompressing CDD files')
                untar(cddExtractDir, os.path.join(args.outdir, cddFileName), False)     # We don't care about returning the updated file name from the untar function
                create_blank_file(os.path.join(cddExtractDir, 'cdd_untar.complete'))    # This relates to ^^
        except Exception as e:
                print('Decompressing file "' + cddFileName + '" to "' + cddExtractDir + ' failed. Check the error log below')
                print(str(e))
                # Delete the file(s) intelligently
                tmpdirContents = os.listdir(cddExtractDir)
                for file in tmpdirContents:
                        if not file.startswith(cddPrefixes) and not file.lower().endswith('.fasta') and file != 'cdd_untar.complete':
                                print('I think I can detect files in "' + cddExtractDir + '" that should not exist (i.e., ' + file + ')')
                                print('Are you using this directory for anything other than this script? Move or delete all files in this directory to resume the program.')
                                quit()

### OPTIONAL DB HANDLING

# SUPERFAMILY
if args.superfamily != None:
        if not os.path.isdir(spfamExtractDir):
                os.mkdir(spfamExtractDir)
        if not os.path.isfile(os.path.join(spfamExtractDir, 'superfamily_untar.complete')):
                try:
                        print('# Decompressing SUPERFAMILY files')
                        args.superfamily = untar(spfamExtractDir, args.superfamily, False)              # This gives us the new location of the decompressed file
                        create_blank_file(os.path.join(spfamExtractDir, 'superfamily_untar.complete'))
                except Exception as e:
                        print('Decompressing file "' + args.superfamily + '" to "' + spfamExtractDir + ' failed. Check the error log below')
                        print(str(e))
                        # Delete the file(s) intelligently
                        tmpdirContents = os.listdir(spfamExtractDir)
                        for file in tmpdirContents:
                                if file != 'hmmlib' and file != 'superfamily_untar.complete':           # SUPERFAMILY's file is currently called 'hmmlib'; the only other file that could be here (but shouldn't if we are handling an exception) is the .complete file
                                        print('I think I can detect files in "' + spfamExtractDir + '" that should not exist (i.e., ' + file + ')')
                                        print('Are you using this directory for anything other than this script? Move or delete all files in this directory to resume the program.')
                                        quit()
        else:
                args.superfamily = untar(spfamExtractDir, args.superfamily, True)                       # True lets us skip the decompression step while still obtaining our file name in the output directory minus compression extension
                                                                                                        # This might be relevant when the user resumes a run that failed without updating their input arguments
# CATH
if args.cath != None:
        if not os.path.isdir(cathExtractDir):
                os.mkdir(cathExtractDir)
        if not os.path.isfile(os.path.join(cathExtractDir, 'cath_untar.complete')):
                try:
                        print('# Decompressing CATH files')
                        args.cath = untar(cathExtractDir, args.cath, False)                             # This gives us the new location of the decompressed file
                        create_blank_file(os.path.join(cathExtractDir, 'cath_untar.complete'))
                except Exception as e:
                        print('Decompressing file "' + args.cath + '" to "' + cathExtractDir + ' failed. Check the error log below')
                        print(str(e))
                        # Delete the file(s) intelligently
                        tmpdirContents = os.listdir(cathExtractDir)
                        for file in tmpdirContents:
                                if not file.endswith('.lib') and file != 'cath_untar.complete':         # CATH's file is currently called 'jackhmmer.S##.hmm3.lib'
                                        print('I think I can detect files in "' + cathExtractDir + '" that should not exist (i.e., ' + file + ')')
                                        print('Are you using this directory for anything other than this script? Move or delete all files in this directory to resume the program.')
                                        quit()
        else:
                args.cath = untar(cathExtractDir, args.cath, True)

### BUILD CDD HMM DB

# Create directory and build our CDD HMM models
if not os.path.isdir(cddHmmDir):
        os.mkdir(cddHmmDir)

print('# Converting CDD fasta to HMM & validating')
hmmbuild_dir(cddExtractDir, cddHmmDir, args.hmmer3dir, args.threads)            # We can afford to spend a bit of time validating that the files are not 0kb (which is all we do if this step has been previously completed)

# Check that all HMMs are present
success = filenum_check(cddExtractDir, cddHmmDir)
if success == False:
        print('Not all MSAs present in the "' + cddExtractDir + '" directory appear to have been converted into HMMs successfully. Recommend that you re-run the program to automatically fix this and identify the problem file(s).')
        quit()
print('# Validated that HMM building appears to have worked.')

# Concatenate individual CDD HMMs
print('# Concatenating individual CDD HMM models into a single large HMM file')
dbFileName = os.path.join(args.outdir, 'CDD.hmm')                               # If we don't add SUPERFAMILY/CATH to this file, this will be our final database name
if not os.path.isfile(os.path.join(args.outdir, 'cdd_hmm.complete')) or not os.path.isfile(os.path.join(args.outdir, 'CDD.hmm')):
        concat_hmms(cddHmmDir, args.outdir, 'CDD.hmm')
        create_blank_file(os.path.join(args.outdir, 'cdd_hmm.complete'))
else:
        print('## A "CDD.hmm" file is already present in outdir "' + args.outdir + '". I won\'t overwrite this file.')
        print('## Specify a new outdir, move the file elsewhere, or delete this file if it is outdated.')

### CONCATENATE ADDITIONAL DBs

# Convert SUPERFAMILY to 3.1 if relevant
if args.superfamily != None:
        if not args.superfamily.endswith('_3.1'):                               # If this already has _3.1 suffix, I'll assume someone already converted it using this program and are re-running this program using this file as an input
                if not os.path.isfile(os.path.join(spfamExtractDir, 'supfam_conv.complete')) or not os.path.isfile(args.superfamily + '_3.1'):  # args.superfamily already has the full path leading to it
                        print('# Making SUPERFAMILY HMM file compatible with HMMER 3.1+')
                        convert_hmm_db(args.superfamily, args.superfamily + '_3.1')
                        create_blank_file(os.path.join(spfamExtractDir, 'supfam_conv.complete'))
                        args.superfamily = args.superfamily + '_3.1'
                else:
                        args.superfamily = args.superfamily + '_3.1'            # If we are resuming a run which failed in a step below and the user doesn't update their input, we need to make sure we use the converted file

# Combine files
if args.superfamily != None or args.cath != None:
        # Format our function inputs and concatenated HMM name
        additionalHmms = []
        concatFileName = 'CDD'
        if args.superfamily != None:
                additionalHmms.append(args.superfamily)
                concatFileName += '_SUPFAM'
        if args.cath != None:
                additionalHmms.append(args.cath)
                concatFileName += '_CATH'
        concatFileName += '.hmm'
        print('# Concatenating additional database(s) into final HMM file inclusive of:')
        if args.superfamily != None:
                print('# > SUPERFAMILY')
        if args.cath != None:
                print('# > CATH')
        # Create concatenated file if relevant
        if not os.path.isfile(os.path.join(args.outdir, 'added_concat.complete')) or not os.path.isfile(os.path.join(args.outdir, concatFileName)):
                concat_hmms(additionalHmms, args.outdir, concatFileName)
                create_blank_file(os.path.join(args.outdir, 'added_concat.complete'))
        else:
                print('## A "' + concatFileName + '" file is already present in outdir "' + args.outdir + '". I won\'t overwrite this file.')
                print('## Specify a new outdir, move the file elsewhere, or delete this file if it is outdated.')
        dbFileName = os.path.join(args.outdir, concatFileName)                  # This will overwrite the database name from above

### PRESS FINAL FILE
if not os.path.isfile(os.path.join(args.outdir, 'press.complete')) or (not os.path.isfile(dbFileName + '.h3f') and not os.path.isfile(dbFileName + '.h3i') and not os.path.isfile(dbFileName + '.h3m') and not os.path.isfile(dbFileName + '.h3p')):
        print('# Running hmmpress to prepare database for use')
        hmmpress(dbFileName, args.hmmer3dir)
        create_blank_file(os.path.join(args.outdir, 'press.complete'))

# All done!
print('# Finished formatting a .hmm file representing the CDD (and SUPERFAMILY/CATH if specified).')
print('# This is called "' + os.path.basename(dbFileName) + '" and is located at "' + os.path.dirname(dbFileName) + '".')
