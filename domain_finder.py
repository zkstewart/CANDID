#! python3

# Import external packages
import argparse, os, time, platform
# Import classes from included script folder
from domfind import domfind, domclust, benchparse, peakdetect

#### USER INPUT SECTION
usage = """Usage: <fasta file> <output directory> <threads> [-options]
----
%(prog)s will find novel globular domain regions using all-against-all MMSeqs2 search.
%(prog)s can be run by providing command-line arguments, or by providing a text file with these arguments.
After running this program once with command-line arguments, a text file will
automatically be generated which can be used in future uses. This text file will
be used by default henceforth - any deviations to this text file provided in the
command-line will be saved back into this text file. The only values that must be
provided even with a text file input are the name of the fasta file, the output
directory, and the number of threads.
----
%(prog)s will look inside the output directory (assumed to be located within the same
directory as this script file) for the results of previous runs. These files will be
used by default and will result in skipping previously completed sections of this program.
----
Finally, this script requires CD-HIT, MMSeqs2, signalP, seg, COILS, HMMER3, as well as Python 3.x
and Python 2.7 to be installed. The directories of these programs must be specified in this script.
Additionally, a HMM database file of domains is required to exclude known models from discovery.
A program 'dbb_db_download.py' is provided to help generate this, but you can create your own.
Specify the full path to this file.
----
Run this main script using Python 3.x. If running on a Windows system, install Cygwin including the
main packages (e.g., Perl, interpretters). A comprehensive list of Cygwin packages required is difficult to collate.
----
IMPORTANT: Note that any directories referred to in this program should NOT include any spaces in their name
(e.g., 'output directory' must be 'output_directory'), and that HMMER version 3.1 or above is required.
"""

# Required
p = argparse.ArgumentParser(description=usage, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("fasta", type = str, help="Specify the fasta file to use for all-against-all PSI-BLAST")
p.add_argument("outdir", type = str, help="Specify the name of the output directory")
p.add_argument("threads", type = int, help="Specify the number of threads to use for any steps that multi-threading is enabled for")
# Opts 1: Directory locations
p.add_argument("-mmseqs2dir", dest="mmseqs2dir", type = str,
                  help="Specify the directory where the MMseqs2 executable is located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-cdhitdir", dest="cdhitdir", type = str,
                  help="Specify the directory where CD-HIT executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-hmmer3dir", dest="hmmer3dir", type = str,
                  help="Specify the directory where HMMER3 executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-signalpdir", dest="signalpdir", type = str,
                  help="Specify the directory where signalp executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-segdir", dest="segdir", type = str,
                  help="Specify the directory where seg executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-coilsdir", dest="coilsdir", type = str,
                  help="Specify the directory where the pscoils .py file is located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-mafftdir", dest="mafftdir", type = str,
                  help="Specify the directory where the mafft executables on a Unix system or .bat file on a Windows is located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-cygwindir", dest="cygwindir", type = str,
                  help="If running this script on a Windows system, Cygwin is required. Specify the location of the /bin directory here. If running on other systems, or if this is already in your PATH, you can leave this blank.")
p.add_argument("-python2dir", dest="python2dir", type = str,
                  help="Specify the python2.7 directory that contains python.exe. If this version of Python is the default in your PATH, you can leave this blank.")
# Opts 2: CD-HIT parameters
p.add_argument("-cdc", dest="cdc", type = float,
                  help="This command is equivalent to the -c option that will be provided to CD-HIT.")
p.add_argument("-cdn", dest="cdn", type = int,
                  help="This command is equivalent to the -n option that will be provided to CD-HIT.")
p.add_argument("-cdg", dest="cdg", type = int,
                  help="This command is equivalent to the -G option that will be provided to CD-HIT.")
p.add_argument("-cdas", dest="cdas", type = float,
                  help="This command is equivalent to the -aS option that will be provided to CD-HIT.")
p.add_argument("-cdal", dest="cdal", type = float,
                  help="This command is equivalent to the -aL option that will be provided to CD-HIT.")
p.add_argument("-cdm", dest="cdm", type = int,
                  help="This command is equivalent to the -M option that will be provided to CD-HIT.")
# Opts 3: SignalP parameter
p.add_argument("-signalporg", dest="signalporg", type = str, choices = ['euk', 'gram-', 'gram+', 'EUK', 'GRAM-', 'GRAM+'],
                  help="Specify the type of organism for SignalP. Refer to the SignalP manual if unsure what this means.")
# Opts 4: HMMER3 parameters
p.add_argument("-hmmdb", dest="hmmdb", type = str,
                  help="Specify the full path to the hmm database file to use for HMMER3 domain prediction. It is recommended you use the complementary 'hmm_db_download.py' program to generate this.")
p.add_argument("-hmmeval", dest="hmmeval", type = float,
                  help="Specify the e-value cut-off to enforce for removing known domains from sequences. Default recommended == 1")
p.add_argument("-hmmevalnov",dest="hmmevalnov", type = float,
                  help="Specify the e-value cut-off to enforce for clustering novel domains from sequences. This should be equal to or stricter than that enforced for removing known domains.")
p.add_argument("-skip",dest="skip", type = str, choices = ['cath', 'superfamily', 'both', 'noskip'],
                  help="Optional ability to ignore domain predictions from CATH and/or SUPERFAMILY databases if they are present in your HMM database.")        # Consider whether this should remain in the file versions
# Opts 5: MMseqs2 parameters
p.add_argument("-mms2eval", dest="mms2eval", type = float,
                  help="Specify the e-value cut-off to enforce for returning MMseqs2 hits. Default recommended == 1")
# Opts 5: Alignment free algorithms
p.add_argument("-alf", dest="alf", type = str, choices = ['braycurtis', 'google', 'canberra'],
                  help="Specify the alignment-free algorithm to employ. Recommended to use braycurtis.")
p.add_argument("-reduce", dest="reduce", choices = ['n', 'N', '11', '15'],
                  help="If you wish to supply a reduced protein alphabet to the alignment-free computation step, specify whether this alphabet should be reduced to 15 characters or 11 (11 was used in Alfree benchmark). Recommended not to use ('n').")
# Opts 6: HDBSCAN algorithm parameters
p.add_argument("-minsize", dest="minsize", type = int,
                  help="Dictates the minimum cluster size argument provided to HDBSCAN. Higher numbers will result in identification of domains that occur more frequently in your data. Recommended to use 3 (at least).")
p.add_argument("-minsample", dest="minsample", type = int,
                  help="Dictates the minimum sample size argument provided to HDBSCAN. Higher numbers will increase the strictness with which HDBSCAN clusters sequences, typically resulting in less sensitivity but higher specificity. Recommended to use 2.")
p.add_argument("-leaf", dest="leaf", type = str, choices = ['y', 'n', 'Y', 'N'],
                  help="Changes the HDBSCAN algorithm to use 'leaf' clustering rather than 'excess of mass'. Should result in a higher number of smaller-sized groups being identified. Not using leaf ('n') is recommended.")
p.add_argument("-singleclust", dest="singleclust", type = str, choices = ['y', 'n', 'Y', 'N'],
                  help="Changes the HDBSCAN algorithm to allow the discovery of only a single cluster. By default HDBSCAN recommends that you do not allow this ('n'). If you believe there may be a single novel domain in your data, you can set this to 'y'.")
# Opts 7: Iteration control
p.add_argument("-numiters", dest="numiters", type = int,
                  help="Optionally specify an upper limit on the amount of times this program will iterate through HMMER domain discovery. Recommended value depends on size of dataset. Specifying numiters == 10 means this program will end after 10 iterations OR when no new regions are found, whichever comes first. Numiters == 0 means the program will stop iterating only when no new domains are found.")
# Opts 8: Various parameters (not intended to be changed, but can be overwrote by command-line arguments)
p.add_argument("-cleanAA", dest="cleanAA", type = int, default = 30,
                  help="This value acts as a 'magic number' for many operations; it is based on 30AA being the expected minimum length of a true domain. This value is not intended to be changed; experienced users may wish to do so, however.")
#p.add_argument("-cooccur", dest="cooccur", type = float, default = 0.75,
#                  help="This command is not intended to be changed; experienced users may wish to change this, however.")      ## This value is part of the 'parse_joiner' function. Currently it is not being used; future updates may support a co-occurrence identifier to find potentially incorrectly separated domains (realistically, this would probably be a separate program...)
p.add_argument("-benchmark", dest="benchmark", choices = ['y', 'n', 'Y', 'N'], default = 'n',
                  help="This setting is used specifically for testing. DELETE BEFORE PROGRAM IS SHARED.")
p.add_argument("-rejects", dest="rejects", type = str,
                  help="This setting is used specifically for testing. DELETE BEFORE PROGRAM IS SHARED.")
cmdargs = p.parse_args()

# Create output directory if needed
outputDir = os.path.basename(cmdargs.outdir)
if not os.path.isdir(os.path.join(os.getcwd(), outputDir)):
        os.mkdir(os.path.join(os.getcwd(), outputDir))

# Get the basename of the fasta file
fasta_base = '.'.join(os.path.basename(cmdargs.fasta).split('.')[0:-1])
raw_fasta = cmdargs.fasta                                                       # We just use this value once when running CD-HIT, since we use '.fasta' for all relevant output files but the user input may have a different suffix. It also allows us to easily read in fasta files from directories not in the current working dir.

# Merge text file and command-line arguments if applicable
args = vars(cmdargs)
param_name = os.path.join(os.getcwd(), outputDir, outputDir + '_parameters.txt')
changes = 'n'
if os.path.isfile(param_name):
        with open(param_name, 'r') as txtparams:
                for line in txtparams:
                        line = line.rstrip('\n').rstrip('\r')
                        if line.startswith('#') or line == ' ' or line == '':
                                continue
                        sl = line.split('->')
                        argname = sl[0]
                        argvalue = sl[1].lstrip(' ').rstrip(' ')
                        # Check if this argument is used by the program
                        try:
                                args[argname]
                        except KeyError:                                                                        # This is OK, it just means that there is a line in the parameters file that doesn't correspond to an argument in this program
                                continue                              
                        # Compare text to command-line
                        if argvalue != '' and args[argname] != None and argvalue != str(args[argname]):         # Need to compare them as strings since the text-file doesn't alter the value type like argparse does for some arguments
                                changes = 'y'
                        elif argvalue != '' and args[argname] == None:
                                args[argname] = argvalue

# Validate user inputs to make sure things are sensible
## None values
close = 'n'
for key, value in args.items():
        # Special handling of cygwin argument
        if platform.system() != 'Windows' and value == None:
                continue
        if value == None:
                print(key + ' is missing an argument. Edit your parameter file, or provide this on the command-line.')
                close = 'y'
if close == 'y':
        quit()
## Possible new checks

# Create or update the text file of parameters if necessary [TO-DO: FIX ISSUES WITH INTEGRATING TEXT-FILE PARAMETERS AND CMD-LINE PARAMETERS]
if not os.path.isfile(param_name) or changes == 'y':
        param_text = []
        for key, value in args.items():
                param_text.append(key + '->' + str(value))
        param_text = '\n'.join(param_text)
        with open(param_name, 'w') as param_out:
                param_out.write(param_text)

# Check if the necessary programs are installed and can be reached
## MMSEQS2
if not os.path.isfile(os.path.join(args['mmseqs2dir'], 'mmseqs')) and not os.path.isfile(os.path.join(args['mmseqs2dir'], 'mmseqs.exe')):
        print('I cannot find "mmseqs" at the location provided (' + args['mmseqs2dir'] + ')')
        quit()
## CD-HIT
if not os.path.isfile(os.path.join(args['cdhitdir'], 'cd-hit')) and not os.path.isfile(os.path.join(args['cdhitdir'], 'cd-hit.exe')):
        print('I cannot find "cd-hit" at the location provided (' + args['cdhitdir'] + ')')
        quit()
## HMMER3.1
if not os.path.isfile(os.path.join(args['hmmer3dir'], 'hmmpress')) and not os.path.isfile(os.path.join(args['hmmer3dir'], 'hmmpress.exe')):
        print('I cannot find "hmmpress" at the location provided (' + args['hmmer3dir'] + ')')
        quit()
if not os.path.isfile(os.path.join(args['hmmer3dir'], 'hmmsearch')) and not os.path.isfile(os.path.join(args['hmmer3dir'], 'hmmsearch.exe')):
        print('I cannot find "hmmsearch" at the location provided (' + args['hmmer3dir'] + ')')
        quit()
## SignalP4.1
if not os.path.isfile(os.path.join(args['signalpdir'], 'signalp')) and not os.path.isfile(os.path.join(args['signalpdir'], 'signalp.exe')):
        print('I cannot find "signalp" at the location provided (' + args['signalpdir'] + ')')
        quit()
## SEG
if not os.path.isfile(os.path.join(args['segdir'], 'seg')) and not os.path.isfile(os.path.join(args['segdir'], 'seg.exe')):
        print('I cannot find "seg" at the location provided (' + args['segdir'] + ')')
        quit()
## COILS
if not os.path.isfile(os.path.join(args['coilsdir'], 'psCoils.py')):
        print('I cannot find "psCoils.py" at the location provided (' + args['coilsdir'] + ')')
        quit()
## MAFFT
if platform.system() == 'Windows':
        if not os.path.isfile(os.path.join(args['mafftdir'], 'mafft.bat')):
                print('You are using a Windows OS, but I cannot find "mafft.bat" at the location provided (' + args['mafftdir'] + ')')
                quit()
else:
        if not os.path.isfile(os.path.join(args['mafftdir'], 'mafft')) and not os.path.isfile(os.path.join(args['mafftdir'], 'mafft.exe')):
                print('You are using a non-Windows OS, but I cannot find "mafft" at the location provided (' + args['mafftdir'] + ')')
                quit()
## CYGWIN
if platform.system() == 'Windows':
        if not os.path.isfile(os.path.join(args['cygwindir'], 'bash.exe')):
                print('You are using a Windows OS, but I cannot find "bash" at the location provided (' + args['cygwindir'] + ')')
                quit()
## PYTHON2
if not os.path.isfile(os.path.join(args['python2dir'], 'python')) and not os.path.isfile(os.path.join(args['python2dir'], 'python.exe')):
        print('I cannot find "python" at the location provided (' + args['python2dir'] + ')')
        quit()

print('All programs can be found and the provided arguments appear to be sound. If the programs are not installed properly, however, errors will occur.')

#### DATA PREPARATION
print('### PROGRAM START ###')
print(time.ctime())

### RUN CD-HIT
if not os.path.isfile(os.path.join(os.getcwd(), outputDir, fasta_base + '_cdhit.fasta')):
        domfind.runcdhit(args, outputDir, raw_fasta, fasta_base)

### CHUNK CD-HIT FOR THREADING
# Check if this is necessary
chunking = 'n'
if int(args['threads']) > 1:
        chunk_names = []
        for i in range(0, int(args['threads'])):
                chunk_names.append(os.path.join(os.getcwd(), outputDir, fasta_base + '_cdhit_chunk' + str(i+1) + '.fasta'))
        for name in chunk_names:
                if not os.path.isfile(name):
                        chunking = 'y'
                        break
        # Check one step ahead to see if a previous run was performed with more threads
        if os.path.isfile(os.path.join(os.getcwd(), outputDir, fasta_base + '_cdhit_chunk' + str(i+2) + '.fasta')):
                print('There are more "cdhit_chunk#" chunk files in the output directory than should exist given the number of threads provided in your argument')
                print('This probably means you have leftover files from a previous run which used more threads. Delete all of these "cdhit_chunk#" files to re-run with less threads, or alternatively provide a -threads argument equivalent to the number of chunk files in the output directory')
                quit()
        if chunking == 'y':
                domfind.chunk_fasta(args, outputDir, fasta_base + '_cdhit')
                          
### RUN HMMER3
if not os.path.isfile(os.path.join(os.getcwd(), outputDir, fasta_base + '_hmmer.results')):
        domfind.runhmmer3(args, outputDir, fasta_base)
if args['benchmark'] == 'n':
        if not os.path.isfile(os.path.join(os.getcwd(), outputDir, fasta_base + '_hmmerParsed.results')):
                domfind.hmmerparse(args, args['hmmeval'], outputDir, fasta_base)
else:
        if not os.path.isfile(os.path.join(os.getcwd(), outputDir, fasta_base + '_hmmerParsed.results')):
                benchparse.benchparse(args, outputDir, fasta_base)
if not os.path.isfile(os.path.join(os.getcwd(), outputDir, fasta_base + '_domCut.fasta')):
        domfind.hmmercutter(args, outputDir, fasta_base)

### RUN SIGNALP
if not os.path.isfile(os.path.join(os.getcwd(), outputDir, fasta_base + '_signalp.fasta')):
        domfind.runsignalp(args, outputDir, fasta_base)

### RUN SEG AND COILS
if not os.path.isfile(os.path.join(os.getcwd(), outputDir, fasta_base + '_segcoils.fasta')):
        domfind.runsegandcoils(args, outputDir, fasta_base)

# BENCHMARK
if args['benchmark'] == 'y':
        benchparse.reject_novelty(args, outputDir, fasta_base)

### FINAL PREP CLEAN-UP
if not os.path.isfile(os.path.join(os.getcwd(), outputDir, fasta_base + '_clean.fasta')):
        domfind.cleanseqs(args, outputDir, fasta_base)

#### MMSEQS2 OPERATIONS

### MAKE MMSEQS2 DB
if not os.path.isdir(os.path.join(os.getcwd(), outputDir, 'mms2tmp')):
        os.mkdir(os.path.join(os.getcwd(), outputDir, 'mms2tmp'))
#else:
#        shutil.rmtree(os.path.join(os.getcwd(), outputDir, 'mms2tmp'), ignore_errors=False, onerror=None)       # MMseqs2 has an annoying habit of using files in the temporary folder incorrectly if you change; the only way to prevent errors is to make sure the temporary folder is empty.
#        os.mkdir(os.path.join(os.getcwd(), outputDir, 'mms2tmp'))
        
if not os.path.isfile(os.path.join(os.getcwd(), outputDir, fasta_base + '_mmseqs2DB')) and not os.path.isfile(os.path.join(os.getcwd(), outputDir, fasta_base + '_mmseqs2DB.index')):
        domfind.makemms2db(args, outputDir, fasta_base)  

### RUN MMSEQS2
if not os.path.isfile(os.path.join(os.getcwd(), outputDir, fasta_base + '_mmseqs2SEARCH.m8')):
        domfind.runmms2(args, outputDir, fasta_base)     

### PARSE MMSEQS2
if not os.path.isfile(os.path.join(os.getcwd(), outputDir, fasta_base + '_unclustered_domains.fasta')):
        domfind.parsemms2_peaks(args, outputDir, fasta_base)

#### DOMAIN CLUSTERING
rejects_list = []
current_doms = ''
iterate = 0
loopCount = 0                                                   # Loop count functions provides an alternative exit condition to the iteration loop. If a user specifies iterate 1, it will perform the loop twice and then exit (i.e., it will 'iterate' once). If a user specifies 0, the loop will only exit when iterate == 2, which means no new domain regions were found for 2 loops.
if not os.path.isfile(os.path.join(os.getcwd(), outputDir, fasta_base + '_clustered_domains.fasta')):
        while iterate < 2 and loopCount <= int(args['numiters']):             # This loop works by feeding iterate 0 into the hmmer_grow function. If no new region is found, iterate becomes 1. This is then fed back into hmmer_grow, and if there is still no new regions found, iterate becomes 2 and the while loop ends. We do this because, even if no new region is found when iterate == 0, changes to domain regions according to HMMER results still occur, and this may influence the clustering of the next round which may subsequently result in new regions being discovered.
                alf_matrix1, alf_matrix2, idlist = domclust.alfree_matrix(args, outputDir, fasta_base + '_unclustered_domains.fasta')
                group_dict, rejects_list = domclust.cluster_hdb(args, alf_matrix1, alf_matrix2, idlist, rejects_list)
                if len(group_dict) == 0:
                        print('No potential novel domains were found. Program end.')
                        quit()
                domclust.mafft_align(args, outputDir, fasta_base, group_dict)
                domclust.cluster_hmms(args, outputDir, fasta_base)
                domclust.hmmer3_doms(args, outputDir, fasta_base, '_clean.fasta')
                domfind.hmmerparse(args, args['hmmevalnov'], os.path.join(outputDir, 'tmp_alignments'), 'tmp')
                domclust.parse_domtblout(os.path.join(outputDir, 'tmp_alignments', 'tmp_hmmer.results'), args['hmmevalnov'], os.path.join(outputDir, 'tmp_alignments', 'tmp_hmmerParsed.results'))
                iterate = domclust.hmmer_grow(args, outputDir, fasta_base, group_dict, rejects_list, iterate)
                if int(args['numiters']) != 0:       # This prevents the numiters exit condition from activating
                        loopCount += 1
        # Provide informative loop exit text
        if loopCount > int(args['numiters']):
                print('Program finished successfully after iterating ' + str(loopCount-1) + ' time(s).')
        else:
                print('Program finished successfully after no new domain regions were able to be found.')
print(time.ctime())
### 

