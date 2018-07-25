#! python3

# Import external packages
import argparse, os, time, platform
# Import classes from included script folder
from domfind import domfind, domclust, benchparse, peakdetect

# Define functions for later use

# Argument validation
def validate_args(args):
        # Ensure that no None values exist
        for key, value in vars(args).items():
                if value == None and key != 'config':   # Config is the only argument that can be None right now
                        print(key + ' argument is not specified. Make sure to specify this on the command-line or ensure it is present within your config file.')
                        quit()
        # Validate file locations
        if not os.path.isfile(args.fasta):
                print('I am unable to locate the input FASTA file (' + args.fasta + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.hmmdb):
                if args.hmmdb == None:
                        print('The hmmdb argument was not specified on the command-line or in a config file (if config was specified on command-line or discoverable within outdir)')
                        print('Make sure you provide the location to this file and try again.')
                        quit()
                else:
                        print('I am unable to locate the input HMM database file (' + args.hmmdb + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        # Validate file location specifications
        ## MMSEQS2
        if args.mmseqs2dir != '':
                if not os.path.isfile(os.path.join(args.mmseqs2dir, 'mmseqs')) and not os.path.isfile(os.path.join(args.mmseqs2dir, 'mmseqs.exe')):
                        print('I cannot find "mmseqs" or "mmseqs.exe" at the location provided (' + args.mmseqs2dir + ')')
                        quit()
        ## CD-HIT
        if args.cdhitdir != '':
                if not os.path.isfile(os.path.join(args.cdhitdir, 'cd-hit')) and not os.path.isfile(os.path.join(args.cdhitdir, 'cd-hit.exe')):
                        print('I cannot find "cd-hit" or "cd-hit.exe" at the location provided (' + args.cdhitdir + ')')
                        quit()
        ## HMMER3.1
        if args.hmmer3dir != '':
                if not os.path.isfile(os.path.join(args.hmmer3dir, 'hmmpress')) and not os.path.isfile(os.path.join(args.hmmer3dir, 'hmmpress.exe')):
                        print('I cannot find "hmmpress" or "hmmpress.exe" at the location provided (' + args.hmmer3dir + ')')
                        quit()
                if not os.path.isfile(os.path.join(args.hmmer3dir, 'hmmsearch')) and not os.path.isfile(os.path.join(args.hmmer3dir, 'hmmsearch.exe')):
                        print('I cannot find "hmmsearch" or "hmmsearch.exe" at the location provided (' + args.h3dir + ')')
                        quit()
        ## SignalP4.1
        if args.signalpdir != '':
                if not os.path.isfile(os.path.join(args.signalpdir, 'signalp')) and not os.path.isfile(os.path.join(args.signalpdir, 'signalp.exe')):
                        print('I cannot find "signalp" at the location provided (' + args.signalpdir + ')')
                        quit()
        ## SEG
        if args.segdir != '':
                if not os.path.isfile(os.path.join(args.segdir, 'seg')) and not os.path.isfile(os.path.join(args.segdir, 'seg.exe')):
                        print('I cannot find "seg" at the location provided (' + args.segdir + ')')
                        quit()
        ## COILS
        if args.coilsdir != '':
                if not os.path.isfile(os.path.join(args.coilsdir, 'psCoils.py')):
                        print('I cannot find "psCoils.py" at the location provided (' + args.coilsdir + ')')
                        quit()
        ## MAFFT
        if args.mafftdir != '':
                if platform.system() == 'Windows':
                        if not os.path.isfile(os.path.join(args.mafftdir, 'mafft.bat')):
                                print('You are using a Windows OS, but I cannot find "mafft.bat" at the location provided (' + args.mafftdir + ')')
                                quit()
                else:
                        if not os.path.isfile(os.path.join(args.mafftdir, 'mafft')) and not os.path.isfile(os.path.join(args.mafftdir, 'mafft.exe')):
                                print('You are using a non-Windows OS, but I cannot find "mafft" or "mafft.exe" at the location provided (' + args.mafftdir + ')')
                                quit()
        ## CYGWIN
        if args.cygwindir != '':
                if platform.system() == 'Windows':
                        if not os.path.isfile(os.path.join(args.cygwindir, 'bash.exe')):
                                print('You are using a Windows OS, but I cannot find "bash.exe" at the location provided (' + args.cygwindir + ')')
                                quit()
        ## PYTHON2
        if args.python2dir != '':
                if not os.path.isfile(os.path.join(args.python2dir, 'python')) and not os.path.isfile(os.path.join(args.python2dir, 'python.exe')):
                        print('I cannot find "python" or "python.exe" at the location provided (' + args.python2dir + ')')
                        quit()
        # Validate program execution is successful
        program_execution_check(os.path.join(args.mmseqs2dir, 'mmseqs -h'))
        program_execution_check(os.path.join(args.cdhitdir, 'cd-hit -h'))
        program_execution_check(os.path.join(args.hmmer3dir, 'hmmpress -h'))
        program_execution_check(os.path.join(args.hmmer3dir, 'hmmsearch -h'))
        program_execution_check(os.path.join(args.segdir, 'seg'))
        program_execution_check(os.path.join(args.python2dir, 'python -h'))
        py2Exe = python_version_check(args.python2dir)
        program_execution_check(py2Exe + ' ' + os.path.join(args.coilsdir, 'psCoils.py -h'))
        
        
        if platform.system() == 'Windows':
                 program_execution_check(os.path.join(args.mafftdir, 'mafft.bat -h'))
        else:
                program_execution_check(os.path.join(args.mafftdir, 'mafft -h'))
        if platform.system() == 'Windows':
                 program_execution_check(os.path.join(args.cygwindir, 'bash.exe'))
        if platform.system() == 'Windows':
                cygwin_program_execution_check(args.outdir, args.cygwindir, args.signalpdir, 'signalp -h')
        else:
                program_execution_check(os.path.join(args.signalpdir, 'signalp -h'))
        
        
        # Validate integer arguments
        intArgs = ['threads', 'cdn', 'cdg', 'cdm', 'minsize', 'minsample', 'numiters', 'cleanAA']
        for entry in intArgs:
                if type(vars(args)[entry]) == str:                                      # If these values aren't strings, they were specified on command-line
                        if "'" in vars(args)[entry] or '"' in vars(args)[entry]:        # Handle a plausible mistake the user may make
                                print('You have quotation marks surrounding your value for the ' + entry + ' argument. Make sure these don\'t exist on your command-line input or in the config file and try again.')
                                quit()
                        if not vars(args)[entry].isdigit():                             # If they are a string, we need them to be able to become an integer; this checks for that and stops program execution if the check fails
                                print(entry + ' argument must be an integer. You specified ' + vars(args)[entry] + ' on command-line or in the config file. Fix this and try again.')
                                quit()
        # Validate float arguments
        floatArgs = ['cdc', 'cdas', 'cdal', 'hmmeval', 'hmmevalnov', 'mms2eval']
        for entry in floatArgs:
                if type(vars(args)[entry]) == str:
                        if "'" in vars(args)[entry] or '"' in vars(args)[entry]:
                                print('You have quotation marks surrounding your value for the ' + entry + ' argument. Make sure these don\'t exist on your command-line input or in the config file and try again.')
                                quit()
                        try:
                                float(vars(args)[entry])
                        except:
                                print(entry + ' argument must be able to become a float value (i.e., a number with decimal places). You specified ' + vars(args)[entry] + ' on command-line or in the config file. Fix this and try again.')
                                quit()
        # Validate arguments with choice specification
        if args.signalporg not in ['euk', 'gram-', 'gram+', 'EUK', 'GRAM-', 'GRAM+']:
                print('signalporg must be a value within the below list. Fix this in your config file and try again.')  # We know this mistake is in the config file since the command-line will enforce this choice
                print(['euk', 'gram-', 'gram+', 'EUK', 'GRAM-', 'GRAM+'])
                quit()
        if args.skip not in ['cath', 'superfamily', 'both', 'noskip']:
                print('skip must be a value within the below list. Fix this in your config file and try again.')
                print(['cath', 'superfamily', 'both', 'noskip'])
                quit()
        if args.alf not in ['braycurtis', 'google', 'canberra']:
                print('alf must be a value within the below list. Fix this in your config file and try again.')
                print(['braycurtis', 'google', 'canberra'])
                quit()
        if args.reduce not in ['n', '11', '15']:
                print('reduce must be a value within the below list. Fix this in your config file and try again.')
                print(['n', '11', '15'])
                quit()

def program_execution_check(cmd):
        import subprocess
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        if cmderr.decode("utf-8") != '' and not cmderr.decode("utf-8").startswith('Usage'):     # Need this extra check for seg since it puts its usage information into stderr rather than stdout
                print('Failed to execute program "' + cmd + '". Is this executable in the location specified, or does the executable exist? I won\'t be able to run properly if I can\'t execute this program.')
                print('Program closing now.')
                raise Exception

def cygwin_program_execution_check(outDir, cygwinDir, exeDir, exeFile):
        import subprocess
        # Format script for cygwin execution
        scriptText = '"' + os.path.join(exeDir, exeFile) + '" -h'
        scriptFile = file_name_gen('tmpscript', '.sh')
        with open(os.path.join(outDir, scriptFile), 'w') as fileOut:
                fileOut.write(scriptText)
        # Format cmd for execution
        cmd = os.path.join(cygwinDir, 'bash') + ' -l -c ' + os.path.join(outDir, scriptFile).replace('\\', '/')
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        if cmderr.decode("utf-8") != '':
                print('Failed to execute program "' + cmd + '". Is this executable in the location specified, or does the executable exist? I won\'t be able to run properly if I can\'t execute this program.')
                print('Program closing now.')
                raise Exception

def python_version_check(pythonDir):
        import subprocess, os
        if pythonDir == '':
                cmd = 'python --version'
        else:
                cmd = os.path.join(pythonDir, 'python') + ' --version'
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        cmderr = cmderr.decode("utf-8")
        if 'Python 3.' in cmderr:
                try:
                        if pythonDir == '':
                                program_execution_check('python2')
                                return 'python2'
                        else:
                                program_execution_check(os.path.join(pythonDir, 'python2'))
                                return os.path.join(pythonDir, 'python2')
                except:
                        print('I\'m unable to find a Python 2 executable.')
                        print('You didn\'t provide a path so I looked in your PATH for "python", and this was a Python 3.X version.')
                        print('I then checked for "python2" and didn\'t find it... Long story short, you need to specify an explicit path to a Python 2 version if it isn\'t discoverable in your path.')
                        print('Program closing now.')
                        quit()
        elif 'Python 2.' not in cmderr:
                try:
                        if pythonDir == '':
                                program_execution_check('python2')
                                return 'python2'
                        else:
                                program_execution_check(os.path.join(pythonDir, 'python2'))
                                return os.path.join(pythonDir, 'python2')
                except:
                        print('I\'m unable to find a Python 2 executable.')
                        print('You didn\'t provide a path so I looked in your PATH for "python", and this was a version that did not appear to be a 3.X or 2.X version.')
                        print('I then checked for "python2" and didn\'t find it... Long story short, you need to specify an explicit path to a Python 2 version if it isn\'t discoverable in your path.')
                        print('Program closing now.')
                        quit()
        else:
                if pythonDir == '':
                        return 'python'
                else:
                        return os.path.join(pythonDir, 'python')

def output_arg_handling(args):
        # Handle output directory
        if not os.path.isdir(args.outdir):
                print('You have not specified an existing outdir, so that means we are starting a new run.')
                os.mkdir(args.outdir)
        else:
                # Config file
                print('You have specified an existing outdir, so that means we are resuming a run.')
                if args.config == None:
                        # If we're resuming a run and we haven't specified a config file, find the config file present in the outdir; a standardised naming scheme is used, so we should be able to find it
                        configFiles = []
                        for file in os.listdir(args.outdir):
                                if os.path.basename(args.outdir) + '_run' in file and '.config' in file:
                                        configFiles.append(file)
                        if configFiles == [] and os.listdir(args.outdir) != []: # If the outdir is empty, then the program probably died during validation or the user made the directory before running this program
                                print('You didn\'t specify a config file on command-line, and I couldn\'t find a file which matches the expected naming scheme within this directory.')
                                print('A file which resembles "' + os.path.basename(args.outdir) + '_run*.config should be in this directory.')
                                print('Since I couldn\'t find such a file, I\'m going to exit now. To fix this problem, make sure this file is in the directory, or explicitly refer to a config file on the command-line.')
                                quit()
                        if os.listdir(args.outdir) != []:
                                configFiles.sort(key = lambda x: -int(x.rsplit('_run', maxsplit=1)[1].split('.config')[0]))
                                args.config = os.path.join(args.outdir, configFiles[0])
                                print('You didn\'t specify a config file on command-line, so I\'m going to use the most recent config file present in this directory.')
                                print('This looks like "' + args.config + '"... does this seem right to you? If it isn\'t, specify the config file explicitly on the command-line.')
                else:
                        # If we're resuming a run and we HAVE specified a config file, warn the user that any changes could have unpredictable results
                        print('You specified a config file on the command-line which I was able to find. However, note that if this config file differs to the one used in the original program run, unexpected results may occur.')
                        print('I\'m going to assume you know what you\'re doing (even though you should just allow me to find the config file within the outdir). If weird errors occur, remember this message when you\'re scratching your head.')
        return args

def default_parameter_dict(inputKey):
        '''There are a handful of parameters not in here since they are "mandatory" and should
        have been provided either on the command-line or discoverable within the outdir. These
        include the fasta and outdir (mandatory command-line arguments), the hmmdb (which needs
        to be specified on command-line or in the config file) as well as the config file itself
        which will have been specified on the command-line, discovered by the output_arg_handling
        function, or the user specified all arguments on the command-line (which we'll check shortly)
        '''
        defaultParams = {'threads': 1, 'mmseqs2dir': '', 'cdhitdir': '',
                         'hmmer3dir': '', 'signalpdir': '', 'segdir': '',
                         'mafftdir': '', 'cygwindir': '', 'cdc': 0.4,
                         'cdn': 2, 'cdg': 0, 'cdas': 0.9, 'cdal': 0.6, 'cdm': 1000,
                         'signalporg': 'euk', 'hmmeval': 1e-1, 'hmmevalnov': 1e-1,
                         'skip': 'noskip', 'mms2eval': 1e-1, 'alf': 'google', 'reduce': 'n',
                         'minsize': 3,'minsample': 2, 'leaf': False, 'singleclust': False,
                         'numiters': 0, 'cleanAA': 30,
                         'fasta': None, 'outdir': None, 'hmmdb': None, 'config': None, 'coilsdir': None, 'python2dir': None}      # This lets us know that we shouldn't be specifying defaults for these arguments
        if inputKey in defaultParams:
                return defaultParams[inputKey]
        else:
                return None

def config_file_args_update(args, configFile):
        # Ensure that configFile is discoverable
        if not os.path.isfile(configFile):
                print('I am unable to locate the input config file (' + configFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Parse config file and update args where applicable
        with open(configFile, 'r') as fileIn:
                for line in fileIn:
                        line = line.rstrip('\r\n')
                        # Skip empty lines; they shouldn't be here, but no harm done (yet)
                        if line.startswith('#') or line == ' ' or line  == '':
                                continue
                        sl = line.split('=')
                        argname = sl[0].rstrip(' ')
                        argvalue = sl[1].lstrip(' ')
                        # Check if this argument is used by the program
                        if default_parameter_dict(argname) == None:
                                print('Your config file is formatted incorrectly. An argument "' + argname + '" is present in this file but not recognised by this program.')
                                print('In the interest of ensuring this program works correctly, you need to make sure this config file is accurate and only contains parameters for the operation of this program.')
                                print('I\'m going to exit now and leave the rest up to you.')
                                quit()
                        # Make sure this argument belongs in a config file
                        elif default_parameter_dict(argname) == False:
                                print('Your config file is formatted incorrectly. An argument "' + argname + '" is present in this file which should either be specified explicitly on the command-line, or contains a reference to a config file within the config file itself.')
                                print('In the interest of ensuring this program works correctly, you need to make sure this config file is accurate and does not refer to mandatory arguments that do not belong in the config file.')
                                print('I\'m going to exit now and leave the rest up to you.')
                                quit()
                        # Update our args value if necessary
                        if vars(args)[argname] == None: # By doing this, we allow command-line arguments to overrule config file arguments
                                vars(args)[argname] = argvalue
        return args

def defaults_args_update(args):
        for key, value in vars(args).items():
                if value == None:
                        default = default_parameter_dict(key)
                        if type(default) != None:                     # The other 'defaults' that are booleans are to be ignored
                                vars(args)[key] = default
        return args

## General purpose arguments
def file_name_gen(prefix, suffix):
        ongoingCount = 2
        while True:
                if not os.path.isfile(prefix + '1' + suffix):
                        return prefix + '1' + suffix
                elif os.path.isfile(prefix + str(ongoingCount) + suffix):
                        ongoingCount += 1
                else:
                        return prefix + str(ongoingCount) + suffix

#### USER INPUT SECTION
usage = """Usage: <fasta file> <output directory> <threads> [-options]
----
%(prog)s will find novel globular domain regions using all-against-all MMSeqs2 search.
%(prog)s can be run by providing command-line arguments, by providing a config file with 
these arguments, or with a combination of the two. Arguments provided on the command-line
will overrule any within your config file, so beware of this behaviour!
----
Outputs from this program will be stored within the specified outdir. Within this
outdir will be a config file for your run; re-running this program and specifying the same
outdir will attempt to resume the program run. This program to look at the files in this
folder and determine where it exited, and it will also make use of the config file in this directory.
Thus, you do not need to specify a config file when resuming a previous run. It is important
to note that if you choose to resume a run and specify command-line or config file arguments
which differ to those used in the initial program run, unexpected results may occur.
----
This program requires CD-HIT, MMseqs2, signalP, seg, COILS, HMMER3, as well as Python 3.x
and Python 2.7 to be installed. The directories of these programs must be specified in this script
or discoverable within your system's PATH variable.
----
Additionally, a HMM database file of domains is required to exclude known models from discovery.
A program 'hmm_db_download.py' is provided to help generate this, but you can create your own.
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
p.add_argument("fasta", type = str, help="Specify the fasta file from which domains will be identified.")
p.add_argument("outdir", type = str, help="Specify the name of the output directory; this will be created if it doesn't exist, and if it does, we will attempt to resume any previous runs based on the files within (you do not need to specify any of the below arguments in this case).")
# Opts 0: Basic program requirements which do not need to be specified in all cases
p.add_argument("-config", dest="config", type = str,
                  help="Specify the config file used for parameter settings; if not provided and you are beginning a new run, defaults will be used; if you are resuming a previous run, this program will read in the parameter file within the specified 'outdir' directory.")
p.add_argument("-threads", dest="threads", type = int,
                  help="Specify the number of threads to use for any steps that are multi-thread capable.") # I think we technically use processes and not threads, but let's not get into pettiness
# Opts 1: Directory locations
p.add_argument("-mms2dir", dest="mmseqs2dir", type = str,
                  help="Specify the directory where the MMseqs2 executable is located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-cdhdir", dest="cdhitdir", type = str,
                  help="Specify the directory where CD-HIT executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-h3dir", dest="hmmer3dir", type = str,
                  help="Specify the directory where HMMER3 executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-sigpdir", dest="signalpdir", type = str,
                  help="Specify the directory where signalp executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-segdir", dest="segdir", type = str,
                  help="Specify the directory where seg executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-mafftdir", dest="mafftdir", type = str,
                  help="Specify the directory where the mafft executables on a Unix system or .bat file on a Windows is located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-cwdir", dest="cygwindir", type = str,
                  help="If running this script on a Windows system, Cygwin is required. Specify the location of the /bin directory here. If running on other systems, or if this is already in your PATH, you can leave this blank.")
p.add_argument("-py2dir", dest="python2dir", type = str,
                  help="Specify the python2.7 directory that contains python.exe. This MUST be specified here or in your config file.")
p.add_argument("-coilsdir", dest="coilsdir", type = str,
                  help="Specify the directory where the pscoils .py file is located. This MUST be specified here or in your config file.")
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
p.add_argument("-sigporg", dest="signalporg", type = str, choices = ['euk', 'gram-', 'gram+', 'EUK', 'GRAM-', 'GRAM+'],
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
p.add_argument("-reduce", dest="reduce", choices = ['n', '11', '15'],
                  help="If you wish to supply a reduced protein alphabet to the alignment-free computation step, specify whether this alphabet should be reduced to 15 characters or 11 (11 was used in Alfree benchmark). Recommended not to use ('n').")
# Opts 6: HDBSCAN algorithm parameters
p.add_argument("-minsize", dest="minsize", type = int,
                  help="Dictates the minimum cluster size argument provided to HDBSCAN. Higher numbers will result in identification of domains that occur more frequently in your data. Recommended to use 3 (at least).")
p.add_argument("-minsample", dest="minsample", type = int,
                  help="Dictates the minimum sample size argument provided to HDBSCAN. Higher numbers will increase the strictness with which HDBSCAN clusters sequences, typically resulting in less sensitivity but higher specificity. Recommended to use 2.")
p.add_argument("-leaf", dest="leaf", action = "store_true", default = False,
                  help="Changes the HDBSCAN algorithm to use 'leaf' clustering rather than 'excess of mass'. Should result in a higher number of smaller-sized groups being identified. Not using leaf is recommended.")
p.add_argument("-singleclust", dest="singleclust",  action = "store_true", default = False,
                  help="Changes the HDBSCAN algorithm to allow the discovery of only a single cluster. By default HDBSCAN recommends that you do not allow this. If you believe there may only be a single novel domain in your data, you can provide this argument.")
# Opts 7: Iteration control
p.add_argument("-numiters", dest="numiters", type = int,
                  help="Optionally specify an upper limit on the amount of times this program will iterate through HMMER domain discovery. Recommended value depends on size of dataset. Specifying numiters == 10 means this program will end after 10 iterations OR when no new regions are found, whichever comes first. Numiters == 0 means the program will stop iterating only when no new domains are found.")
# Opts 8: Various parameters (not intended to be changed, but can be overwrote by command-line arguments)
p.add_argument("-cleanAA", dest="cleanAA", type = int, default = 30,
                  help="This value acts as a 'magic number' for many operations; it is based on 30AA being the expected minimum length of a true domain. This value is not intended to be changed; experienced users may wish to do so, however.")
#p.add_argument("-cooccur", dest="cooccur", type = float, default = 0.75,
#                  help="This command is not intended to be changed; experienced users may wish to change this, however.")      ## This value is part of the 'parse_joiner' function. Currently it is not being used; future updates may support a co-occurrence identifier to find potentially incorrectly separated domains (realistically, this would probably be a separate program...)
p.add_argument("-benchmark", dest="benchmark", action = "store_true", default = False,
                  help="This setting is used specifically for testing. DELETE BEFORE PROGRAM IS SHARED.")
# Opts 9: Alternative program operations
p.add_argument("-generate_config", dest="generate_config", action = "store_true", default = False,
                  help="Instead of running this program, instead generate a .config file within the specified outdir; this will be a combination of default parameters plus any you specify here on the command-line.")

args = p.parse_args()

# Handle output directory for a new or existing run
args = output_arg_handling(args)

# Combine command-line & config file arguments, filling in defaults for any gaps
if args.config != None:
        args = config_file_args_update(args, args.config)

# Update any blank arguments with defaults
args = defaults_args_update(args)

# Validate our arguments to ensure they are accurate
validate_args(args)

# Generate a config file within the output directory

if args.generate_config:
        print('Since -generate_config was provided, I am now stopping program execution after config file generation within ' + args.outdir)
        quit()

stophere

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

print('All programs can be found and the provided arguments appear to be sound. If the programs are not installed properly, however, errors will occur.')

#### DATA PREPARATION
print('### PROGRAM START ###')
print(time.ctime())


## SET UP THE WORKING DIRECTORY
#file_name_gen(prefix, suffix)

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
        #domfind.parsemms2_peaks(args, outputDir, fasta_base)
        domfind.parsemms2_nccheck(args, outputDir, fasta_base)

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

