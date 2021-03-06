# Assistant functions
def consecutive_character_coords(inputString, character, base, outType):
        # Parse the index positions of the specified character to derive start-stop coordinates of character stretches
        charIndices = []
        for i in range(len(inputString)):
                if inputString[i] == character:
                        if base == 0:
                                charIndices.append(i)
                        elif base == 1:
                                charIndices.append(i+1)
                        else:
                                print('Base value is wrong. Need to fix.')
                                quit()
        charCoords = []
        for i in range(len(charIndices)):
                if i == 0:
                        charStart = charIndices[i]
                        charStretch = 0         # This acts 0-based, a charStretch of 0 means it's 1 character long
                elif i != len(charIndices) - 1:
                        if charIndices[i] == charIndices[i-1] + 1:
                                charStretch += 1
                        else:
                                # Save
                                if outType == 'coords':
                                        charCoords.append(str(charStart) + '-' + str(charStart + charStretch)) # Note that this does not act like a Python range(), it is everything up to AND including the final index
                                elif outType == 'pairs':
                                        charCoords.append([charStart, charStart + charStretch])
                                # Other stuff
                                charStretch = 0
                                charStart = charIndices[i]
                else:
                        if charIndices[i] == charIndices[i-1] + 1:
                                charStretch += 1
                        # Save
                        if outType == 'coords':
                                charCoords.append(str(charStart) + '-' + str(charStart + charStretch))
                        elif outType == 'pairs':
                                charCoords.append([charStart, charStart + charStretch])
                        # Other stuff
                        charStretch = 0
                        charStart = charIndices[i]
                        if charIndices[i] != charIndices[i-1] + 1:
                                charStretch = 0
                                charStart = charIndices[i]
                                # Save
                                if outType == 'coords':
                                        charCoords.append(str(charStart) + '-' + str(charStart + charStretch))
                                elif outType == 'pairs':
                                        charCoords.append([charStart, charStart + charStretch])
        return charCoords

def thread_file_name_gen(prefix, threadNum):
        import os
        ongoingCount = 0
        while True:
                if not os.path.isfile(prefix + threadNum):
                        return prefix + threadNum
                elif os.path.isfile(prefix + threadNum + '.' + str(ongoingCount)):
                        ongoingCount += 1
                else:
                        return prefix + threadNum + '.' + str(ongoingCount)

def clean_seqs(fastaFile, length, outputFileName):
        # Set up
        from Bio import SeqIO
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        with open(outputFileName, 'w') as fileOut:
                for record in records:
                        sequence = str(record.seq)
                        cleanSeq = sequence.replace('x', '')
                        if len(cleanSeq) < int(length):
                                continue
                        # Output
                        fileOut.write('>' + record.description + '\n' + sequence + '\n')

def chunk_fasta(outputDir, inputFasta, fastaSuffix, threads):
        import os, math
        from Bio import SeqIO
        # Quickly handle 1 thread arguments
        if threads == 1:
                return [inputFasta]
        # Count number of sequences in file
        with open(inputFasta, 'r') as inFile:
                numSeqs = 0
                for line in inFile:
                        if line.startswith('>'):
                                numSeqs += 1
        # Find out where we are chunking the file
        rawNum = numSeqs / threads                              # In cases where threads > numSeqs, rawNum will be less than 1. numRoundedUp will equal the number of threads, and so we'll end up rounding these to 1. Yay!
        numRoundedUp = round((rawNum % 1) * threads, 0)         # By taking the decimal place and multiplying it by the num of threads, we can figure out how many threads need to be rounded up to process every sequence
        chunkPoints = []
        ongoingCount = 0
        for i in range(threads):
                if i+1 <= numRoundedUp:                 # i.e., if two threads are being rounded up, we'll round up the first two loops of this
                        chunkPoints.append(math.ceil(rawNum) + ongoingCount)    # Round up the rawNum, and also add our ongoingCount which corresponds to the number of sequences already put into a chunk
                        ongoingCount += math.ceil(rawNum)
                else:
                        chunkPoints.append(math.floor(rawNum) + ongoingCount)
                        ongoingCount += math.floor(rawNum)
                if ongoingCount >= numSeqs:             # Without this check, if we have more threads than sequences, we can end up with "extra" numbers in the list (e.g., [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]).
                        break                           # This doesn't actually affect program function, but for aesthetic reasons and for clarity of how this function works, I prevent this from occurring.
        # Perform the chunking
        ongoingCount = 0    # This will keep track of what sequence number we are on
        records = SeqIO.parse(open(inputFasta, 'r'), 'fasta')
        fileNames = []
        for i in range(threads):
                # Flow control
                if ongoingCount == numSeqs: # This lets us stop making new files if we have more threads than we do sequences
                        break
                # Generate the file name
                chunkName = os.path.join(outputDir, os.path.basename(inputFasta).rsplit('.', maxsplit=1)[0] + fastaSuffix + str(i+1) + '.fasta')
                fileNames.append(chunkName)
                # Write sequences to chunk file
                with open(chunkName, 'w') as fileOut:
                        for record in records:      # We'll run out of records if we get to a point where ongoingCount == numSeqs
                                fileOut.write('>' + record.description + '\n' + str(record.seq) + '\n')
                                ongoingCount += 1
                                if ongoingCount in chunkPoints:
                                        break
        return fileNames

def coord_cutter(fastaFile, coordDict, outputFileName):
        # Set up
        from Bio import SeqIO
        # Remove domain regions from fasta
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        with open(outputFileName, 'w') as outFile:
                for record in records:
                        # Processing steps
                        seqName = record.description
                        if seqName not in coordDict:
                                outFile.write('>' + seqName + '\n' + str(record.seq) + '\n')
                        else:
                                # Mask the domain regions
                                currSeq = str(record.seq)
                                coords = coordDict[seqName]
                                if coords != []:                        # Empty lists are fine, they won't trigger our coord loop so the sequence will be unaltered
                                        if type(coords[0]) == int:      # If this is true, our structure looks something like [1,2]; this function expects a list of lists e.g., [[1,2]]
                                                coords = [coords]       # By doing this, we can make this function compatible with other coord dict structures which only ever have 1 entry per sequence
                                for coord in coords:
                                        currSeq = currSeq[0:coord[0]-1] + ('x' * (coord[1] + 1 - coord[0])) + currSeq[coord[1]:]        # -1 to coord[0] to make it 0-based; +1 to coord[1] since a domain range of 1-1 still has a length of 1;
                                outFile.write('>' + seqName + '\n' + currSeq + '\n')

def fasta_domain_extract(coordDict, fastaFile, outputFileName, minSize):
        # Setup
        from Bio import SeqIO
        # Load input file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Produce output
        with open(outputFileName, 'w') as fileOut:
                for record in records:
                        if record.description in coordDict:
                                seqid = record.description
                        elif record.id in coordDict:    # POSSIBLE PROBLEM: May need to be more strict with ID parsing consistency especially since MMseqs2 does change sequence IDs...
                                seqid = record.id
                        else:
                                continue
                        seq = str(record.seq)
                        ranges = coordDict[seqid]
                        ongoingCount = 1
                        for i in range(len(ranges)):
                                # Minimum size cut-off enforcement
                                if ranges[i][1] - ranges[i][0] + 1 < minSize:
                                        continue
                                # Extract sequence
                                tmpDomain = seq[ranges[i][0]-1:ranges[i][1]]
                                fileOut.write('>' + record.description + '_Domain_' + str(ongoingCount) + '_' + str(ranges[i][0]) + '-' + str(ranges[i][1]) + '\n' + tmpDomain + '\n')
                                ongoingCount += 1

### CD-HIT
def run_cdhit(cdhitDir, outputDir, inputFasta, outputFasta, params):
        import os, subprocess
        # Format command
        cmd = os.path.join(cdhitDir, 'cd-hit') + ' -i ' + inputFasta + ' -o ' + os.path.join(outputDir, outputFasta) + ' -c {} -n {} -G {} -aS {} -aL {} -M {} -T {}'.format(*params)
        # Run CD-HIT
        run_cdhit = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        cdout, cderr = run_cdhit.communicate()
        if cderr.decode("utf-8") != '':
                raise Exception('CD-HIT Error text below' + str(cderr.decode("utf-8")))

### HMMER3
def run_hmmer3(hmmer3dir, hmmDB, outputDir, threads, evalue, inputFasta, outputFileName):
        import os, subprocess
        # Check if we need to run hmmpress
        if not os.path.isfile(hmmDB + '.h3f') and not os.path.isfile(hmmDB + '.h3i') and not os.path.isfile(hmmDB + '.h3m') and not os.path.isfile(hmmDB + '.h3p'):
                cmd = os.path.join(hmmer3dir, 'hmmpress') + ' -f "' + hmmDB + '"'
                run_hmmpress = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                hmmout, hmmerr = run_hmmpress.communicate()
                if hmmerr.decode("utf-8") != '':
                        raise Exception('hmmpress error text below' + str(hmmerr.decode("utf-8")))
        # Run HMMER3
        cmd = os.path.join(hmmer3dir, 'hmmsearch') + ' --cpu ' + str(threads) + ' -E ' + str(evalue) + ' --domtblout ' + outputFileName + ' "' + hmmDB + '" "' + inputFasta + '"'
        run_hmmer3 = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        hmmout, hmmerr = run_hmmer3.communicate()
        if hmmerr.decode("utf-8") != '':
                raise Exception('hmmsearch error text below' + str(hmmerr.decode("utf-8")))

def hmmer_parse_domfind(domtbloutFile, evalueCutoff, skip):
        # Set up
        import os
        domDict = {}                            # We need to use a dictionary for later sorting since hmmsearch does not produce output that is ordered in the way we want to work with. hmmscan does, but it is SIGNIFICANTLY slower.
        # Main function
        with open(domtbloutFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip unnecessary lines
                        if line.startswith('#') or line == '' or line == ' ' or line == '\n' or line == '\r\n':
                                continue
                        # Parse line and skip if evalue is not significant
                        sl = line.rstrip('\r\n').split()
                        evalue = float(sl[12])
                        if evalue > float(evalueCutoff):
                                continue
                        # Get relevant details
                        pid = sl[0]
                        dstart = int(sl[17])
                        dend = int(sl[18])
                        # Special handling for domain IDs
                        if sl[3].startswith('cath'):
                                did = sl[3]
                        else:
                                did = os.path.basename(sl[3])                   # Some of the databases will have the full path to the domain ID, so we want to get rid of those. We handle cath especially since cath's domain IDs don't exhibit this behaviour but do incorporate a '/' in their ID
                        # Optional skipping of CATH/SUPERFAMILY databases
                        if skip != 'noskip':
                                if skip == 'cath' and sl[3].startswith('cath'):
                                        continue
                                elif skip == 'superfamily':
                                        try:
                                                int(did)                # SUPERFAMILY is the only database that has purely integer domain IDs
                                                continue
                                        except ValueError:
                                                None
                                elif skip == 'both':
                                        if sl[3].startswith('cath'):
                                                continue
                                        try:
                                                int(did)
                                                continue
                                        except ValueError:
                                                None
                        # Add into domain dictionary
                        if pid not in domDict:
                                domDict[pid] = [[did, dstart, dend, evalue]]
                        else:
                                domDict[pid].append([did, dstart, dend, evalue])
        return domDict

### SIGNALP (Loads in CD-HIT results for running program, then modifies HMMER3 cut results)
def run_signalp(signalpdir, cygwindir, outputDir, outputFileName, organism, fileNames):
        # Set up
        import threading, os
        # Define functions integral to this one
        def signalp_thread(signalpdir, cygwindir, outputDir, organism, fastaFile, resultNames):
                import os, subprocess, platform
                # Get the full fasta file location
                fastaFile = os.path.abspath(fastaFile)
                # Format signalP script text
                sigpResultFile = os.path.join(outputDir, thread_file_name_gen('tmp_sigpResults_' + os.path.basename(fastaFile), ''))
                scriptText = '"' + os.path.join(signalpdir, 'signalp') + '" -t ' + organism + ' -f short -n "' + sigpResultFile + '" "' + fastaFile + '"'
                # Generate a script for use with cygwin (if on Windows)
                if platform.system() == 'Windows':
                        sigpScriptFile = os.path.join(outputDir, thread_file_name_gen('tmp_sigpScript_' + os.path.basename(fastaFile), '.sh'))
                        with open(sigpScriptFile, 'w') as fileOut:
                                fileOut.write(scriptText.replace('\\', '/'))
                # Run signalP depending on operating system
                if platform.system() == 'Windows':
                        cmd = os.path.join(cygwindir, 'bash') + ' -l -c "' + sigpScriptFile.replace('\\', '/') + '"'
                        runsigP = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                        sigpout, sigperr = runsigP.communicate()
                        os.remove(sigpScriptFile)       # Clean up temporary file
                else:
                        runsigP = subprocess.Popen(scriptText, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                        sigpout, sigperr = runsigP.communicate()
                # Process output
                okayLines = ['is an unknown amino amino acid', 'perl: warning:', 'LC_ALL =', 'LANG =', 'are supported and installed on your system']
                for line in sigperr.decode("utf-8").split('\n'):
                        # If sigperr indicates null result, create an output file we can skip later
                        if line.rstrip('\n') == '# No sequences predicted with a signal peptide':
                                with open(sigpResultFile, 'w') as fileOut:
                                        fileOut.write(line)
                                break
                        # Check if this line has something present within okayLines
                        okay = 'n'
                        for entry in okayLines:
                                if entry in line or line == '':
                                        okay = 'y'
                                        break
                        if okay == 'y':
                                continue
                        # If nothing matches the okayLines list, we have a potentially true error
                        else:
                                raise Exception('SignalP error occurred when processing file name ' + fastaFile + '. Error text below\n' + sigperr.decode("utf-8"))
                # Store the result file name in a mutable object so we can retrieve it after joining
                resultNames.append(sigpResultFile)
        # Main function
        # Run signalP on each of the input files
        resultNames = []
        processing_threads = []
        for name in fileNames:
                build = threading.Thread(target=signalp_thread, args=(signalpdir, cygwindir, outputDir, organism, name, resultNames))
                processing_threads.append(build)
                build.start()
        # Wait for all threads to end
        for process_thread in processing_threads:
                process_thread.join()
        # Join and parse signalP results files
        combinedFile = ''
        for name in resultNames:
                with open(name, 'r') as fileIn:
                        for line in fileIn:
                                if not line.endswith('\n'):
                                        line += '\n'
                                combinedFile += line
        # Clean up temporary files
        for name in resultNames:
                os.remove(name)
        # Write main output file
        with open(outputFileName, 'w') as fileOut:
                fileOut.write(combinedFile)

def parse_sigp_results(sigpFile):
        sigPredictions = {}
        with open(sigpFile, 'r') as fileIn:
                for line in fileIn:
                        if line.startswith('#'):
                                continue
                        sl = line.split('\t')
                        sigPredictions[sl[0]] = [int(sl[3]), int(sl[4])]
        # Return signalP prediction dictionary
        return sigPredictions

### SEG AND COILS (Loads in signalP masked file)
def run_seg(segdir, outputDir, fileNames, outputFileName):
        import threading, os
        # Define functions integral to this one
        def seg_thread(segdir, fastaFile, resultNames):
                # Set up
                import subprocess
                # Get the full fasta file location & derive our output file name
                fastaFile = os.path.abspath(fastaFile)
                segResultFile = os.path.join(outputDir, thread_file_name_gen('tmp_segResults_' + os.path.basename(fastaFile), ''))
                # Format seg command and run
                cmd = os.path.join(segdir, 'seg') + ' "' + fastaFile + '" -x > ' + '"' + segResultFile + '"'
                #print(cmd)
                runseg = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                segout, segerr = runseg.communicate()
                # Process output
                if segerr.decode("utf-8") != '':
                        raise Exception('SEG error text below\n' + segerr.decode("utf-8"))
                # Store the result file name in a mutable object so we can retrieve it after joining
                resultNames.append(segResultFile)
        # Main function
        # Run seg on each of the input files
        processing_threads = []
        resultNames = []        # Use a mutable list here so we can retrieve the file names in the absence of being able to return these through the threaded function
        for name in fileNames:
                build = threading.Thread(target=seg_thread, args=(segdir, name, resultNames))
                processing_threads.append(build)
                build.start()
        # Wait for all threads to end
        for process_thread in processing_threads:
                process_thread.join()
        # Join seg results files
        combinedFile = ''
        for name in resultNames:
                with open(name, 'r') as fileIn:
                        for line in fileIn:
                                if line == '\r\n' or line == '\n':
                                        continue
                                if not line.endswith('\n'):
                                        line += '\n'
                                combinedFile += line
        # Clean up temporary files
        for name in resultNames:
                os.remove(name)
        # Write main output file
        with open(outputFileName, 'w') as fileOut:
                fileOut.write(combinedFile)

def parse_seg_results(segFile):
        # Set up
        from Bio import SeqIO
        # Parse seg results files
        segPredictions = {}
        segRecords = SeqIO.parse(open(segFile, 'r'), 'fasta')
        for record in segRecords:
                seqid = record.description
                seq = str(record.seq)
                xCoords = consecutive_character_coords(seq, 'x', 1, 'pairs')
                segPredictions[seqid] = xCoords
        # Return seg prediction dictionary
        return segPredictions

def run_coils(coilsdir, py2dir, fileNames, outputFileName):
        # Set up
        import threading
        # Define functions integral to this one
        def coils_thread(coilsdir, py2dir, fastaFile, coilsResults):
                import os, subprocess
                # Get the full fasta file location & derive our output file name
                fastaFile = os.path.abspath(fastaFile)
                # Format coils command & run
                cmd = '"' + os.path.join(py2dir, 'python') + '" "' + os.path.join(coilsdir, 'psCoils.py') + '" -f "' + fastaFile + '"'
                runcoils = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
                coilsout, coilserr = runcoils.communicate()
                # Process output
                if coilserr.decode("utf-8") != '':
                        raise Exception('Coils error text below\n' + coilserr.decode("utf-8"))
                # Store the result file name in a mutable object so we can retrieve it after joining
                coilsResults.append(coilsout.decode("utf-8"))
        # Main function
        # Run coils on each of the input files
        processing_threads = []
        coilsResults = []        # Use a mutable list here so we can retrieve the file names in the absence of being able to return these through the threaded function
        for name in fileNames:
                build = threading.Thread(target=coils_thread, args=(coilsdir, py2dir, name, coilsResults))
                processing_threads.append(build)
                build.start()
        # Wait for all threads to end
        for process_thread in processing_threads:
                process_thread.join()
        # Write results to file
        with open(outputFileName, 'w') as fileOut:
                for result in coilsResults:
                        result = result.replace('\r', '')
                        for line in result.split('\n'):
                                if line == '':
                                        continue
                                fileOut.write(line + '\n')

def parse_coils_results(coilsFile, fastaFiles):
        # Set up
        from Bio import SeqIO
        coilsPredictions = {}
        # Main function
        coilsResults = open(coilsFile, 'r').read()
        # Split result by headers - each header corresponds to a sequence's result
        result = coilsResults.split(' Pos A Hep Score   Prob    Gcc     Gg    Pred (Loop=L Coiledcoil=C)')
        while '' in result:     # There should only be one entry corresponding to this at the very start of the result list
                del result[result.index('')]
        for i in range(len(result)):
                # Build a sequence consisting of L's (loops) and C's (coils) in addition to the original sequence
                coilSeq = ''
                protSeq = ''
                for row in result[i].split('\n'):
                        if row == '':
                                continue
                        sr = row.split()
                        coilSeq += sr[7]
                        protSeq += sr[1]
                # Extract coil coordinates
                cCoords = consecutive_character_coords(coilSeq, 'C', 1, 'pairs')
                # Add to our coilsPredictions dictionary
                coilsPredictions[protSeq] = cCoords     # psCoils doesn't provide ordered results, so we need to match protein sequences to their coil results
        # Associate coils results to sequence IDs
        coilsIDs = {}
        for name in fastaFiles:
                records = SeqIO.parse(open(name, 'r'), 'fasta')
                for record in records:
                        coilsIDs[record.description] = coilsPredictions[str(record.seq)]
        # Return coils prediction dictionary
        return coilsIDs

## FUNCTIONS RELATING TO MMseqs2
def makemms2db(mmseqs2dir, query, target, which):
        import os, subprocess
        # Format command
        dbname1 = query + '_queryDB'
        cmd1 = os.path.join(mmseqs2dir, 'mmseqs') + ' createdb "' + query + '" "' + dbname1 + '"'
        if target != None:
                dbname2 = target + '_targetDB'
                cmd2 = os.path.join(mmseqs2dir, 'mmseqs') + ' createdb "' + target + '" "' + dbname2 + '"'
        # Query DB generation
        if which == 'query' or which == 'both':
                run_makedb = subprocess.Popen(cmd1, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                makedbout, makedberr = run_makedb.communicate()
                if makedberr.decode("utf-8") != '':
                        raise Exception('Make MMseqs2 query db error text below\n' + makedberr.decode("utf-8"))
        # Run target DB generation if target != query
        if query != target and target != None:  # This lets us use this function when we know we don't have a target; we can just specify None as a more intuitive way of producing the expected behaviour of only working with a single query
                if which == 'target' or which == 'both':
                        run_makedb = subprocess.Popen(cmd2, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                        makedbout, makedberr = run_makedb.communicate()
                        if makedberr.decode("utf-8") != '':
                                raise Exception('Make MMseqs2 target db error text below\n' + makedberr.decode("utf-8"))

def indexmms2(mmseqs2dir, query, target, tmpdir, threads, which):
        import os, subprocess
        # Format command
        dbname1 = query + '_queryDB'
        cmd1 = os.path.join(mmseqs2dir, 'mmseqs') + ' createindex "' + dbname1 + '" "' + tmpdir + '" --threads ' + str(threads)
        if target != None:
                dbname2 = target + '_targetDB'
                cmd2 = os.path.join(mmseqs2dir, 'mmseqs') + ' createindex "' + dbname2 + '" "' + tmpdir + '" --threads ' + str(threads)
        # Run query index
        if which == 'query' or which == 'both':
                run_index = subprocess.Popen(cmd1, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                indexout, indexerr = run_index.communicate()
                if indexerr.decode("utf-8") != '':
                        raise Exception('Indexing MMseqs2 query db error text below\n' + indexerr.decode("utf-8"))
        # Run target DB indexing if target != query
        if query != target and target != None:
                if which == 'target' or which == 'both':
                        run_index = subprocess.Popen(cmd2, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                        indexout, indexerr = run_index.communicate()
                        if indexerr.decode("utf-8") != '':
                                raise Exception('Indexing MMseqs2 target db error text below\n' + indexerr.decode("utf-8"))

def runmms2(mmseqs2dir, query, target, tmpdir, searchName, params):
        import os, subprocess
        # Format command
        dbname1 = query + '_queryDB'
        if query != target and target != None:
                dbname2 = target + '_targetDB'
        else:
                dbname2 = query + '_queryDB'
        cmd = os.path.join(mmseqs2dir, 'mmseqs') + ' search "' + dbname1 + '" "' + dbname2 + '" "' + searchName + '" "' + tmpdir + '" -e {} --threads {} --num-iterations {} -s {} --alt-ali {}'.format(*params)
        print('#' + cmd)
        # Run query
        run_mms2 = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        mms2out, mms2err = run_mms2.communicate()
        if mms2err.decode("utf-8") != '':
                raise Exception('MMseqs2 search error text below\n' + mms2err.decode("utf-8"))

def mms2tab(mmseqs2dir, query, target, searchName, threads):
        import os, subprocess
        # Get file details
        dbname1 = query + '_queryDB'
        if query != target and target != None:
                dbname2 = target + '_targetDB'
        else:
                dbname2 = query + '_queryDB'
        # Create tab-delim BLAST-like output
        cmd = os.path.join(mmseqs2dir, 'mmseqs') + ' convertalis "' + dbname1 + '" "' + dbname2 + '" "' + searchName + '" "' + searchName + '.m8" --threads ' + str(threads)
        run_mms2 = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        mms2out, mms2err = run_mms2.communicate()
        if mms2err.decode("utf-8") != '':
                raise Exception('MMseqs2 tabular output generation error text below\n' + mms2err.decode("utf-8"))

def mms2sort_all(searchName, outName):
        from itertools import groupby
        groups = []
        # Parse file
        grouper = lambda x: x.split('\t')[0]
        with open(searchName, 'r') as fileIn:
                for key, group in groupby(fileIn, grouper):
                        group = list(group)
                        # Sort group [If len(group) == 1 it doesn't matter]
                        for i in range(len(group)):
                                group[i] = group[i].rstrip('\n').split('\t')
                        group.sort(key = lambda x: (float(x[10]),-float(x[11])))
                        # Get a representative E-value for the group
                        groupEvalue = 100000 # Handles groups with only self hit
                        for i in range(len(group)):
                                if group[i][0] != group[i][1] and groupEvalue == 100000:
                                        groupEvalue = float(group[i][10])
                                group[i] = '\t'.join(group[i])
                        # Store group for overall sorting
                        groups.append([groupEvalue, group])
        # Sort overall groups
        groups.sort(key = lambda x: x[0])
        # Put in output
        with open(outName, 'w') as fileOut:
                for entry in groups:
                        fileOut.write("\n".join(entry[1]) + '\n')

def parsemms2tab_to_array(mms2Table, fastaFile, lowLenCutoff):
        # Set up
        from Bio import SeqIO
        from itertools import groupby
        seqArrays = {}
        grouper = lambda x: x.split('\t')[0]
        # Get the sequence lengths for each query
        records = SeqIO.parse(open(fastaFile, 'rU'), 'fasta')
        for record in records:
                seqId = record.id
                seqLen = len(str(record.seq))
                seqArrays[seqId] = [0]*seqLen         # This can be converted into a numpy array later
        # Load in file as a groupby iterator
        with open(mms2Table, 'r') as fileIn:
                # Parse pblast file to pull out coverage per position arrays
                for key, group in groupby(fileIn, grouper):
                        for entry in group:
                                if entry == '\n':
                                        continue
                                line = entry.split('\t')
                                # Extract details
                                qCoord = (int(line[6]), int(line[7]))
                                tCoord = (int(line[8]), int(line[9]))
                                # Allow discovery of internal repeats by looking to see if there is overlap between self hits
                                if line[0] == line[1]:
                                        ovl = min(qCoord[1], tCoord[1]) - max(qCoord[0], tCoord[0]) + 1         # +1 since we're working 1-based here; 300-250==50, but this range is inclusive of 250 and 300, so it should be 51
                                        if ovl / (qCoord[1] - qCoord[0] + 1) > 0.05 and ovl > 5:                # If there is more than 5% or 5 AA overlap (to account for shorter domains), we reject it as a hit
                                                continue
                                # Handle any hits which reach here by checking overlap length
                                if (qCoord[1] - qCoord[0] + 1) > lowLenCutoff:                                  # This will mean we only find putative globular domains
                                        for i in range(qCoord[0]-1, qCoord[1]):
                                                seqArrays[line[0]][i] += 1
                                        for i in range(tCoord[0]-1, tCoord[1]):
                                                seqArrays[line[1]][i] += 1
        return seqArrays

def parse_array_peaks(seqArrays, minPlateauSize):
        # Setup
        import math
        import numpy as np
        domDict = {}
        covRatio = 0.33         # This is an arbitrary value which controls the coverage cutoff for extending and chaining plateaus together
        # Define functions integral to this one
        def plateau_extens(array, plateaus, coverages, covRatio):
                import math
                ## We can handle plateau extension without causing incorrect overlap by checking for the first of two occurrences. 1: we find a point where the length cutoff becomes enforced, or 2: the coverage starts to increase again, which means we're heading towards another peak (which wasn't collapsed into this plateau).
                for i in range(len(plateaus)):
                        covCutoff = math.ceil(coverages[i] * covRatio)
                        # Look back
                        prevCov = array[plateaus[i][0]]
                        newStart = plateaus[i][0]
                        for x in range(plateaus[i][0]-1, -1, -1):
                                indexCov = array[x]
                                # Increasing check
                                if indexCov > prevCov:  # This means we're leading up to another peak, and should stop extending this plateau
                                        break
                                # Decreasing cut-off
                                if indexCov >= covCutoff:
                                        newStart = x
                                        prevCov = indexCov
                                        continue
                                else:
                                        break
                        plateaus[i][0] = newStart
                        # Look forward
                        prevCov = array[plateaus[i][1]]
                        newEnd = plateaus[i][1]
                        for x in range(plateaus[i][1]+1, len(array)):
                                indexCov = array[x]
                                # Increasing check
                                if indexCov > prevCov:  # This means we're leading up to another peak, and should stop extending this plateau
                                        break
                                # Decreasing cut-off
                                if indexCov >= covCutoff:
                                        newEnd = x
                                        prevCov = indexCov
                                        continue
                                else:
                                        break
                        plateaus[i][1] = newEnd
                return plateaus
        
        def plateau_chain_seeds(array, plateaus, coverages, covRatio):
                # Sort plateaus based on coverages
                plateaus = [x for _,x in sorted(zip(coverages,plateaus), key=lambda pair: -pair[0])]            # https://stackoverflow.com/questions/6618515/sorting-list-based-on-values-from-another-list
                coverages.sort(reverse = True)                                                                  # We don't need to do any fancy sorting for coverages, any duplicates will sit where they should
                # Perform pairwise comparison loop to merge plateaus on the basis of coverage
                for y in range(len(plateaus)-1):
                        z = y + 1
                        while True:
                                # Exit condition
                                if z >= len(plateaus):                                  # len(plateaus)-1 would correspond to the final entry, this means we've gone at least one step further beyond
                                        break
                                # Extract depression details;                           # S=start,E=end,A=array... I use 'depression' to refer to dips in the coverage array (i.e., we're going away from a peak so it will be decreasing)
                                depressS = min(plateaus[y][1], plateaus[z][1]) + 1      # Since we know these values should never overlap, this provides a way to determine the 'gap'/'depression' region between the peaks; +1 to look at the first position AFTER the peak i.e., the first base of the depression
                                depressE = max(plateaus[y][0], plateaus[z][0]) - 1      # -1 to look at the first position BEFORE the next peak i.e., the last base of the depression
                                depressA = array[depressS:depressE + 1]
                                # Use covRatio to determine if chaining is correct
                                maxVal = max([coverages[y], coverages[z]])              # We determine chaining based on the highest coverage plateau for a few reasons. Most importantly, the higher coverage plateau is most likely to be a real domain, so we should prioritise it and make sure it doesn't join to lower coverage regions unless they meet our arbitrary chaining limits
                                cutoff = math.ceil(maxVal * covRatio)                   # We round up to handle low numbers properly. For example, 2*0.5=1.0. We don't need to round this, and 1.0 is good here if depressL<= cleanAA. However, 2*0.75=1.5. Rounding down would be 1, but this would mean we chain probably separate domains together incorrectly. Thus, by rounding up, we are more strict.
                                belowCutoff = np.where(depressA < cutoff)
                                aboveMax = np.where(depressA > maxVal)
                                if len(belowCutoff[0]) > 0 or len(aboveMax[0]) > 0:     # We don't chain these two plateaus together if we have any positions with less coverage than our cutoff/have regions with no coverage OR if we have a peak inbetween the two plateaus we are currently comparing
                                        z += 1
                                else:
                                        # Merge this pair's plateaus/coverages values
                                        plateaus[y] = [min(plateaus[y][0], plateaus[z][0]), max(plateaus[y][1], plateaus[z][1])]
                                        coverages[y] = maxVal                                   # It's important that we make coverages == maxVal since this will prevent us from chaining plateaus together in steps (e.g., 5 coverage plateau merges with a 4 cov, then that 4 cov merges with a 3, etc.)
                                        del plateaus[z]
                                        del coverages[z]
                # Unsort coverages based on plateaus & return
                coverages = [x for _,x in sorted(zip(plateaus,coverages), key=lambda pair: pair[0])]
                plateaus.sort() 
                return plateaus, coverages
        
        def one_based_index_fix(coords):
                for i in range(len(coords)):
                        coords[i][0] += 1
                        coords[i][1] += 1
                return coords
        
        # Main function
        for key, vlist in seqArrays.items():
                # Check if we got any hits
                if sum(vlist) == 0:            # This means we had no accepted hits (i.e., those which pass E-value and self-hit checks)
                        continue
                elif len(set(vlist)) == 1:     # This means we only had completely overlapping hits against this sequence  ## CHECK THE OUTCOME OF THIS ##
                        continue
                # Convert to an array
                array = np.array(vlist + [0])                                   # We need to add a 0 to the array at the end since the peakdetect algorithm doesn't work if a plateau runs to the end of the array... I didn't make the function and don't understand it, so I'll just work around it
                # Find peak indices                                             # Return format is as [[peakLocation, peakValue]]
                maxindices, minindices = peakdet(array,0.5)                     # 0.5 means that we'll detect peaks that are at least > 0.5 than the preceding number; in effect, this means we'll detect "peaks" which are 1 number higher than their neighbour which will be all "peaks" in our array since we solely use integers
                # Get our accurate array now
                array = np.array(vlist)
                # Get plateau regions
                plateaus = []           # Plateaus lists the coordinate ranges of our identified peaks
                coverages = []          # Coverages acts as a paired list to plateaus, listing the coverage number of the plateau
                for maximum in maxindices:
                        index = maximum[0]
                        coverage = maximum[1]
                        # Look forward                          # The peakdet values are always at the start of the plateau, so we don't need to look back, we just need to look forward to find where the plateau ends]
                        plat = None
                        for i in range(index, len(array)):
                                if array[i] == coverage:
                                        continue
                                else:
                                        plat = [index,i-1]      # This is 0-indexed, and we -1 since we want the previous i value
                                        break
                        if plat == None:                        # This acts as a check for plateaus that run to the end of the sequence
                                plat = plat = [index,i]         # We don't -1 here since i will be equal to the last position of the sequence (in 0-based notation)
                        plateaus.append(plat)
                        coverages.append(coverage)
                # Chain and clean up multiple plateaus
                if len(plateaus) > 1:
                        plateaus, coverages = plateau_chain_seeds(array, plateaus, coverages, covRatio)
                # Extend plateaus following slightly modified chaining rules and add to our domDict for later output generation
                plateaus = plateau_extens(array, plateaus, coverages, covRatio)
                plateaus = one_based_index_fix(plateaus)        # This is necessary since we're going back from working with 0-based numpy arrays to a 1-based BLAST-like format
                # Ensure the algorithm is working properly [TESTING: Can be stripped out once I'm sure the algorithm is solid, it'll just slow things down otherwise]
                for i in range(len(plateaus)-1):
                        if plateaus[i][1] > plateaus[i+1][0] and plateaus[i+1][1] > plateaus[i][0]:
                                print('parse_array_peaks: Algorithm problem, this should never happen')
                                print(plateaus)
                                print(coverages)
                                print(vlist)
                # Remove plateaus that do not meet minimum length cut-off
                for i in range(len(plateaus)-1, -1, -1):
                        platLen = plateaus[i][1] - plateaus[i][0] + 1           # +1 since [100,100] is a plateau with length of 1
                        if platLen < minPlateauSize:
                                del plateaus[i]
                # Add to main dictionary
                domDict[key] = plateaus
        return domDict

def peakdet(v, delta, x = None):
        import sys
        from numpy import NaN, Inf, arange, isscalar, asarray, array
        """
        Converted from MATLAB script at http://billauer.co.il/peakdet.html
        
        Returns two arrays
        
        function [maxtab, mintab]=peakdet(v, delta, x)
        %PEAKDET Detect peaks in a vector
        %                [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
        %                maxima and minima ("peaks") in the vector V.
        %                MAXTAB and MINTAB consists of two columns. Column 1
        %                contains indices in V, and column 2 the found values.
        %          
        %                With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
        %                in MAXTAB and MINTAB are replaced with the corresponding
        %                X-values.
        %
        %                A point is considered a maximum peak if it has the maximal
        %                value, and was preceded (to the left) by a value lower by
        %                DELTA.
        
        % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
        % This function is released to the public domain; Any use is allowed.
        
        """
        maxtab = []
        mintab = []
           
        if x is None:
                x = arange(len(v))
        
        v = asarray(v)
        
        if len(v) != len(x):
                sys.exit('Input vectors v and x must have same length')
        
        if not isscalar(delta):
                sys.exit('Input argument delta must be a scalar')
        
        if delta <= 0:
                sys.exit('Input argument delta must be positive')
        
        mn, mx = Inf, -Inf
        mnpos, mxpos = NaN, NaN
        
        lookformax = True
        
        for i in arange(len(v)):
                this = v[i]
                if this > mx:
                        mx = this
                        mxpos = x[i]
                if this < mn:
                        mn = this
                        mnpos = x[i]
                
                if lookformax:
                        if this < mx-delta:
                                maxtab.append((mxpos, mx))
                                mn = this
                                mnpos = x[i]
                                lookformax = False
                else:
                        if this > mn+delta:
                                mintab.append((mnpos, mn))
                                mx = this
                                mxpos = x[i]
                                lookformax = True
        
        return array(maxtab), array(mintab)
