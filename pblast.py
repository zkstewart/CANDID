# Class assistants
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

class domfind:
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
        
        ### HMMER3
        def run_hmmer3(hmmer3dir, hmmDB, outputDir, threads, evalue, inputFasta):
                import os, subprocess
                # Check if we need to run hmmpress
                if not os.path.isfile(hmmDB + '.h3f') and not os.path.isfile(hmmDB + '.h3i') and not os.path.isfile(hmmDB + '.h3m') and not os.path.isfile(hmmDB + '.h3p'):
                        print('Running hmmpress to format your database...')
                        cmd = os.path.join(hmmer3dir, 'hmmpress') + ' -f "' + hmmDB + '"'
                        run_hmmpress = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                        hmmout, hmmerr = run_hmmpress.communicate()
                        if hmmerr.decode("utf-8") != '':
                                raise Exception('hmmpress error text below' + str(hmmerr.decode("utf-8")))
                # Run HMMER3
                fastaBase = os.path.basename(inputFasta).rsplit('.', maxsplit=1)[0]
                cmd = os.path.join(hmmer3dir, 'hmmsearch') + ' --cpu ' + str(threads) + ' -E ' + str(evalue) + ' --domtblout ' + os.path.join(outputDir, fastaBase + '_hmmer.results') + ' "' + hmmDB + '" "' + inputFasta + '"'
                print('Running hmmsearch to detect known domains...')
                run_hmmer3 = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                hmmout, hmmerr = run_hmmer3.communicate()
                if hmmerr.decode("utf-8") != '':
                        raise Exception('hmmsearch error text below' + str(hmmerr.decode("utf-8")))

        def hmmer_parse(domtbloutFile, evalueCutoff, skip):
                # Set up
                import os
                domDict = {}                            # We need to use a dictionary for later sorting since hmmsearch does not produce output that is ordered in the way we want to work with. hmmscan does, but it is SIGNIFICANTLY slower.
                print('Parsing hmmsearch output...')
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
        
        def hmmer_dict_to_file(domDict, outputFileName):
                with open(outputFileName, 'w') as fileOut:
                        for key, value in domDict.items():
                                for domregion in value:
                                        fileOut.write(key + '\t' + str(domregion[1]) + '\t' + str(domregion[2]) + '\t' + domregion[0] + '\n')

        def hmmer_coord_parse(parsedHmmerFile):
                # Set up
                from itertools import groupby
                # Declare function integral to this function
                def coord_merge(coordList, coord):
                        # Merge the new coord into the current coordList
                        merged = 'n'
                        if coord != None:
                                for i in range(len(coordList)):
                                        pair1 = coordList[i]
                                        pair2 = coord
                                        # Detect overlap
                                        if pair1[1] >= pair2[0] and pair2[1] >= pair1[0]:
                                                # Merge coords
                                                start = min([pair1[0], pair2[0]])
                                                end = max([pair1[1], pair2[1]])
                                                coordList[i] = [start, end]
                                                merged = 'y'
                                                break
                        # If we didn't merge this coord into an existing one, add it into the list
                        if merged == 'n' and coord != None:
                                coordList.append(coord)
                        # If we did merge it, re-process the coordList to merge anything else that needs it
                        else:
                                for x in range(len(coordList)-1,-1,-1):
                                        if x != 0:
                                                pair1 = coordList[x]
                                                pair2 = coordList[x-1]
                                                # Detect overlap
                                                if pair1[1] >= pair2[0] and pair2[1] >= pair1[0]:
                                                        # Merge coords
                                                        start = min([pair1[0], pair2[0]])
                                                        end = max([pair1[1], pair2[1]])
                                                        coordList[x-1] = [start, end]
                                                        # Cull entry
                                                        del coordList[x]
                        # Sort coordList and return
                        coordList.sort()
                        return coordList
                # Main function
                hmmerDict = {}
                grouper = lambda x: x.split('\t')[0]
                with open(parsedHmmerFile, 'r') as fileIn:
                        for key, group in groupby(fileIn, grouper):
                                coords = []
                                for line in group:
                                        sl = line.rstrip('\r\n').split('\t')
                                        coord_merge(coords, [int(sl[1]), int(sl[2])])
                                # Retain results for subsequent cutting in a dictionary
                                hmmerDict[key] = coords
                return hmmerDict
        
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

        ######## FUNCTIONS RELATING TO MMseqs2
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
                print(cmd)
                # Run query
                run_mms2 = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                mms2out, mms2err = run_mms2.communicate()
                if mms2err.decode("utf-8") != '':
                        raise Exception('MMseqs2 search error text below\n' + mms2err.decode("utf-8"))
        
        def mms2tab(mmseqs2dir, query, target, tmpdir, searchName, threads):
                import os, subprocess
                # Get file details
                dbname1 = query + '_queryDB'
                if query != target and target != None:
                        dbname2 = target + '_targetDB'
                else:
                        dbname2 = query + '_queryDB'
                # Create tab-delim BLAST-like output
                cmd = os.path.join(mmseqs2dir, 'mmseqs') + ' convertalis "' + dbname1 + '" "' + dbname2 + '" "' + searchName + '" "' + searchName + '.m8" "' + tmpdir + '" --threads ' + str(threads)
                run_mms2 = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                mms2out, mms2err = run_mms2.communicate()
                if mms2err.decode("utf-8") != '':
                        raise Exception('MMseqs2 tabular output generation error text below\n' + mms2err.decode("utf-8"))

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

        def parsemms2_peaks(seqArrays, lowLenCutoff):
                # Setup
                import math
                import numpy as np
                from .peakdetect import peakdetect
                domDict = {}
                covRatio = 0.75         # This is an arbitrary value which controls the coverage cutoff for extending and chaining plateaus together
                # Define functions integral to this one
                def plateau_extens(array, lowLenCutoff, plateaus, coverages, covRatio):
                        import math
                        ## We can handle plateau extension without causing incorrect overlap by checking for the first of two occurrences. 1: we find a point where the length cutoff becomes enforced, or 2: the coverage starts to increase again, which means we're heading towards another peak (which wasn't collapsed into this plateau).
                        for i in range(len(plateaus)):
                                covCutoff = math.ceil(coverages[i] * covRatio)
                                # Look back
                                ongoingCount = 0                        # This measures how long our extension is so we can apply the two arbitrary values when appropriate
                                prevCov = array[plateaus[i][0]]
                                newStart = plateaus[i][0]
                                for x in range(plateaus[i][0]-1, -1, -1):
                                        ongoingCount += 1
                                        indexCov = array[x]
                                        # Increasing check
                                        if indexCov > prevCov:  # This means we're leading up to another peak, and should stop extending this plateau
                                                break
                                        # Decreasing cut-off
                                        if indexCov >= covCutoff:
                                                newStart = x
                                                continue
                                        else:
                                                break
                                plateaus[i][0] = newStart
                                # Look forward
                                ongoingCount = 0
                                prevCov = array[plateaus[i][1]]
                                newEnd = plateaus[i][1]
                                for x in range(plateaus[i][1]+1, len(array)):
                                        ongoingCount += 1
                                        indexCov = array[x]
                                        # Increasing check
                                        if indexCov > prevCov:  # This means we're leading up to another peak, and should stop extending this plateau
                                                break
                                        # Decreasing cut-off
                                        if indexCov >= covCutoff:
                                                newEnd = x
                                                continue
                                        else:
                                                break
                                plateaus[i][1] = newEnd
                        return plateaus

                def plateau_chain(array, plateaus, coverages, covRatio):
                        for i in range(len(plateaus)-1):                        # S=start,E=end,A=array... I use 'depress' here to refer to the 'depression' in the coverage array (i.e., we're going away from a peak so it will be decreasing)
                                depressS = plateaus[i][1] + 1                   # +1 to look at the first position AFTER the peak i.e., the first base of the depression
                                depressE = plateaus[i+1][0] -1                  # -1 to look at the first position BEFORE the next peak i.e., the last base of the depression
                                depressA = array[depressS:depressE + 1]
                                # Use covRatio to determine if chaining is correct
                                maxVal = max([coverages[i], coverages[i+1]])    # We determine chaining based on the highest coverage plateau for a few reasons. Most importantly, the higher coverage plateau is most likely to be a real domain, so we should prioritise it and make sure it doesn't join to lower coverage regions unless they meet our arbitrary chaining limits
                                cutoff = math.ceil(maxVal * covRatio)           # We round up to handle low numbers properly. For example, 2*0.5=1.0. We don't need to round this, and 1.0 is good here if depressL<= cleanAA. However, 2*0.75=1.5. Rounding down would be 1, but this would mean we chain probably separate domains together incorrectly. Thus, by rounding up, we are more strict.
                                belowCutoff = np.where(depressA < cutoff)
                                if len(belowCutoff[0]) > 0:                     # We don't chain these two plateaus together if we have any positions with less coverage than our cutoff/have regions with no coverage
                                        continue                                
                                else:
                                        # Overwrite previous plateaus/coverages values to produce redundancy
                                        plateaus[i] = [plateaus[i][0], plateaus[i+1][1]]        # We can clean up redundancy later rather than dealing with the consequences of deletion during this loop
                                        plateaus[i+1] = [plateaus[i][0], plateaus[i+1][1]]
                                        coverages[i] = maxVal                                   # It's important that we make coverages == maxVal since this will prevent us from chaining plateaus together in steps (e.g., 5 coverage plateau merges with a 4 cov, then that 4 cov merges with a 3, etc.)
                                        coverages[i+1] = maxVal
                        return plateaus, coverages

                def one_based_index_fix(coords):        # This is necessary since we're going back from working with 0-based numpy arrays to a 1-based BLAST-like format
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
                        maxindices, minindices = peakdetect.peakdet(array,0.5)          # 0.5 means that we'll detect peaks that are at least > 0.5 than the preceding number; in effect, this means we'll detect "peaks" which are 1 number higher than their neighbour which will be all "peaks" in our array since we solely use integers
                        # Get our accurate array now
                        array = np.array(vlist)
                        # Get plateau regions
                        plateaus = []           # Plateaus lists the coordinate ranges of our identified peaks
                        coverages = []          # Coverages acts as a paired list to plateaus, listing the coverage number of the plateau
                        for maximum in maxindices:
                                index = maximum[0]
                                coverage = maximum[1]
                                # Look forward                          # The peakdetect.peakdet values are always at the start of the plateau, so we don't need to look back, we just need to look forward to find where the plateau ends]
                                for i in range(index, len(array)):
                                        if array[i] == coverage:
                                                continue
                                        else:
                                                plat = [index,i-1]      # This is 0-indexed, and we -1 since we want the previous i value
                                                break
                                plateaus.append(plat)
                                coverages.append(coverage)
                        # Chain and clean up multiple plateaus
                        if len(plateaus) > 1:
                                plateaus, coverages = plateau_chain(array, plateaus, coverages, covRatio)
                                while True:
                                        overlap = 'n'
                                        for y in range(len(plateaus)-1):        # Plateaus should always maintain its sorting so we can process it in pairwise steps
                                                if plateaus[y+1][0] > plateaus[y][1]:
                                                        continue
                                                elif plateaus[y+1][0] == plateaus[y][1]:
                                                        plateaus[y] = [plateaus[y][0], plateaus[y+1][1]]
                                                        del plateaus[y+1]
                                                        del coverages[y+1]
                                                        overlap = 'y'
                                                        break
                                                else:
                                                        print('This should never happen. Something is wrong with the algorithm.')
                                                        stopplateauclean
                                        if overlap == 'n':          # Exit condition if we make it through the 'for y' loop without encountering any overlaps
                                                break
                        # Extend plateaus following slightly modified chaining rules and add to our domDict for later output generation
                        plateaus = plateau_extens(array, lowLenCutoff, plateaus, coverages, covRatio)
                        plateaus = one_based_index_fix(plateaus)
                        domDict[key] = plateaus
                return domDict

        def parsemms2_nccheck(domDict, lowLenCutoff):
                for key, value in domDict.items():
                        # Remove redundancy
                        value = list(set(value))
                        # Skip processing if only one domain range
                        if len(value) == 1:
                                domDict[key] = value
                                continue
                        value.sort()
                        # Split and collapse domain regions into non-overlapping sequences
                        while True:
                                # Exit condition if multiple domains have been reduced to a single domain
                                if len(value) == 1:
                                        break
                                # Exit condition if no more overlaps are present
                                overlap = 'n'
                                for y in range(len(value)-1):
                                        if value[y+1][0] > value[y][1]:
                                                continue
                                        else:
                                                overlap = 'y'
                                                break
                                if overlap == 'n':
                                        break
                                # Continue processing if overlaps are still present   
                                for i in range(0, len(value)-1):
                                        # Remove redundant domain ranges
                                        if value[i] == value[i+1]:
                                                del value[i+1]
                                                break
                                        # Check for overlaps
                                        if value[i+1][0] > value[i][1]:                                 # If these two domain regions do not overlap...
                                                continue
                                        #### N Check
                                        if value[i][0] != value[i+1][0]:                # If these two domain regions do not start at the same position...
                                                # Get N overhang region
                                                n_value = (value[i][0], value[i+1][0]-1)

                                                # Check if we should split this N overhang off as its own potential domain, or delete it
                                                if n_value[1] - n_value[0] >= lowLenCutoff:
                                                        # Remove N overhang from current value
                                                        value[i] = (n_value[1] + 1, value[i][1])
                                                        # Insert N overhang as new domain
                                                        value.insert(i, n_value)
                                                        # Check if the trimmed current value has been shortened too much
                                                        if value[i+1][1] - value[i+1][0] < lowLenCutoff:
                                                                del value[i+1]
                                                                continue        # Only need to continue here if we insert and delete the old val, since it doesn't disrupt our 'for i' loop
                                                        break                   # We break here if we insert and don't delete the old value, since it will disrupt our 'for i' loop
                                                        
                                                else:
                                                        # Remove N overhang from current value
                                                        value[i] = (value[i+1][0], value[i][1])
                                                        # Check if the trimmed current value has been shortened too much
                                                        if value[i][1] - value[i][0] < lowLenCutoff:
                                                                del value[i]
                                                                break           # We break here if we delete the old value since it disrupts our 'for i' loop. Otherwise, we can continue straight into the C check process
                                        ### C Check
                                        # Contingency to see if the C values are equivalent if we continue here straight from the N check
                                        if value[i][1] == value[i+1][1]:
                                                del value[i]
                                                break
                                        # Find the longer sequence
                                        val1_range = value[i][1] - value[i][0]
                                        val2_range = value[i+1][1] - value[i+1][0]
                                        longer_value = max(val1_range, val2_range)
                                        if longer_value == val1_range:
                                                longer_value = value[i]
                                                shorter_value = value[i+1]
                                                tmp_mark = 0 # Use a temporary marker so we can tell whether the longer_range was == x or == x + 1
                                        else:
                                                longer_value = value[i+1]
                                                shorter_value = value[i]
                                                tmp_mark = 1 # As above
                                        c_value = (shorter_value[1]+1, longer_value[1])
                                        # Check if we should split this C overhang off as its own potential domain, or delete it
                                        if c_value[1] - c_value[0] >= lowLenCutoff:
                                                # Delete the longer sequence
                                                if tmp_mark == 0:
                                                        del value[i]
                                                else:
                                                        del value[i+1]
                                                # Insert C overhang as new domain
                                                value.insert(i+1, c_value)
                                                # Re-sort our domain regions
                                                value.sort()
                                                break
                                        else:
                                                # Delete the longer sequence and restart loop
                                                if tmp_mark == 0:
                                                        del value[i]
                                                else:
                                                        del value[i+1]
                                                break
                        # Update domDict
                        domDict[key] = value
                return domDict

        def fasta_domain_extract(domDict, fastaFile, outputFileName):
                # Setup
                from Bio import SeqIO
                # Load input file
                records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
                # Produce output
                with open(outputFileName, 'w') as fileOut:
                        for record in records:
                                if record.description in domDict:
                                        seqid = record.description
                                elif record.id in domDict:  # POSSIBLE PROBLEM: May need to be more strict with ID parsing consistency especially since MMseqs2 does change sequence IDs...
                                        seqid = record.id
                                else:
                                        continue
                                seq = str(record.seq)
                                ranges = domDict[seqid]
                                for i in range(len(ranges)):
                                        tmpDomain = seq[ranges[i][0]-1:ranges[i][1]]
                                        fileOut.write('>' + record.description + '_Domain_' + str(i+1) + '_' + str(ranges[i][0]) + '-' + str(ranges[i][1]) + '\n' + tmpDomain + '\n')
