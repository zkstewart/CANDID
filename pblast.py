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
                import os
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
                chunkSize = int(numSeqs / threads) + (numSeqs % threads > 0)        # Need to round up
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
                                        if ongoingCount % chunkSize == 0:
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
        
        def hmmer_cutter(fastaFile, hmmerCoordDict, outputFileName):
                # Set up
                from Bio import SeqIO
                # Remove domain regions from fasta
                records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
                with open(outputFileName, 'w') as outFile:
                        for record in records:
                                # Processing steps
                                seqName = record.id
                                if seqName not in hmmerCoordDict:
                                        outFile.write('>' + seqName + '\n' + str(record.seq) + '\n')
                                else:
                                        # Mask the domain regions
                                        currSeq = str(record.seq)
                                        for coord in hmmerCoordDict[seqName]:
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
                
                def thread_file_name_gen(prefix, threadNum):
                        ongoingCount = 0
                        while True:
                                if not os.path.isfile(prefix + threadNum):
                                        return prefix + threadNum
                                elif os.path.isfile(prefix + threadNum + '.' + str(ongoingCount)):
                                        ongoingCount += 1
                                else:
                                        return prefix + threadNum + '.' + str(ongoingCount)
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
        
        def mask_fasta_by_sigp(sigpDict, fastaFile, outputFileName):
                # Set up
                from Bio import SeqIO
                # Mask signal peptides
                records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
                with open(outputFileName, 'w') as fileOut:
                        for record in records:
                                seqid = record.id
                                seq = str(record.seq)
                                if seqid in sigpDict:
                                        coord = sigpDict[seqid]
                                        seq = seq[0:coord[0]-1] + ('x' * (coord[1] + 1 - coord[0])) + seq[coord[1]:]        # -1 to coord[0] to make it 0-based; +1 to coord[1] since a domain range of 1-1 still has a length of 1;
                                        fileOut.write('>' + seqid + '\n' + seq + '\n')
                                else:
                                        fileOut.write('>' + seqid + '\n' + seq + '\n')

        ### SEG AND COILS (Loads in signalP masked file)
        def runsegandcoils(args, outdir, basename):
                import os, subprocess, platform, io
                from Bio import SeqIO
                ### RUN SEG
                print('Masking LCRs from sequences...')
                segout = os.path.join(os.getcwd(), outdir, basename + '_seg.fasta')
                if not os.path.isfile(segout):
                        cmd = os.path.join(args['segdir'], 'seg') + ' "' + os.path.join(os.getcwd(), outdir, basename + '_signalp.fasta') + '" -x > ' + '"' + segout + '"'
                        run_seg = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                        run_seg.wait()
                        segErr = run_seg.stderr.read().decode("utf-8")
                        if segErr != '':
                                raise Exception('SEG error text below\n' + segErr)
                ### RUN COILS
                print('Masking coils from sequences...')
                # Format cmd
                coils_outname = os.path.join(os.getcwd(), outdir, basename + '_segcoils.fasta')
                cmd = '"' + os.path.join(args['python2dir'], 'python') + '" "' + os.path.join(args['coilsdir'], 'psCoils.py') + '" -f "' + segout + '"'
                #run_coils = subprocess.check_output(cmd, shell = True, stderr = subprocess.PIPE)
                run_coils = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                coilsout, coilserr = run_coils.communicate()
                if coilserr.decode("utf-8") != '':
                        raise Exception('COILS error text below\n' + coilserr.decode("utf-8"))         
                # Parse coils results
                output = coilsout.decode("utf-8").split(' Pos A Hep Score   Prob        Gcc         Gg        Pred (Loop=L Coiledcoil=C)')
                parsed_coils = {}   # Need to use a dictionary since coils doesn't appear to process sequences in the order that they occur in a fasta file?
                for pred in output:
                        if pred == '': continue
                        tempcoil = ''
                        tempseq = ''
                        for row in pred.replace('\r', '').split('\n'):
                                if row.endswith('L') or row.endswith('C'):
                                        temprow = row.split()
                                        tempseq += temprow[1]
                                        tempcoil += temprow[-1]
                        parsed_coils[tempseq] = tempcoil
                # Mask coiled coil domains
                coils_seqids = []
                coils_seqs = []
                records = SeqIO.parse(open(segout, 'rU'), 'fasta')
                for record in records:
                        seqid = record.id
                        seq = str(record.seq)
                        coilpred = parsed_coils[seq]
                        for i in range(len(coilpred)):
                                if coilpred[i] == 'C':
                                        seq = seq[:i] + 'x' + seq[i+1:]
                        coils_seqids.append('>' + seqid)
                        coils_seqs.append(seq)
                # Fix up capital X's introduced by seg
                coils_output = []
                for i in range(len(coils_seqs)):
                        coils_seqs[i] = coils_seqs[i].replace('X', 'x')
                        coils_output.append(coils_seqids[i] + '\n' + coils_seqs[i])
                # Create combined seg + COILS output file
                outfile = open(coils_outname, 'w')
                outfile.write('\n'.join(coils_output))
                outfile.close()

        def cleanseqs(args, outdir, basename):
                import os
                from Bio import SeqIO
                # Load in sequence file
                preclean_file = os.path.join(os.getcwd(), outdir, basename + '_segcoils.fasta')
                records = SeqIO.parse(open(preclean_file, 'rU'), 'fasta')
                # Create output file
                outname = os.path.join(os.getcwd(), outdir, basename + '_clean.fasta')
                outfile = open(outname, 'w')
                # Begin cleaning file
                print('Cleaning up fasta file before running PSI-BLAST...')
                outputText = []
                ongoingCount = 0
                for record in records:
                        # Processing steps
                        seqid = record.id
                        seq = str(record.seq)
                        aaSeq = seq.replace('x', '')
                        if len(aaSeq) >= int(args['cleanAA']):
                                outputText.append('>' + seqid + '\n' + seq)
                        ongoingCount += 1
                        # Make backup files to reduce memory usage
                        if ongoingCount%100000 == 0 and os.path.isfile(outname) == False:
                                output = open(outname, 'w')
                                output.write('\n'.join(outputText))
                                output.close()
                                print('Backup made after ' + str(ongoingCount) + ' sequences cleaned.')
                                outputText = ''
                        elif ongoingCount%100000 == 0 and os.path.isfile(outname) == True:
                                output = open(outname, 'a')
                                output.write('\n'.join(outputText))
                                output.close()
                                print('Backup made after ' + str(ongoingCount) + ' sequences cleaned...')
                                outputText = ''

                # Dump the last few results after the script has finished, or create the output if there were less than 10,000 sequences
                if os.path.isfile(outname) == False:
                        output = open(outname, 'w')
                        output.write('\n'.join(outputText))
                        output.close()
                        print('Final save made after ' + str(ongoingCount) + ' sequences cleaned...')
                elif os.path.isfile(outname) == True:
                        output = open(outname, 'a')
                        output.write('\n'.join(outputText))
                        output.close()
                        print('Final save made after ' + str(ongoingCount) + ' sequences cleaned...')

        ######## FUNCTIONS RELATING TO MMseqs2
        def makemms2db(args, outdir, basename):
                import os, subprocess, platform
                # Format command
                fasta = os.path.join(os.getcwd(), outdir, basename + '_clean.fasta')
                db = os.path.join(os.getcwd(), outdir, basename + '_mmseqs2DB')
                tmpdir = os.path.join(os.getcwd(), outdir, 'mms2tmp')
                if platform.system() == 'Windows':
                        cmd = os.path.join(args['mmseqs2dir'], 'mmseqs.exe') + ' createdb "' + fasta + '" "' + db + '"'
                else:
                        cmd = os.path.join(args['mmseqs2dir'], 'mmseqs') + ' createdb "' + fasta + '" "' + db + '"'
                print('Making MMseqs2 db...')
                run_makedb = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                makedbout, makedberr = run_makedb.communicate()
                if makedberr.decode("utf-8") != '':
                        raise Exception('Make MMseqs2 db error text below\n' + makedberr.decode("utf-8"))
                # Run further indexing
                if platform.system() == 'Windows':
                        #cmd = os.path.join(args['mmseqs2dir'], 'mmseqs.exe') + ' createindex "' + db + '" "' + db + '" "' + tmpdir + '" --threads ' + str(args['threads'])
                        cmd = os.path.join(args['mmseqs2dir'], 'mmseqs.exe') + ' createindex "' + db + '" "' + tmpdir + '" --threads ' + str(args['threads'])
                else:
                        #cmd = os.path.join(args['mmseqs2dir'], 'mmseqs') + ' createindex "' + db + '" "' + db + '" "' + tmpdir + '" --threads ' + str(args['threads'])
                        cmd = os.path.join(args['mmseqs2dir'], 'mmseqs') + ' createindex "' + db + '" "' + tmpdir + '" --threads ' + str(args['threads'])
                print('Indexing MMseqs2 db...')
                run_index = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                indexout, indexerr = run_index.communicate()
                if indexerr.decode("utf-8") != '':
                        raise Exception('Indexing MMseqs2 db error text below\n' + indexerr.decode("utf-8"))

        def runmms2(args, outdir, basename):
                import os, subprocess, platform
                # Format command
                if platform.system() == 'Windows':
                        db = os.path.join(outdir, basename + '_mmseqs2DB')          # MMS2 behaves a bit weirdly with Cygwin at this step. In order to get it to work, we need to use . to represent the cwd.
                        tmpdir = os.path.join(outdir, 'mms2tmp')
                        outprefix = os.path.join(outdir, basename + '_mmseqs2SEARCH')
                        evalue = args.mms2eval
                        cmd = os.path.join(args['mmseqs2dir'], 'mmseqs.exe') + ' search .\\' + db + ' .\\' + db + ' .\\' + outprefix + ' .\\' + tmpdir + ' --threads ' + str(args['threads']) + ' --num-iterations 4 -s 8'  # I considered adding the ability to change these values as an option, but MMSeqs2 is so fast that I don't think it's worth trading off accuracy for some speed improvement since it runs counter to the intent of this program to find novel domains not previously identified.
                else:
                        db = os.path.join(os.getcwd(), outdir, basename + '_mmseqs2DB')
                        tmpdir = os.path.join(os.getcwd(), outdir, 'mms2tmp')
                        outprefix = os.path.join(os.getcwd(), outdir, basename + '_mmseqs2SEARCH')
                        cmd = os.path.join(args['mmseqs2dir'], 'mmseqs') + ' search "' + db + '" "' + db + '" "' + outprefix + '" "' + tmpdir + '" --threads ' + str(args['threads']) + ' --num-iterations 4 -s 8'
                print('Running MMseqs2 profile iteration...')
                run_mms2 = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                mms2out, mms2err = run_mms2.communicate()
                if mms2err.decode("utf-8") != '':
                        raise Exception('MMseqs2 profile iteration error text below\n' + mms2err.decode("utf-8"))
                # Create tab-delim BLAST-like output
                if platform.system() == 'Windows':
                        querydb = os.path.join(os.getcwd(), outdir, basename + '_mmseqs2DB')
                        searchdb = os.path.join(os.getcwd(), outdir, basename + '_mmseqs2SEARCH')
                        tmpdir = os.path.join(os.getcwd(), outdir, 'mms2tmp')
                        #cmd = os.path.join(args['mmseqs2dir'], 'mmseqs.exe') + ' convertalis ".\\' + db + '" ".\\' + db + '" ".\\' + outprefix + '.m8" ".\\' + tmpdir + '" --threads ' + str(args['threads'])
                        cmd = os.path.join(args['mmseqs2dir'], 'mmseqs.exe') + ' convertalis "' + querydb + '" "' + querydb + '" "' + searchdb + '" "' + searchdb + '.m8" "' + tmpdir + '" --threads ' + str(args['threads'])
                else:
                        fasta = os.path.join(os.getcwd(), outdir, basename + '_clean.fasta')
                        db = os.path.join(os.getcwd(), outdir, basename + '_mmseqs2DB')
                        tmpdir = os.path.join(os.getcwd(), outdir, 'mms2tmp')
                        cmd = os.path.join(args['mmseqs2dir'], 'mmseqs') + ' convertalis "' + querydb + '" "' + querydb + '" "' + searchdb + '" "' + searchdb + '.m8" "' + tmpdir + '" --threads ' + str(args['threads'])
                print('Generating MMseqs2 tabular output...')
                run_mms2 = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                mms2out, mms2err = run_mms2.communicate()
                if mms2err.decode("utf-8") != '':
                        raise Exception('MMseqs2 tabular output generation error text below\n' + mms2err.decode("utf-8"))

        def parsemms2_peaks(args, outdir, basename):
                ### TO-DO: PARSE PBLAST FILE BEFORE PLATEAU ALGORITHM TO INCORPORATE QUERY AND HIT RESULTS ###
                import os, time, math
                import numpy as np
                from Bio import SeqIO
                from itertools import groupby
                from .peakdetect import peakdetect
                def plateau_extens(args, plateaus, values, arbitrary1, arbitrary2):
                        ## We can handle plateau extension without causing incorrect overlap by checking for the first of two occurrences. 1: we find a point where the length cutoff becomes enforced, or 2: the coverage starts to increase again, which means we're heading towards another peak (which wasn't collapsed into this plateau).
                        for i in range(len(plateaus)):
                                cutoff1 = math.ceil(values[i] * arbitrary1)
                                cutoff2 = math.ceil(values[i] * arbitrary2)
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
                                        # Low coverage check
                                        elif ongoingCount <= args['cleanAA']:
                                                if indexCov >= cutoff1:
                                                        newStart = x
                                                        continue
                                                else:
                                                        break
                                        else:
                                                if indexCov >= cutoff2:         # We start getting more strict now that we're beyond our expected domain region length
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
                                        # Low coverage check
                                        elif ongoingCount <= args['cleanAA']:
                                                if indexCov >= cutoff1:
                                                        newEnd = x
                                                        continue
                                                else:
                                                        break
                                        else:
                                                if indexCov >= cutoff2:         # We start getting more strict now that we're beyond our expected domain region length
                                                        newEnd = x
                                                        continue
                                                else:
                                                        break
                                plateaus[i][1] = newEnd
                        return plateaus
                                        
                # Prep for later output generation
                infile = os.path.join(os.getcwd(), outdir, basename + '_cdhit.fasta')
                #records = SeqIO.parse(open(infile, 'rU'), 'fasta')
                outname = os.path.join(os.getcwd(), outdir, basename + '_unclustered_domains.fasta')
                # Get the sequence lengths for each query
                seqArrays = {}
                fileIn = os.path.join(os.getcwd(), outdir, basename + '_cdhit.fasta')
                records = SeqIO.parse(open(fileIn, 'rU'), 'fasta')
                for record in records:
                        seqId = record.id
                        seqLen = len(str(record.seq))
                        seqArrays[seqId] = [0]*seqLen         # This can be converted into a numpy array later
                # Load in the PBLAST file as a groupby iterator
                mms2name = os.path.join(os.getcwd(), outdir, basename + '_mmseqs2SEARCH.m8')
                grouper = lambda x: x.split('\t')[0]
                domDict = {}        # This will hold ranges of potential domains associated with sequence IDs as key
                with open(mms2name, 'r') as mms2file, open(outname, 'w') as fileOut:
                        # Parse pblast file to pull out coverage per position arrays
                        for key, group in groupby(mms2file, grouper):
                                if key == '\n':
                                        continue
                                for entry in group:
                                        # Obtain data from this mms2 hit
                                        line = entry.split('\t')
                                        qname = line[0]
                                        hname = line[1]
                                        if qname == hname:        # i.e., skip self-hits
                                                continue
                                        qstart = int(line[6])
                                        hstart = int(line[8])
                                        qend = int(line[7])
                                        hend = int(line[9])
                                        # Update the numpy array-like list
                                        for i in range(qstart-1, qend):
                                                seqArrays[qname][i] += 1
                                        for i in range(hstart-1, hend):
                                                seqArrays[hname][i] += 1
                        # Handle plateaus
                        for key in seqArrays.keys():
                                # Check if we got any hits
                                if sum(seqArrays[key]) == 0:                        # This means we only had self-hits
                                        continue
                                elif len(set(seqArrays[key])) == 1:         # This means we only had completely overlapping hits against this sequence  ## CHECK THE OUTCOME OF THIS ##
                                        continue
                                # Convert to an array
                                array = np.array(seqArrays[key])
                                # Find peak indices
                                maxindices, minindices = peakdetect.peakdet(array,0.5)
                                # Get plateau regions
                                plateaus = []
                                values = []
                                arbitrary1 = 0.50   # ADD THESE AS ARGUMENTS IF ACCEPTED INTO FINAL SCRIPT VERSION
                                arbitrary2 = 0.75   # ADD THESE AS ARGUMENTS IF ACCEPTED INTO FINAL SCRIPT VERSION
                                for maximum in maxindices:
                                        index = maximum[0]
                                        value = maximum[1]
                                        # Look forward          [the peakdetect.peakdet values are always at the start of the plateau, so we don't need to look back, we just need to look forward to find where the plateau ends]
                                        for i in range(index, len(array)):
                                                if array[i] == value:
                                                        continue
                                                else:
                                                        plat = [index,i-1]         # This is 0-indexed, and we -1 since we want the previous i value
                                                        break
                                        plateaus.append(plat)
                                        values.append(value)
                                if len(plateaus) == 1:                          ### ADD EXTENSION: MAKE PLATEAU EXTENSION A FUNCTION ###
                                        plateaus = plateau_extens(args, plateaus, values, arbitrary1, arbitrary2)
                                        domDict[key] = plateaus
                                        continue
                                # Chain plateaus together
                                for i in range(len(plateaus)-1):                # S=start,E=end,L=length
                                        depressS = plateaus[i][1] + 1
                                        depressE = plateaus[i+1][0] -1                  # +1/-1 to start/end respectively  since the plateau ranges are the actual regions of overlap, i.e., 1->3 means 1, 2, and 3. Thus, 4 is the first character that does not overlap.
                                        depressL = depressE - depressS + 1          # +1 here to get the proper length of the gap region (i.e., a gap of 1 amino acid might look like 3->3, 3-3 = 0, so we need to +1
                                        depressR = array[depressS:depressE + 1]
                                        if depressL <= args['cleanAA']:
                                                # Use arbitrary value 1 to determine if chaining is correct
                                                maxVal = max([values[i], values[i+1]])  # We determine chaining based on the highest coverage plateau for a few reasons. Most importantly, the higher coverage plateau is most likely to be a real domain, so we should prioritise it and make sure it doesn't join to lower coverage regions unless they meet our arbitrary chaining limits
                                                cutoff = math.ceil(maxVal * arbitrary1) # We round up to handle low numbers properly. For example, 2*0.5=1.0. We don't need to round this, and 1.0 is good here if depressL<= cleanAA. However, 2*0.75=1.5. Rounding down would be 1, but this would mean we chain probably separate domains together incorrectly. Thus, by rounding up, we are more strict.
                                                belowCutoff = np.where(depressR < cutoff)
                                                if len(belowCutoff[0]) > 0:
                                                        continue                                                        # We don't chain these two plateaus together if we have any positions with less coverage than our cutoff/have regions with no coverage
                                                else:
                                                        plateaus[i] = [plateaus[i][0], plateaus[i+1][1]]
                                                        plateaus[i+1] = [plateaus[i][0], plateaus[i+1][1]]
                                                        values[i] = maxVal
                                                        values[i+1] = maxVal                                # We want to prevent extending a plateau significantly beyond where it should like what could happen if we don't update the value here
                                        else:
                                                # Use arbitrary value 2 to determine if chaining is correct [we treat gaps longer than the expected minimum domain length more strictly]
                                                maxVal = max([values[i], values[i+1]])  # We determine chaining based on the highest coverage plateau for a few reasons. Most importantly, the higher coverage plateau is most likely to be a real domain, so we should prioritise it and make sure it doesn't join to lower coverage regions unless they meet our arbitrary chaining limits
                                                cutoff = math.ceil(maxVal * arbitrary2)
                                                belowCutoff = np.where(depressR < cutoff)
                                                if len(belowCutoff[0]) > 0:
                                                        continue                                                        # We don't chain these two plateaus together if we have any positions with less coverage than our cutoff
                                                else:
                                                        plateaus[i] = [plateaus[i][0], plateaus[i+1][1]]
                                                        plateaus[i+1] = [plateaus[i][0], plateaus[i+1][1]]
                                                        values[i] = maxVal
                                                        values[i+1] = maxVal
                                # Clean up the plateaus
                                while True:
                                        if len(plateaus) == 1:  # Exit condition if we reduced the number of plateaus to 1, thus meaning there cannot be overlaps
                                                break
                                        overlap = 'n'
                                        for y in range(len(plateaus)-1):
                                                if plateaus[y+1][0] > plateaus[y][1]:
                                                        continue
                                                else:
                                                        plateaus[y] = [plateaus[y][0], plateaus[y+1][1]]
                                                        del plateaus[y+1]
                                                        del values[y+1]
                                                        overlap = 'y'
                                                        break
                                        if overlap == 'n':          # Exit condition if we make it through the 'for y' loop without encountering any overlaps
                                                break
                                # Extend plateaus following slightly modified chaining rules
                                plateaus = plateau_extens(args, plateaus, values, arbitrary1, arbitrary2)
                                # Add results to our domDict for later output generation
                                domDict[key] = plateaus

                        # Format output fasta file
                        records = SeqIO.parse(open(infile, 'rU'), 'fasta')
                        for record in records:
                                seqid = record.id
                                if seqid in domDict:
                                        seq = str(record.seq)
                                        ranges = domDict[seqid]
                                        for i in range(len(ranges)):
                                                tmpDomain = seq[ranges[i][0]:ranges[i][1]+1]
                                                fileOut.write('>' + seqid + '_Domain_' + str(i+1) + '_' + str(ranges[i][0]+1) + '-' + str(ranges[i][1]+1) + '\n' + tmpDomain + '\n')

        def parsemms2_nccheck(args, outdir, basename):
                import os
                from Bio import SeqIO
                from itertools import groupby
                # Load in file as a groupby iterator
                print('Parsing MMseq2 output (N/C-check protocol)...')
                mms2_file = os.path.join(os.getcwd(), outdir, basename + '_mmseqs2SEARCH.m8')
                grouper = lambda x: x.split('\t')[0]
                # Pull out potential domain ranges
                domDict = {}
                with open(mms2_file, 'r') as fileIn:
                        for key, group in groupby(fileIn, grouper):
                                for entry in group:
                                        if entry == '\n':
                                                continue
                                        line = entry.split('\t')
                                        #### NEW ADDITION: Allow for the discovery of internal repeats. Previously we skipped all self-hits. Now, we look to see if there is overlap. If there is less than 5% or 5 AA overlap (to account for shorter domains), we accept it as a hit ###
                                        if line[0] == line[1]:
                                                #print('Active 1')
                                                range1 = set(range(int(line[6]), int(line[7])))
                                                range2 = set(range(int(line[8]), int(line[9])))
                                                if not (len(range1-range2) / len(range1) > 0.95 or len(range1) - len(range1-range2) < 5):
                                                        #print('Active 2')
                                                        continue
                                        if line[0] == line[1]:
                                                print('Active?')
                                        ### NEW ADDITION ###
                                        else:
                                                if not int(line[7]) + 1 - int(line[6]) < args['cleanAA']:                 # Added in later: This will mean we only find putative globular domains. Another function may be made to find regions specifically shorter than the argument length which are putative motifs.
                                                        if line[0] in domDict:
                                                                domDict[line[0]].append((int(line[6]), int(line[7])))
                                                        else:
                                                                domDict[line[0]] = [(int(line[6]), int(line[7]))]
                                                        if line[1] in domDict:
                                                                domDict[line[1]].append((int(line[8]), int(line[9])))
                                                        else:
                                                                domDict[line[1]] = [(int(line[8]), int(line[9]))]
                if len(domDict) == 0:
                        print('No potential novel domain regions were found from PSI-BLAST. Program end.')
                        quit()
                # Split overlapping domain ranges
                for key, value in domDict.items():
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
                                                if n_value[1] - n_value[0] >= args['cleanAA']:
                                                        # Remove N overhang from current value
                                                        value[i] = (n_value[1] + 1, value[i][1])
                                                        # Insert N overhang as new domain
                                                        value.insert(i, n_value)
                                                        # Check if the trimmed current value has been shortened too much
                                                        if value[i+1][1] - value[i+1][0] < args['cleanAA']:
                                                                del value[i+1]
                                                                continue        # Only need to continue here if we insert and delete the old val, since it doesn't disrupt our 'for i' loop
                                                        break                   # We break here if we insert and don't delete the old value, since it will disrupt our 'for i' loop
                                                        
                                                else:
                                                        # Remove N overhang from current value
                                                        value[i] = (value[i+1][0], value[i][1])
                                                        # Check if the trimmed current value has been shortened too much
                                                        if value[i][1] - value[i][0] < args['cleanAA']:
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
                                        if c_value[1] - c_value[0] >= args['cleanAA']:
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
                                                                        
                # Format output fasta file
                infile = os.path.join(os.getcwd(), outdir, basename + '_cdhit.fasta')
                records = SeqIO.parse(open(infile, 'rU'), 'fasta')
                outfile = open(os.path.join(os.getcwd(), outdir, basename + '_unclustered_domains.fasta'), 'w')
                for record in records:
                        seqid = record.id
                        if seqid in domDict:
                                seq = str(record.seq)
                                ranges = domDict[seqid]
                                for i in range(len(ranges)):
                                        tmpDomain = seq[ranges[i][0]-1:ranges[i][1]]
                                        outfile.write('>' + seqid + '_Domain_' + str(i+1) + '_' + str(ranges[i][0]) + '-' + str(ranges[i][1]) + '\n' + tmpDomain + '\n')
                outfile.close()
