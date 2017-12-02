class domfind:
    ### CD-HIT
    def runcdhit(args, outdir, raw_fasta, basename):
        import os, subprocess
        # Format command
        cdhit_command = os.path.join(args['cdhitdir'], 'cd-hit') + ' -i ' + raw_fasta + ' -o ' + os.path.join(os.getcwd(), outdir, basename + '_cdhit.fasta') + ' -c ' + str(args['cdc']) + ' -n ' + str(args['cdn']) + ' -G ' + str(args['cdg']) + ' -aS ' + str(args['cdas']) + ' -aL ' + str(args['cdal']) + ' -M ' + str(args['cdm']) + ' -T ' + str(args['threads'])
        # Run CD-HIT
        print('Running CD-HIT...')
        run_cdhit = subprocess.Popen(cdhit_command, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        cdout, cderr = run_cdhit.communicate()
        if cderr.decode("utf-8") != '':
            raise Exception('CD-HIT Error text below' + str(cderr.decode("utf-8")))

    def chunk_fasta(args, outdir, basename):
        import os
        from Bio import SeqIO
        print('Chunking ' + basename + '_cdhit.fasta file to allow multi-threading operations...')
        with open(os.path.join(os.getcwd(), outdir, basename + '.fasta'), 'r') as fasta_file:
            num_seqs = 0
            for line in fasta_file:
                if line.startswith('>'):
                    num_seqs += 1
        # Check user argument to make sure it is sensible
        if num_seqs < int(args['threads']):
            print('You are using more threads than the number of sequences in your fasta file (' + str(num_seqs) + '). This isn\'t productive, and this program can\'t handle empty files. Specify a lower number of threads.')
            quit()
        # Find out where we are chunking the file
        chunk_size = int(num_seqs / int(args['threads']))
        cutoff_points = []
        for i in range(1, int(args['threads'])):
            cutoff_points.append(chunk_size * i)
        # Perform the chunking
        cdhit_file = os.path.join(os.getcwd(), outdir, basename + '.fasta')
        records = SeqIO.parse(open(cdhit_file, 'rU'), 'fasta')
        outfile = open(os.path.join(os.getcwd(), outdir, basename + '_chunk1.fasta'), 'w')
        ongoingCount = 1
        ongoingThread = 1
        for record in records:
            if ongoingCount not in cutoff_points:
                outfile.write('>' + record.id + '\n' + str(record.seq) + '\n')
                ongoingCount += 1
            else:
                outfile.write('>' + record.id + '\n' + str(record.seq) + '\n')
                outfile.close()
                ongoingThread += 1
                ongoingCount += 1
                outfile = open(os.path.join(os.getcwd(), outdir, basename + '_chunk' + str(ongoingThread)  + '.fasta'), 'w')
        outfile.close()
        
    ### HMMER3
    def runhmmer3(args, outdir, basename):
        import os, subprocess
        # Check if we need to run hmmpress
        if not os.path.isfile(args['hmmdb'] + '.h3f') and not os.path.isfile(args['hmmdb'] + '.h3i') and not os.path.isfile(args['hmmdb'] + '.h3m') and not os.path.isfile(args['hmmdb'] + '.h3p'):
            print('Running hmmpress to format your database...')
            cmd = os.path.join(args['hmmer3dir'], 'hmmpress') + ' -f "' + args['hmmdb'] + '"'
            run_hmmpress = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
            hmmout, hmmerr = run_hmmpress.communicate()
            if hmmerr.decode("utf-8") != '':
                raise Exception('hmmpress error text below' + str(hmmerr.decode("utf-8")))
        # Run HMMER3
        cmd = os.path.join(args['hmmer3dir'], 'hmmsearch') + ' --cpu ' + str(args['threads']) + ' -E ' + args['hmmeval'] + ' --domtblout ' + os.path.join(os.getcwd(), outdir, basename + '_hmmer.results') + ' "' + args['hmmdb'] + '" "' + os.path.join(os.getcwd(), outdir, basename + '_cdhit.fasta') + '"'
        print('Running hmmsearch to detect known domains...')
        run_hmmer3 = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        hmmout, hmmerr = run_hmmer3.communicate()
        if hmmerr.decode("utf-8") != '':
            raise Exception('hmmsearch error text below' + str(hmmerr.decode("utf-8")))
        
    def hmmerparse(args, evalarg, outdir, basename):
        # Script modified from one written by Andrzej Zielezinski (http://www.staff.amu.edu.pl/~andrzejz/)
        import os
        print('Parsing hmmsearch output...')
        # Load in file
        with open(os.path.join(os.getcwd(), outdir, basename + '_hmmer.results'), 'r') as fh:
            domdict = {}
            # Loop through file
            for line in fh:
                if line.startswith('#'): continue  # Skip header and footer of a domtblout file.
                if line == '' or line == '\n': continue     # Skip blank lines, shouldn't exist, but can't hurt
                sl = line.split()
                evalue = float(sl[12])
                if evalue > float(evalarg):
                    continue
                pid = sl[0]
                if sl[3].startswith('cath'):
                    did = sl[3]
                else:
                    did = os.path.basename(sl[3])           # Some of the databases will have the full path to the domain ID, so we want to get rid of those. We handle cath especially since cath's domain IDs don't exhibit this behaviour but do incorporate a '/' in their ID
                # Optional skipping of CATH/SUPERFAMILY databases
                if args['skip'] != 'noskip':
                    if args['skip'] == 'cath' and sl[3].startswith('cath'):
                        continue
                    elif args['skip'] == 'superfamily':
                        try:
                            int(did)        # SUPERFAMILY is the only database that has purely integer domain IDs
                            continue
                        except ValueError:
                            doNothing = ''
                    elif args['skip'] == 'both':
                        if sl[3].startswith('cath'):
                            continue
                        try:
                            int(did)
                            continue
                        except ValueError:
                            doNothing = ''
                # End of optional skipping
                dstart = int(sl[17])
                dend = int(sl[18])
                # Add into domain dictionary    # We need to use a dictionary for later sorting since hmmsearch does not produce output that is ordered in the way we want to work with. hmmscan does, but it is SIGNIFICANTLY slower.
                if pid not in domdict:          # This unfortunately serves as the only aspect of this script that should pose potential problems for memory use. To my knowledge there is no way to re-order hmmsearch output into the same format as hmmscan that is both fast and memory efficient.
                    domdict[pid] = [[dstart, dend, did]]
                else:
                    domdict[pid].append([dstart, dend, did])
        # Produce output
        ongoingCount = 0
        outname = os.path.join(os.getcwd(), outdir, basename + '_hmmerParsed.results')
        outputText = []
        with open(outname, 'w') as outFile:
            for key, value in domdict.items():
                value.sort()
                for domregion in value:
                    outFile.write(key + '\t' + str(domregion[0]) + '\t' + str(domregion[1]) + '\t' + domregion[2] + '\n')
        
    def hmmercutter(args, outdir, basename):
        import os
        from Bio import SeqIO
        from itertools import groupby
        print('Masking known domains from sequences...')
        # Get relevant inputs
        fileName = os.path.join(os.getcwd(), outdir, basename + '_cdhit.fasta')
        outputFileName = os.path.join(os.getcwd(), outdir, basename + '_domCut.fasta')
        pfamName = os.path.join(os.getcwd(), outdir, basename + '_hmmerParsed.results')
        # Load in files
        sequences = SeqIO.parse(open(fileName, 'rU'), 'fasta')
        pfams = open(pfamName, 'r')
        # Load the pfam domtblout file as a groupby iterator and sort out overlapping domains
        grouper = lambda x: x.split('\t')[0]
        pfamDict = {}
        for key, group in groupby(pfams, grouper):
            coords = []
            for entry in group:
                line = entry.rstrip('\n').split('\t')
                coords.append((int(line[1]), int(line[2])))
            # Remove identical coordinates and re-sort the list
            coords = list(set(coords))
            coords.sort()
            # Join overlaps together into a single region
            overlapping = 'y'
            while True:
                if len(coords) == 1 or overlapping == 'n':
                    break
                for y in range(len(coords)-1):
                    if coords[y+1][0] > coords[y][1]+1 and y != len(coords)-2:     # Adding +1 to coords[y][1] means we'll collapse domain ranges like (100-297, 298-400) into a single value
                        continue
                    elif coords[y+1][0] <= coords[y][1]+1:
                        highest = max(coords[y][1], coords[y+1][1])
                        coords[y] = [coords[y][0], highest] ### ERROR? ### WAS [coords[0][0], highest]
                        del coords[y+1]
                        break
                    else:                                                           # We need the y != check above since we need to set an exit condition when no more overlaps are present. The if/elif will always trigger depending on whether there is/is not an overlap UNLESS it's the second last entry and there is no overlap. In this case we finally reach this else clause, and we trigger an exit.
                        overlapping = 'n'
                        break
            # Retain results for subsequent cutting in a dictionary
            pfamDict[key] = coords

        # Remove domain regions from fasta
        with open(outputFileName, 'w') as outFile:
            for sequence in sequences:
                # Processing steps
                seqName = sequence.id
                if seqName not in pfamDict:
                    outFile.write('>' + seqName + '\n' + str(sequence.seq) + '\n')
                else:
                    # Mask the domain regions
                    currSeq = str(sequence.seq)
                    for coord in pfamDict[seqName]:
                        currSeq = currSeq[0:coord[0]-1] + ('x' * (coord[1] + 1 - coord[0])) + currSeq[coord[1]:]
                    outFile.write('>' + seqName + '\n' + currSeq + '\n')

    ### SIGNALP (Loads in CD-HIT results for running program, then modifies HMMER3 cut results)
    def runsignalp(args, outdir, basename):
        import os, subprocess, platform, threading
        from Bio import SeqIO
        def signalp_threading(args, outdir, basename, thread):
            # Format signalP script file for use on unix systems and with cygwin
                script_text = '"' + os.path.join(args['signalpdir'], 'signalp') + '" -t ' + args['signalporg'] + ' -f short -n "' + os.path.join(os.getcwd(), outdir, basename + '_signalp_chunk' + str(thread) + '.results"') + ' "' + os.path.join(os.getcwd(), outdir, basename + '_cdhit_chunk' + str(thread) + '.fasta"') # Loads in CD-HIT results here intentionally
                if platform.system() == 'Windows':
                    signalpscript_location = os.path.join(os.getcwd(), outdir, basename + '_signalpScript' + str(thread) + '.sh')
                    signalpscript = open(signalpscript_location, 'w')
                    signalpscript.write(script_text.replace('\\', '/'))
                    signalpscript.close()
                # Run signalP depending on operating system
                if platform.system() == 'Windows':
                    cmd = os.path.join(args['cygwindir'], 'bash') + ' -l -c ' + signalpscript_location.replace('\\', '/')
                    runsigP = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                    sigpout, sigperr = runsigP.communicate()
                    for line in sigperr.decode("utf-8").split('\n'):
                        if line.rstrip('\n') == '# No sequences predicted with a signal peptide':
                            with open(os.path.join(os.getcwd(), outdir, basename + '_signalp_chunk' + str(thread) + '.results'), 'w') as null_out:
                                null_out.write(line)
                            break
                        elif not 'is an unknown amino amino acid' in line and not line == '':
                            print(line + '<')
                            print('--')
                            raise Exception('SignalP error (thread ' + str(thread) + ') text below\n' + str(sigperr.decode("utf-8")))
                else:
                    runsigP = subprocess.Popen(script_text, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                    sigpout, sigperr = runsigP.communicate()
                    for line in sigperr.decode("utf-8").split('\n'):
                        if line.rstrip('\n') == '# No sequences predicted with a signal peptide':
                            with open(os.path.join(os.getcwd(), outdir, basename + '_signalp_chunk' + str(thread) + '.results'), 'w') as null_out:
                                null_out.write(line)
                            break
                        elif not 'is an unknown amino amino acid' in line and not line == '':
                            print(line + '<')
                            print('--')
                            raise Exception('SignalP error (thread ' + str(thread) + ') text below\n' + str(sigperr.decode("utf-8")))
        def signalp_1thread(args, outdir, basename):
            if not os.path.isfile(os.path.join(os.getcwd(), outdir, basename + '_signalp.results')):
                # Format signalP script file for use on unix systems and with cygwin
                script_text = '"' + os.path.join(args['signalpdir'], 'signalp') + '" -t ' + args['signalporg'] + ' -f short -n "' + os.path.join(os.getcwd(), outdir, basename + '_signalp.results"') + ' "' + os.path.join(os.getcwd(), outdir, basename + '_cdhit.fasta"') # Loads in CD-HIT results here intentionally
                if platform.system() == 'Windows':
                    signalpscript_location = os.path.join(os.getcwd(), outdir, basename + '_signalpScript.sh')
                    signalpscript = open(signalpscript_location, 'w')
                    signalpscript.write(script_text.replace('\\', '/'))
                    signalpscript.close()
                # Run signalP depending on operating system
                if platform.system() == 'Windows':
                    cmd = os.path.join(args['cygwindir'], 'bash') + ' -l -c ' + signalpscript_location.replace('\\', '/')
                    runsigP = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                    sigpout, sigperr = runsigP.communicate()
                    for line in sigperr.decode("utf-8").split('\n'):
                        if line.rstrip('\n') == '# No sequences predicted with a signal peptide':
                            with open(os.path.join(os.getcwd(), outdir, basename + '_signalp.results'), 'w') as null_out:
                                null_out.write(line)
                            break
                        elif not 'is an unknown amino amino acid' in line and not line == '':
                            print(line + '<')
                            print('--')
                            raise Exception('SignalP error text below\n' + str(sigperr.decode("utf-8")))
                else:
                    runsigP = subprocess.Popen(script_text, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                    sigpout, sigperr = runsigP.communicate()
                    for line in sigperr.decode("utf-8").split('\n'):
                        if line.rstrip('\n') == '# No sequences predicted with a signal peptide':
                            with open(os.path.join(os.getcwd(), outdir, basename + '_signalp.results'), 'w') as null_out:
                                null_out.write(line)
                            break
                        elif not 'is an unknown amino amino acid' in line and not line == '':
                            print(line + '<')
                            print('--')
                            raise Exception('SignalP error text below\n' + str(sigperr.decode("utf-8")))
    
        print('Masking signal peptides from sequences...')
        if not os.path.isfile(os.path.join(os.getcwd(), outdir, basename + '_signalp.results')):
            # Run signalP multi-threaded if enabled
            if int(args['threads']) > 1:
                processing_threads = []
                # Begin the loop
                for i in range(int(args['threads'])):
                    build = threading.Thread(target=signalp_threading, args=(args, outdir, basename, i+1))
                    processing_threads.append(build)
                    build.start()
                    print('........Initiated thread num ' + str(i+1) + ' for signalP operations...')
                # Wait for all threads to end.
                for process_thread in processing_threads:
                    process_thread.join()
                print('........SignalP completed...')
            # Run signalP signal-thread
            else:
                signalp_1thread(args, outdir, basename)
        # Join signalP output files if run multi-threaded
        if int(args['threads']) > 1:
            with open(os.path.join(os.getcwd(), outdir, basename + '_signalp.results'), 'w') as outFile:
                for i in range(int(args['threads'])):
                    sigp_chunk_file = open(os.path.join(os.getcwd(), outdir, basename + '_signalp_chunk' + str(i+1) + '.results'), 'r').read()
                    outFile.write(sigp_chunk_file)
        # Parse signalP output
        sigpred = {}
        with open(os.path.join(os.getcwd(), outdir, basename + '_signalp.results'), 'r') as sigp_results:
            for line in sigp_results:
                if line.startswith('#'):
                    continue
                sl = line.split('\t')
                sigpred[sl[0]] = [int(sl[3]), int(sl[4])]            
        # Mask signal peptides
        sigp_out = []
        records = SeqIO.parse(open(os.path.join(os.getcwd(), outdir, basename + '_domCut.fasta'), 'rU'), 'fasta')
        with open(os.path.join(os.getcwd(), outdir, basename + '_signalp.fasta'), 'w') as sigp_outfile:
            for record in records:
                seqid = record.id
                seq = str(record.seq)
                if seqid in sigpred:
                    sig_location = sigpred[seqid]
                    seq = seq[0:sig_location[0]-1] + ('x' * (sig_location[1] + 1 - sig_location[0])) + seq[sig_location[1]:]
                    sigp_outfile.write('>' + seqid + '\n' + seq + '\n')
                else:
                    sigp_outfile.write('>' + seqid + '\n' + seq + '\n')
    
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
        output = coilsout.decode("utf-8").split(' Pos A Hep Score   Prob    Gcc     Gg    Pred (Loop=L Coiledcoil=C)')
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
            db = os.path.join(outdir, basename + '_mmseqs2DB')      # MMS2 behaves a bit weirdly with Cygwin at this step. In order to get it to work, we need to use . to represent the cwd.
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
                ongoingCount = 0            # This measures how long our extension is so we can apply the two arbitrary values when appropriate
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
                        if indexCov >= cutoff2:     # We start getting more strict now that we're beyond our expected domain region length
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
                        if indexCov >= cutoff2:     # We start getting more strict now that we're beyond our expected domain region length
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
            seqArrays[seqId] = [0]*seqLen     # This can be converted into a numpy array later
        # Load in the PBLAST file as a groupby iterator
        mms2name = os.path.join(os.getcwd(), outdir, basename + '_mmseqs2SEARCH.m8')
        grouper = lambda x: x.split('\t')[0]
        domDict = {}    # This will hold ranges of potential domains associated with sequence IDs as key
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
                    if qname == hname:    # i.e., skip self-hits
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
                if sum(seqArrays[key]) == 0:            # This means we only had self-hits
                    continue
                elif len(set(seqArrays[key])) == 1:     # This means we only had completely overlapping hits against this sequence  ## CHECK THE OUTCOME OF THIS ##
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
                    # Look forward      [the peakdetect.peakdet values are always at the start of the plateau, so we don't need to look back, we just need to look forward to find where the plateau ends]
                    for i in range(index, len(array)):
                        if array[i] == value:
                            continue
                        else:
                            plat = [index,i-1]     # This is 0-indexed, and we -1 since we want the previous i value
                            break
                    plateaus.append(plat)
                    values.append(value)
                if len(plateaus) == 1:              ### ADD EXTENSION: MAKE PLATEAU EXTENSION A FUNCTION ###
                    plateaus = plateau_extens(args, plateaus, values, arbitrary1, arbitrary2)
                    domDict[key] = plateaus
                    continue
                # Chain plateaus together
                for i in range(len(plateaus)-1):        # S=start,E=end,L=length
                    depressS = plateaus[i][1] + 1
                    depressE = plateaus[i+1][0] -1          # +1/-1 to start/end respectively  since the plateau ranges are the actual regions of overlap, i.e., 1->3 means 1, 2, and 3. Thus, 4 is the first character that does not overlap.
                    depressL = depressE - depressS + 1      # +1 here to get the proper length of the gap region (i.e., a gap of 1 amino acid might look like 3->3, 3-3 = 0, so we need to +1
                    depressR = array[depressS:depressE + 1]
                    if depressL <= args['cleanAA']:
                        # Use arbitrary value 1 to determine if chaining is correct
                        maxVal = max([values[i], values[i+1]])  # We determine chaining based on the highest coverage plateau for a few reasons. Most importantly, the higher coverage plateau is most likely to be a real domain, so we should prioritise it and make sure it doesn't join to lower coverage regions unless they meet our arbitrary chaining limits
                        cutoff = math.ceil(maxVal * arbitrary1) # We round up to handle low numbers properly. For example, 2*0.5=1.0. We don't need to round this, and 1.0 is good here if depressL<= cleanAA. However, 2*0.75=1.5. Rounding down would be 1, but this would mean we chain probably separate domains together incorrectly. Thus, by rounding up, we are more strict.
                        belowCutoff = np.where(depressR < cutoff)
                        if len(belowCutoff[0]) > 0:
                            continue                            # We don't chain these two plateaus together if we have any positions with less coverage than our cutoff/have regions with no coverage
                        else:
                            plateaus[i] = [plateaus[i][0], plateaus[i+1][1]]
                            plateaus[i+1] = [plateaus[i][0], plateaus[i+1][1]]
                            values[i] = maxVal
                            values[i+1] = maxVal                # We want to prevent extending a plateau significantly beyond where it should like what could happen if we don't update the value here
                    else:
                        # Use arbitrary value 2 to determine if chaining is correct [we treat gaps longer than the expected minimum domain length more strictly]
                        maxVal = max([values[i], values[i+1]])  # We determine chaining based on the highest coverage plateau for a few reasons. Most importantly, the higher coverage plateau is most likely to be a real domain, so we should prioritise it and make sure it doesn't join to lower coverage regions unless they meet our arbitrary chaining limits
                        cutoff = math.ceil(maxVal * arbitrary2)
                        belowCutoff = np.where(depressR < cutoff)
                        if len(belowCutoff[0]) > 0:
                            continue                            # We don't chain these two plateaus together if we have any positions with less coverage than our cutoff
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
                    if overlap == 'n':      # Exit condition if we make it through the 'for y' loop without encountering any overlaps
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
