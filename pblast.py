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
        
    def hmmerparse(args, outdir, basename):
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
                if evalue > float(args['hmmeval']):
                    continue
                pid = sl[0]
                if sl[3].startswith('cath'):
                    did = sl[3]
                else:
                    did = os.path.basename(sl[3])           # Some of the databases will have the full path to the domain ID, so we want to get rid of those. We handle cath especially since cath's domain IDs don't exhibit this behaviour but do incorporate a '/' in their ID
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
                        coords[y] = [coords[0][0], highest]
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
                    print('Initiated thread num ' + str(i+1) + ' for signalP operations...')
                # Wait for all threads to end.
                for process_thread in processing_threads:
                    process_thread.join()
                print('SignalP completed...')
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
    
    ######## FUNCTIONS RELATING TO PSI-BLAST
    ### MAKE BLAST DB
    def makeblastdb(args, outdir, basename):
        import os, subprocess
        # Format command
        fasta = os.path.join(os.getcwd(), outdir, basename + '_clean.fasta')
        db = os.path.join(os.getcwd(), outdir, basename + '_blastdb')
        cmd = os.path.join(args['blastdir'], 'makeblastdb') + ' -in "' + fasta + '" -dbtype prot -out "' + db + '"'
        if args['parse'] != None:
            cmd += ' -parse_seqids'
        print('Making blast db...')
        run_makedb = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        makedbout, makedberr = run_makedb.communicate()
        if makedberr.decode("utf-8") != '':
            raise Exception('Make blast db error text below\n' + makedberr.decode("utf-8"))  

    ### RUN PSI-BLAST
    def runpblast(args, outdir, basename):
        import os, subprocess, threading
        from domfind import domfind
        def multithread_pblast(args, outdir, basename, thread):
            # Format command
            cmd = os.path.join(args['blastdir'], 'psiblast') + ' -num_iterations 10 -outfmt 6 -soft_masking y -evalue ' + str(args['psieval']) + ' -query ' + os.path.join(os.getcwd(), outdir, basename + '_clean_chunk' + str(thread) + '.fasta') + ' -db ' + os.path.join(os.getcwd(), outdir, basename + '_blastdb') + ' -out ' + os.path.join(os.getcwd(), outdir, basename + '_psireport_chunk' + str(thread) + '.results')
            #print('Running PSI-BLAST...')
            run_psi = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
            psiout, psierr = run_psi.communicate()
            if psierr.decode("utf-8") != '':            
                print('Check the last few lines of the PSI-BLAST stderr output to make sure no actual errors are present (note: warnings are expected)')    # Since PSI-BLAST always returns warnings, we can't use stderr output to find actual errors
                print(psierr.decode("utf-8")[-432:])
        def singlethread_pblast(args, outdir, basename):
            # Format command
            cmd = os.path.join(args['blastdir'], 'psiblast') + ' -num_iterations 10 -outfmt 6 -soft_masking y -evalue ' + str(args['psieval']) + ' -query ' + os.path.join(os.getcwd(), outdir, basename + '_clean.fasta') + ' -db ' + os.path.join(os.getcwd(), outdir, basename + '_blastdb') + ' -out ' + os.path.join(os.getcwd(), outdir, basename + '_psireport.results')
            print('Running PSI-BLAST...')
            run_psi = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
            psiout, psierr = run_psi.communicate()
            if psierr.decode("utf-8") != '':            
                print('Check the last few lines of the PSI-BLAST stderr output to make sure no actual errors are present (note: warnings are expected)')    # Since PSI-BLAST always returns warnings, we can't use stderr output to find actual errors
                print(psierr.decode("utf-8")[-432:])
        # Chunk the clean fasta file if running multithreaded
        if int(args['threads']) > 1:
            chunking = 'n'
            chunk_names = []
            for i in range(0, int(args['threads'])):
                chunk_names.append(os.path.join(os.getcwd(), outdir, basename + '_clean_chunk' + str(i+1) + '.fasta'))
            for name in chunk_names:
                if not os.path.isfile(name):
                    chunking = 'y'
                    break
            if os.path.isfile(os.path.join(os.getcwd(), outdir, basename + '_clean_chunk' + str(i+2) + '.fasta')):
                print('There are more "clean_chunk#" chunk files in the output directory than should exist given the number of threads provided in your argument')
                print('This probably means you have leftover files from a previous run which used more threads. Delete all of these "clean_chunk#" files to re-run with less threads, or alternatively provide a -threads argument equivalent to the number of chunk files in the output directory')
                quit()
            if chunking == 'y':
                domfind.chunk_fasta(args, outdir, basename + '_clean')
            processing_threads = []
            # Begin the loop
            for i in range(int(args['threads'])):
                build = threading.Thread(target=multithread_pblast, args=(args, outdir, basename, i+1))
                processing_threads.append(build)
                build.start()
                print('Initiated thread num ' + str(i+1) + ' for PSI-BLAST operations...')
            # Wait for all threads to end.
            for process_thread in processing_threads:
                process_thread.join()
            print('PSI-BLAST completed...')
        # Run signalP signal-thread
        else:
            singlethread_pblast(args, outdir, basename)        
        # Join output files if running multithreaded
        if int(args['threads']) > 1:
            outfile = open(os.path.join(os.getcwd(), outdir, basename + '_psireport.results'), 'w')
            for i in range(int(args['threads'])):
                psi_chunk_file = open(os.path.join(os.getcwd(), outdir, basename + '_psireport_chunk' + str(i+1) + '.results'), 'r').read()
                outfile.write(psi_chunk_file)
            outfile.close()

    ### PARSE PSI-BLAST
    def parsepblast_peaks(args, outdir, basename):
        import os, time
        import numpy as np
        from Bio import SeqIO
        from itertools import groupby
        # Get the sequence lengths for each query
        seqArrays = {}
        with open(os.path.join(os.getcwd(), outdir, basename + '_cdhit.fasta'), 'r') as fileIn:
            records = SeqIO.parse(open(fileIn, 'rU'), 'fasta')
            for record in records:
                seqId = record.id
                seqLen = len(str(record.seq))
                seqArrays[seqId] = [0]*seqLen     # This can be converted into a numpy array later
        # Load in the PBLAST file as a groupby iterator
        pblastname = os.path.join(os.getcwd(), outdir, basename + '_psireport.results')
        grouper = lambda x: x.split('\t')[0]
        with open(pblastname, 'r') as pblastfile:
            for key, group in groupby(pblastfile, grouper):
                group = list(group).reverse()   # This lets us parse the PBLAST results so that we can obtain the last iteration's result immediately and then skip redundant hits
                alreadyHit = []
                for entry in group:
                    if entry == '\n' or entry == 'Search has CONVERGED!\n':
                        continue
                    # Obtain data from this blast hit
                    line = entry.split('\t')
                    qname = line[0]
                    hname = line[1]
                    if qname == hname or hname in alreadyHit:    # i.e., skip self-hits, and skip redundant hits
                        continue
                    alreadyHit.append(hname)
                    start = int(line[6])
                    end = int(line[7])
                    # Update the numpy array-like list
                    for i in range(start-1, end):
                        seqArrays[key][i] += 1
                # Convert to an array
                array = np.array(seqArrays[key])
            
    def parsepblast_doms(args, outdir, basename):
        import os, time
        from Bio import SeqIO
        from itertools import groupby
        # Load in file as a groupby iterator
        print('Parsing PSI-BLAST output...')
        pblast_file = open(os.path.join(os.getcwd(), outdir, basename + '_psireport.results'), 'r')
        grouper = lambda x: x.split('\t')[0]
        # Pull out potential domain ranges
        domDict = {}
        for key, group in groupby(pblast_file, grouper):
            for entry in group:
                if entry == '\n' or entry == 'Search has CONVERGED!\n':
                    continue
                line = entry.split('\t')
                ### NEW ADDITION: Allow for the discovery of internal repeats. Previously we skipped all self-hits. Now, we look to see if there is overlap. If there is less than 5% or 5 AA overlap (to account for shorter domains), we accept it as a hit ###
                if line[0] == line[1]:
                    range1 = set(range(int(line[6]), int(line[7])))
                    range2 = set(range(int(line[8]), int(line[9])))
                    if not (len(range1-range2) / len(range1) > 0.95 or len(range1) - len(range1-range2) < 5):
                        continue
                ### NEW ADDITION ###
                else:
                    if not int(line[7]) + 1 - int(line[6]) < args['cleanAA']:         # Added in later: This will mean we only find putative globular domains. Another function may be made to find regions specifically shorter than the argument length which are putative motifs.
                        if line[0] in domDict:
                            domDict[line[0]].append((int(line[6]), int(line[7])))
                        else:
                            domDict[line[0]] = [(int(line[6]), int(line[7]))]
                        if line[1] in domDict:
                            domDict[line[1]].append((int(line[8]), int(line[9])))
                        else:
                            domDict[line[1]] = [(int(line[8]), int(line[9]))]
        starttime = time.time()
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
                    if value[i+1][0] > value[i][1]:                 # If these two domain regions do not overlap...
                        continue
                    #### N Check
                    if value[i][0] != value[i+1][0]:        # If these two domain regions do not start at the same position...
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
                                continue    # Only need to continue here if we insert and delete the old val, since it doesn't disrupt our 'for i' loop
                            break           # We break here if we insert and don't delete the old value, since it will disrupt our 'for i' loop
                            
                        else:
                            # Remove N overhang from current value
                            value[i] = (value[i+1][0], value[i][1])
                            # Check if the trimmed current value has been shortened too much
                            if value[i][1] - value[i][0] < args['cleanAA']:
                                del value[i]
                                break       # We break here if we delete the old value since it disrupts our 'for i' loop. Otherwise, we can continue straight into the C check process
                            
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
