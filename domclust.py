class domclust:
    def alfree_matrix(args, outdir, filename):
        import os
        from alfpy import word_pattern, word_vector, word_distance, ncd, word_d2
        from alfpy.utils import seqrecords, distmatrix
        from alfpy.utils.data import seqcontent
        from Bio import SeqIO
        print('Calculating alignment-free matrix...')
        # Read in unclustered domains file
        unclust_doms = open(os.path.join(os.getcwd(), outdir, filename))
        seq_records = seqrecords.read_fasta(unclust_doms)
        unclust_doms.close()
        seq_list = seq_records.seq_list
        length_list = seq_records.length_list
        id_list = seq_records.id_list
        # Optional reduction of protein alphabet
        if str(args['reduce']) == '15':
            murphy_15_tab = {"L":"L","V":"L","I":"L","M":"L","C":"C","A":"A","G":"G","S":"S","T":"T","P":"P","F":"F","Y":"F","W":"W","E":"E","D":"D","N":"N","Q":"Q","K":"K","R":"K","H":"H"}
            for i in range(len(seq_list)):
                newseq = ''
                for letter in seq_list[i]:
                    newseq += murphy_15_tab[letter]
                seq_list[i] = newseq
        elif str(args['reduce']) == '11':
            eleven_tab = seqcontent.get_reduced_alphabet('protein')
            for i in range(len(seq_list)):
                newseq = ''
                for letter in seq_list[i]:
                    if letter in eleven_tab:
                        newseq += eleven_tab[letter]
                    else:
                        newseq += letter
                seq_list[i] = newseq
        ### TEST
        # Compute distance matrix of choice
        if args['alf'] == 'google':
            # Calc for word_size = 2
            p = word_pattern.create(seq_list, word_size=2)
            counts = word_vector.Counts(length_list, p)
            dist = word_distance.Distance(counts, 'google')
            matrix2 = distmatrix.create(id_list, dist)
            # Now 1
            p = word_pattern.create(seq_list, word_size=1)
            counts = word_vector.Counts(length_list, p)
            dist = word_distance.Distance(counts, 'google')
            matrix1 = distmatrix.create(id_list, dist)
        elif args['alf'] == 'canberra':
            # 2
            p = word_pattern.create(seq_list, word_size=2)
            weightmodel = word_vector.WeightModel(seqcontent.get_weights('protein'))
            counts = word_vector.CountsWeight(length_list, p, weightmodel)
            dist = word_distance.Distance(counts, 'canberra')
            matrix2 = distmatrix.create(id_list, dist)
            # 1
            p = word_pattern.create(seq_list, word_size=1)
            weightmodel = word_vector.WeightModel(seqcontent.get_weights('protein'))
            counts = word_vector.CountsWeight(length_list, p, weightmodel)
            dist = word_distance.Distance(counts, 'canberra')
            matrix1 = distmatrix.create(id_list, dist)
        else:
            # 2
            p = word_pattern.create(seq_list, word_size=2)
            counts = word_vector.Counts(length_list, p)
            dist = word_distance.Distance(counts, 'braycurtis')
            matrix2 = distmatrix.create(id_list, dist)
            # 1
            p = word_pattern.create(seq_list, word_size=1)
            counts = word_vector.Counts(length_list, p)
            dist = word_distance.Distance(counts, 'braycurtis')
            matrix1 = distmatrix.create(id_list, dist)
        ### TESTING ###
        #import numpy as np
        #matrix2 = np.array(matrix2)
        #np.set_printoptions(threshold=np.inf)
        #matrix2.tofile('dataset.csv', sep=',', format="%s")
        #text = pyperclip.paste()
        #table = text.split(']\r\n [ ')
        #for i in range(len(table)):
        #    table[i] = table[i].replace('[', '')
        #    table[i] = table[i].replace(']', '')
        #for line in table:
        #    sl = line.split()
        #    out.append('\t'.join(sl))
        #pyperclip.copy('\r\n'.join(out))
        #print(id_list)
        #quit()
        ### TESTING ###
        # Return value
        return matrix1, matrix2, id_list

    def cluster_hdb(args, matrix1, matrix2, idlist, rejects_list):
        import os, hdbscan
        import numpy as np
        #import seaborn as sns
        #import pandas as pd
        print('Clustering putative domain groups...')
        group_dict = {}
        # Run clustering algorithm
        if args['leaf'] == 'y' and args['singleclust'] == 'y':
            clusterer = hdbscan.HDBSCAN(metric='precomputed', cluster_selection_method = 'leaf', min_cluster_size = int(args['minsize']), min_samples = int(args['minsample']), allow_single_cluster = True)
        elif args['leaf'] == 'y' and args['singleclust'] == 'n':
            clusterer = hdbscan.HDBSCAN(metric='precomputed', cluster_selection_method = 'leaf', min_cluster_size = int(args['minsize']), min_samples = int(args['minsample']))
        elif args['leaf'] == 'n' and args['singleclust'] == 'y':
            clusterer = hdbscan.HDBSCAN(metric='precomputed', min_cluster_size = int(args['minsize']), min_samples = int(args['minsample']), allow_single_cluster = True)   
        else:
            clusterer = hdbscan.HDBSCAN(metric='precomputed', min_cluster_size = int(args['minsize']), min_samples = int(args['minsample']))
        clusterer.fit(matrix2.data)     # We look at matrix2 (or word_size == 2) first since it should, theoretically, find more 'similar' groups better than a word_size of 1
        # Pull out domain groups
        clust_groups = clusterer.labels_
        # Sort groups
        for i in range(len(idlist)):
            if clust_groups[i] != -1:
                if clust_groups[i] not in group_dict:
                    group_dict[clust_groups[i]] = [idlist[i]]
                else:
                    group_dict[clust_groups[i]].append(idlist[i])
            else:
                rejects_list.append(idlist[i])
        #### TESTING ####
        # Look through the results of the matrix1 (or word_size == 1) for any groups not found for word_size == 2
        if args['leaf'] == 'y' and args['singleclust'] == 'y':
            clusterer = hdbscan.HDBSCAN(metric='precomputed', cluster_selection_method = 'leaf', min_cluster_size = int(args['minsize']), min_samples = int(args['minsample']), allow_single_cluster = True)
        elif args['leaf'] == 'y' and args['singleclust'] == 'n':
            clusterer = hdbscan.HDBSCAN(metric='precomputed', cluster_selection_method = 'leaf', min_cluster_size = int(args['minsize']), min_samples = int(args['minsample']))
        elif args['leaf'] == 'n' and args['singleclust'] == 'y':
            clusterer = hdbscan.HDBSCAN(metric='precomputed', min_cluster_size = int(args['minsize']), min_samples = int(args['minsample']), allow_single_cluster = True)   
        else:
            clusterer = hdbscan.HDBSCAN(metric='precomputed', min_cluster_size = int(args['minsize']), min_samples = int(args['minsample']))
        clusterer.fit(matrix1.data)
        new_groups = clusterer.labels_
        # Format new group into a dictionary
        new_dict = {}
        for i in range(len(idlist)):
            if new_groups[i] != -1:
                if new_groups[i] not in new_dict:
                    new_dict[new_groups[i]] = [idlist[i]]
                else:
                    new_dict[new_groups[i]].append(idlist[i])
        # Look for groups not found in previous clustering
        oldvals = []
        ongoingCount = 0
        oldlen = len(group_dict)
        for val in group_dict.values():
            oldvals += val
        skip = 'n'
        for key, value in new_dict.items():
            for val in value:
                if val in oldvals:
                    skip = 'y'
                    break
            if skip == 'y':
                skip = 'n'
                continue
            else:
                group_dict[ongoingCount + oldlen] = value
                ongoingCount += 1
                for val in value:
                    del rejects_list[rejects_list.index(val)]       # Now that these are no longer rejects, we delete them from our rejects_list
        return group_dict, rejects_list

    def mafft_align(args, outdir, basename, group_dict):
        import os, platform, shutil, threading
        from Bio import SeqIO
        from Bio.Align.Applications import MafftCommandline
        def run_mafft(args, outdir, basename, startNum, endNum, thread):
            for i in range(startNum, endNum):
                value = group_dict[i]
                tmp_text = []
                for seqid in value:
                    tmp_text.append('>' + records[seqid].id + '\n' + str(records[seqid].seq))
                tmp_text = '\n'.join(tmp_text)
                tmp_name = os.path.join(os.getcwd(), outdir, basename + '_tmpfile_thread' + str(thread) + '.fasta')
                tmp_file = open(tmp_name, 'w')
                tmp_file.write(tmp_text)
                tmp_file.close()
                # Run MAFFT
                if platform.system() == 'Windows':
                    mafft_cline = MafftCommandline(os.path.join(args['mafftdir'], 'mafft.bat'), input=tmp_name)
                    stdout, stderr = mafft_cline()
                else:
                    mafft_cline = MafftCommandline(os.path.join(args['mafftdir'], 'mafft'), input=tmp_name)
                    stdout, stderr = mafft_cline()
                if stdout == '':
                    raise Exception('MAFFT error text below' + str(stderr))
                # Process MAFFT output
                stdout = stdout.split('\n')
                while stdout[-1] == '\n' or stdout[-1] == '' or stdout[-1] == 'Terminate batch job (Y/N)?\n':   # Remove junk, sometimes MAFFT will have the 'Terminate ...' line
                    del stdout[-1]
                stdout = '\n'.join(stdout)
                # Create output alignment files
                out_name = os.path.join(align_out_dir, 'Domain_' + str(i) + '_align.fasta')
                mafft_out = open(out_name, 'w')
                mafft_out.write(stdout)
                mafft_out.close()
            
        print('Aligning domain groups...')
        # Read in fasta file
        unclust_doms = os.path.join(os.getcwd(), outdir, basename + '_unclustered_domains.fasta')
        records = SeqIO.to_dict(SeqIO.parse(open(unclust_doms, 'rU'), 'fasta')) # We should be able to read it in as a dictionary without being excessively memory intensive, as we are now only looking at putative domain regions
        # Clean out the temporary alignment file directory
        align_out_dir = os.path.join(os.getcwd(), outdir, 'tmp_alignments')
        if os.path.isdir(align_out_dir):
            for file in os.listdir(align_out_dir):
                file_path = os.path.join(align_out_dir, file)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                except Exception as e:
                    print(e)
        else:
            os.mkdir(align_out_dir)
        # Set up threading requirements
        dict_size = len(group_dict)
        chunk_size = int(dict_size / int(args['threads']))
        processing_threads = []
        # Begin the loop
        for i in range(int(args['threads'])):
            start = chunk_size * i
            if i+1 != int(args['threads']):
                end = chunk_size * (i+1)
            else:
                end = dict_size
            build = threading.Thread(target=run_mafft, args=(args, outdir, basename, start, end, i+1))
            processing_threads.append(build)
            build.start()
            print('Initiated thread num ' + str(i+1) + ' for MAFFT alignment...')

        # Wait for all threads to end.
        for process_thread in processing_threads:
            process_thread.join()
        print('MAFFT alignment completed')
            
    def cluster_hmms(args, outdir, basename):
        import os, subprocess
        print('Building HMMs for iterative searching...')
        # Get file details
        align_out_dir = os.path.join(os.getcwd(), outdir, 'tmp_alignments')
        msas = os.listdir(align_out_dir)
        hmms = []
        for msa in msas:
            # Format names
            out_name = '.'.join(msa.split('.')[0:-1]) + '.hmm'
            hmms.append(out_name)
            # Format cmd
            cmd = os.path.join(args['hmmer3dir'], 'hmmbuild') + ' "' + os.path.join(align_out_dir, out_name) + '" "' + os.path.join(align_out_dir, msa) + '"'
            # Run HMMBUILD
            run_hmmbuild = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
            hmmout, hmmerr = run_hmmbuild.communicate()
            if hmmerr.decode("utf-8") != '':
                raise Exception('hmmbuild error text below\n' + str(hmmerr.decode("utf-8")) + '\nMake sure that you define the -hmmer3dir argument if this directory is not in your PATH')
        # Concatenate HMMs
        outhmm_name = os.path.join(align_out_dir, 'dom_models.hmm')
        outhmm_file = open(outhmm_name, 'w')
        for hmm in hmms:
            outhmm_file.write(open(os.path.join(align_out_dir, hmm), 'r').read())
        outhmm_file.close()
        # Press HMMs
        cmd = os.path.join(args['hmmer3dir'], 'hmmpress') + ' -f "' + os.path.join(align_out_dir, outhmm_name) + '"'
        run_hmmpress = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        hmmout, hmmerr = run_hmmpress.communicate()
        if hmmerr.decode("utf-8") != '':
            raise Exception('hmmpress error text below' + str(hmmerr.decode("utf-8")))

    def hmmer3_doms(args, outdir, basename):
        import os, subprocess
        # Run HMMER3
        print('Running hmmsearch to grow domain clusters...')
        align_out_dir = os.path.join(os.getcwd(), outdir, 'tmp_alignments')
        outhmm_name = os.path.join(align_out_dir, 'dom_models.hmm')
        cmd = os.path.join(args['hmmer3dir'], 'hmmsearch') + ' --cpu ' + str(args['threads']) + ' -E ' + str(args['hmmevalnov']) + ' --domtblout ' + os.path.join(align_out_dir, 'tmp_hmmer.results') + ' "' + os.path.join(align_out_dir, outhmm_name) + '" "' + os.path.join(os.getcwd(), outdir, basename + '_cdhit.fasta') + '"'
        run_hmmer3 = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        hmmout, hmmerr = run_hmmer3.communicate()
        if hmmerr.decode("utf-8") != '':
            raise Exception('hmmsearch error text below' + str(hmmerr.decode("utf-8")))

    def parse_domtblout(domtbloutFile, evalue, outfile):
        import os
        from itertools import groupby
        # Define functions for later use
        def findMiddle(input_list):             # https://stackoverflow.com/questions/38130895/find-middle-of-a-list
            middle = float(len(input_list))/2
            if middle % 2 != 0:
                return [input_list[int(middle - .5)]]
            else:
                return [input_list[int(middle)], input_list[int(middle-1)]]

        def domain_ovl_resolver(ovlCutoff, seqHits):
                overlapping = 'y'
                while True:
                        if len(seqHits) == 1 or overlapping == 'n':
                                break
                        for y in range(len(seqHits)-1):
                                overlapping = 'y'
                                if seqHits[y+1][1] > seqHits[y][2] and y != len(seqHits)-2:
                                        continue
                                elif seqHits[y+1][1] == seqHits[y][2]:
                                        seqHits[y+1][1] =  seqHits[y+1][1] + 1    # Consistent design choice: the most N-proximal domain gets the extra AA position for no particular reason, we just need to handle this
                                        continue
                                elif seqHits[y+1][1] < seqHits[y][2]:
                                        # Handle identical E-values
                                        if seqHits[y][3] == seqHits[y+1][3]:
                                                len1 = seqHits[y][2] - seqHits[y][1]
                                                len2 = seqHits[y+1][2] - seqHits[y+1][1]
                                                # Handle identical lengths [just pick the first one, follows the consistent design decision we've made thus far and may result in less overlaps]
                                                if len1 == len2:
                                                        del seqHits[y+1]
                                                        break
                                                # Handle differing lengths
                                                else:
                                                        if len1 > len2:
                                                                del seqHits[y+1]
                                                                break
                                                        else:
                                                                del seqHits[y]
                                                                break
                                        # Get best E-value sequence
                                        else:
                                                bestEval = min(seqHits[y][3], seqHits[y+1][3])
                                                if bestEval == seqHits[y][3]:
                                                        # Figure out how much the worse E-value overlaps the best E-value
                                                        worseRange = range(seqHits[y+1][1], seqHits[y+1][2]+1)
                                                        betterRange = range(seqHits[y][1], seqHits[y][2]+1)
                                                        sharedPos = set(worseRange) & set(betterRange)
                                                        worseOvl = (len(worseRange) - len(set(worseRange) - sharedPos))
                                                        worsePerc = worseOvl / len(worseRange)
                                                        if worsePerc > ovlCutoff:
                                                                del seqHits[y+1]
                                                                break
                                                        else:
                                                                overlapping = 'n'       # We need to put this in here as an exit condition if we're looking at the last two entries and they don't overlap past our cutoff point. If we aren't looking at the last two entries, overlapping will == 'y' again when it goes back to the 'for y'... loop. If it IS the last two entries, when we continue we exit the 'for y' loop and immediately check if overlapping == 'n', which it does.
                                                                continue
                                                else:
                                                        # Figure out how much the worse E-value overlaps the best E-value
                                                        worseRange = range(seqHits[y][1], seqHits[y][2]+1)
                                                        betterRange = range(seqHits[y+1][1], seqHits[y+1][2]+1)
                                                        sharedPos = set(worseRange) & set(betterRange)
                                                        worseOvl = (len(worseRange) - len(set(worseRange) - sharedPos)) 
                                                        worsePerc = worseOvl / len(worseRange)
                                                        if worsePerc > ovlCutoff:
                                                                del seqHits[y]
                                                                break
                                                        else:
                                                                overlapping = 'n'
                                                                continue
                                else:
                                        overlapping = 'n'
                                        break
                return seqHits

        # Parse hmmer domblout file
        domDict = {}                            # We need to use a dictionary for later sorting since hmmsearch does not produce output that is ordered in the way we want to work with. hmmscan does, but it is SIGNIFICANTLY slower.
        with open(domtbloutFile, 'r') as fileIn:
                for line in fileIn:
                        if line.startswith('#'):
                                continue        # Skip header and footer of a domtblout file.
                        if line == '' or line == '\n':
                                continue        # Skip blank lines, shouldn't exist, but can't hurt
                        # Parse line and skip if evalue is not significant
                        sl = line.rstrip('\n').rstrip('\r').split()
                        e_value = float(sl[12])
                        if e_value > float(evalue):
                                continue
                        # Get relevant details
                        pid = sl[0]
                        #print(did)
                        #quit()
                        did = sl[3]
                        if not did.startswith('cath'):
                                did = os.path.basename(did)
                        dstart = int(sl[17])
                        dend = int(sl[18])
                        # Add into domain dictionary
                        if pid not in domDict:
                                domDict[pid] = [[did, dstart, dend, e_value]]
                        else:
                                domDict[pid].append([did, dstart, dend, e_value])

        # Delve into parsed hmmer dictionary and sort out overlapping domain hits from different databases
        finalDict = {}
        for key, value in domDict.items():
                value = list(value)
                # Collapse overlaps of identical domains
                uniqueModels = []
                for val in value:
                        uniqueModels.append(val[0])
                uniqueModels = list(set(uniqueModels))
                collapsedIdentical = []
                for model in uniqueModels:
                        modelGroup = []
                        for val in value:
                                if val[0] == model:
                                        modelGroup.append(val)
                        # Begin collapsing process
                        overlapping = 'y'
                        while True:
                                if len(modelGroup) == 1 or overlapping == 'n':
                                        break
                                for y in range(len(modelGroup)-1):
                                        if modelGroup[y+1][1] > modelGroup[y][2] and y != len(modelGroup)-2:
                                                continue
                                        elif modelGroup[y+1][1] == modelGroup[y][2]:
                                                modelGroup[y+1][1] =  modelGroup[y+1][1] + 1    # Consistent design choice: the most N-proximal domain gets the extra AA position for no particular reason, we just need to handle this
                                                continue
                                        elif modelGroup[y+1][1] < modelGroup[y][2]:
                                                # Calculate overlap proportion
                                                range1 = range(modelGroup[y][1], modelGroup[y][2]+1)            # +1 to offset Python counting up-to but not including the last value in a range
                                                range2 = range(modelGroup[y+1][1], modelGroup[y+1][2]+1)
                                                sharedPos = set(range1) & set(range2)
                                                r1Ovl = (len(range1) - len(set(range1) - sharedPos))
                                                r2Ovl = (len(range2) - len(set(range2) - sharedPos))
                                                r1Perc = r1Ovl / len(range1)
                                                r2Perc = r2Ovl / len(range2)
                                                highest = max(r1Perc, r2Perc)
                                                lowest = min(r1Perc, r2Perc)
                                                # Determine the length of the sequence extension of the most-overlapped sequence
                                                if highest == 0.50:
                                                        longest = max(r1Ovl, r2Ovl)
                                                        if longest == r1Ovl:
                                                                shortest = 'r2'
                                                                extension = len(range2) - r2Ovl
                                                        else:
                                                                shortest = 'r1'
                                                                extension = len(range1) - r1Ovl
                                                elif highest == r1Perc:
                                                        shortest = 'r1'
                                                        extension = len(range1) - r1Ovl
                                                else:
                                                        shortest = 'r2'
                                                        extension = len(range2) - r2Ovl
                                                ## Handle the various scenarios indicated by the highest/lowest values
                                                extensCutoff = 20       # This is arbitrary; can vary to test its effects
                                                # Block 1: small overlap of longer sequence
                                                if highest <= 0.20 and lowest <= 0.20:                                                  # small overlap of longer sequence / short overlap of shorter sequence
                                                        ## TRIM BASED ON E-VALUE
                                                        # Find best E-value
                                                        if modelGroup[y][3] < modelGroup[y+1][3]:
                                                                # Trim y+1
                                                                modelGroup[y+1] = [modelGroup[y+1][0], modelGroup[y][2]+1, *modelGroup[y+1][2:]]
                                                                continue
                                                        elif modelGroup[y+1][3] < modelGroup[y][3]:
                                                                # Trim y
                                                                modelGroup[y] = [*modelGroup[y][0:2], modelGroup[y+1][1]-1, modelGroup[y][3]]
                                                                continue
                                                        else:
                                                                # If the two E-value are identical, we just split down the middle!
                                                                splitPos = list(sharedPos)
                                                                splitPos.sort()
                                                                middle = findMiddle(splitPos)
                                                                if len(middle) == 1:
                                                                        modelGroup[y] = [*modelGroup[y][0:2], middle[0], modelGroup[y][3]]         # When we have an odd number, we just give the domain that is most N-proximal the extra AA position - it shouldn't realistically matter.
                                                                        modelGroup[y+1] = [modelGroup[y+1][0], middle[0]+1, *modelGroup[y+1][2:]]
                                                                        continue
                                                                else:
                                                                        modelGroup[y] = [*modelGroup[y][0:2], middle[1], modelGroup[y][3]]      # The findMiddle function returns reversed tuples like (181, 180) when the length of splitPos is even
                                                                        modelGroup[y+1] = [modelGroup[y+1][0], middle[0], *modelGroup[y+1][2:]]
                                                                        continue
                                                elif 0.20 < highest <= 0.80 and extension > extensCutoff and lowest <= 0.20:            # small overlap of longer sequence / intermediate overlap of shorter sequence [contingent on length]
                                                        # SPLIT MIDDLE
                                                        splitPos = list(sharedPos)
                                                        splitPos.sort()
                                                        middle = findMiddle(splitPos)
                                                        if len(middle) == 1:
                                                                modelGroup[y] = [*modelGroup[y][0:2], middle[0], modelGroup[y][3]]
                                                                modelGroup[y+1] = [modelGroup[y+1][0], middle[0]+1, *modelGroup[y+1][2:]]
                                                                continue
                                                        else:
                                                                modelGroup[y] = [*modelGroup[y][0:2], middle[1], modelGroup[y][3]]
                                                                modelGroup[y+1] = [modelGroup[y+1][0], middle[0], *modelGroup[y+1][2:]]
                                                                continue
                                                elif 0.20 < highest <= 0.80 and extension <= extensCutoff and lowest <= 0.20:           # small overlap of longer sequence / intermediate overlap of shorter sequence [contingent on length]
                                                        # JOIN
                                                        firstPos = modelGroup[y][1]
                                                        lastPos = max(modelGroup[y][2], modelGroup[y+1][2])
                                                        highestEval = max(modelGroup[y][3], modelGroup[y+1][3])
                                                        modelGroup[y] = [modelGroup[y][0], firstPos, lastPos, highestEval]
                                                        del modelGroup[y+1]
                                                        break
                                                elif highest > 0.80 and lowest <= 0.20:                                                 # small overlap of longer sequence / large overlap of shorter sequence
                                                        # JOIN
                                                        firstPos = modelGroup[y][1]
                                                        lastPos = max(modelGroup[y][2], modelGroup[y+1][2])
                                                        highestEval = max(modelGroup[y][3], modelGroup[y+1][3])
                                                        modelGroup[y] = [modelGroup[y][0], firstPos, lastPos, highestEval]
                                                        del modelGroup[y+1]
                                                        break
                                                # Block 2: intermediate overlap of longer sequence
                                                elif 0.20 <= highest <= 0.80 and extension > extensCutoff and 0.20 <= lowest <= 0.80:    # intermediate overlap of longer sequence / intermediate overlap of shorter sequence [contingent on length]
                                                        # SPLIT MIDDLE
                                                        splitPos = list(sharedPos)
                                                        splitPos.sort()
                                                        middle = findMiddle(splitPos)
                                                        if len(middle) == 1:
                                                                modelGroup[y] = [*modelGroup[y][0:2], middle[0], modelGroup[y][3]]
                                                                modelGroup[y+1] = [modelGroup[y+1][0], middle[0]+1, *modelGroup[y+1][2:]]
                                                                continue
                                                        else:
                                                                modelGroup[y] = [*modelGroup[y][0:2], middle[1], modelGroup[y][3]]
                                                                modelGroup[y+1] = [modelGroup[y+1][0], middle[0], *modelGroup[y+1][2:]]
                                                                continue
                                                elif 0.20 < highest <= 0.80 and extension <= extensCutoff and 0.20 <= lowest <= 0.80:   # intermediate overlap of longer sequence / intermediate overlap of shorter sequence [contingent on length]
                                                        # JOIN
                                                        firstPos = modelGroup[y][1]
                                                        lastPos = max(modelGroup[y][2], modelGroup[y+1][2])
                                                        highestEval = max(modelGroup[y][3], modelGroup[y+1][3])
                                                        modelGroup[y] = [modelGroup[y][0], firstPos, lastPos, highestEval]
                                                        del modelGroup[y+1]
                                                        break
                                                elif highest > 0.80 and 0.20 <= lowest <= 0.80:                                         # intermediate overlap of longer sequence / large overlap of shorter sequence
                                                        # JOIN
                                                        firstPos = modelGroup[y][1]
                                                        lastPos = max(modelGroup[y][2], modelGroup[y+1][2])
                                                        highestEval = max(modelGroup[y][3], modelGroup[y+1][3])
                                                        modelGroup[y] = [modelGroup[y][0], firstPos, lastPos, highestEval]
                                                        del modelGroup[y+1]
                                                        break
                                                # Block 3: large overlap of longer sequence
                                                elif highest > 0.80 and lowest > 0.80:                                                  # large overlap of longer sequence / large overlap of shorter sequence
                                                        # JOIN
                                                        firstPos = modelGroup[y][1]
                                                        lastPos = max(modelGroup[y][2], modelGroup[y+1][2])
                                                        highestEval = max(modelGroup[y][3], modelGroup[y+1][3])
                                                        modelGroup[y] = [modelGroup[y][0], firstPos, lastPos, highestEval]
                                                        del modelGroup[y+1]
                                                        break
                                                else:
                                                        print('This should never happen. I don\'t know how you did this, but I guess you have a right to complain to me about it... (just mention this error message)')
                                                        print(highest)
                                                        print(lowest)
                                                        print(modelGroup)
                                                        print(y)
                                                        quit()
                                        else:                                                           # We need the y != check above since we need to set an exit condition when no more overlaps are present. The if/elif will always trigger depending on whether there is/is not an overlap UNLESS it's the second last entry and there is no overlap. In this case we finally reach this else clause, and we trigger an exit.
                                                overlapping = 'n'
                                                break
                        # Add corrected individual models to collapsedIdentical list
                        collapsedIdentical += modelGroup
                # Process collapsedIdentical list to get our list of domains annotated against the sequence from each individual database
                ovlCutoff = 0.20        # This is also arbitrary; can test its effects
                if len(collapsedIdentical) == 1:
                        if key not in finalDict:
                                finalDict[key] = collapsedIdentical
                        else:
                                finalDict[key].append(collapsedIdentical)
                else:
                        collapsedIdentical.sort(key = lambda x: (x[1], x[2]))
                        overlapping = 'y'
                        collapsedIdentical = domain_ovl_resolver(ovlCutoff, collapsedIdentical)
                        if key not in finalDict:
                                finalDict[key] = collapsedIdentical
                        else:
                                finalDict[key].append(collapsedIdentical)

        # Generate output
        with open(outfile, 'w') as fileOut:
                for key, value in finalDict.items():
                        for val in value:
                                #fileOut.write(key + '\t' + '\t'.join(list(map(str, val))) + '\n')
                                fileOut.write(key + '\t' + str(val[1]) + '\t' + str(val[2]) + '\t' + str(val[0]) + '\n')

    def hmmer_grow(args, outdir, basename, group_dict, rejects_list, iterate):
        import os, re
        from Bio import SeqIO
        # Define functions
        def seqgrab(args, outdir, basename, seqid, start, stop):
            # Pull out the protein region
            infile = os.path.join(os.getcwd(), outdir, basename + '_cdhit.fasta')
            records = SeqIO.parse(open(infile, 'rU'), 'fasta')
            for record in records:
                if record.id == seqid:
                    aaseq = str(record.seq)
                    break
            aaseq = aaseq[int(start)-1:int(stop)]       # -1 to start to make this 0-indexed
            return aaseq
        def seqnamer(args, outdir, basename, seqid, start, stop):
            # Figure out which domain number this region should be
            infile = os.path.join(os.getcwd(), outdir, basename + '_unclustered_domains.fasta')
            records = SeqIO.parse(open(infile, 'rU'), 'fasta')
            seqnum = 1
            for record in records:
                baseid = '_'.join(record.id.split('_')[0:-3])
                if baseid == seqid:
                    seqnum += 1
            # Add to fasta
            newseqname = seqid + '_Domain_' + str(seqnum) + '_' + start + '-' + stop
            return newseqname
        
        # Let user know if this script is going to go through another iterative round
        print('Checking HMMER results for changes...')
        if iterate == 0:
            iterate = 1
        else:
            iterate = 2
        # Get hmmer output file details
        align_out_dir = os.path.join(os.getcwd(), outdir, 'tmp_alignments')
        hmmer_results = os.path.join(align_out_dir, 'tmp_hmmerParsed.results')
        # Parse the unclustered domains fasta file
        infile = os.path.join(os.getcwd(), outdir, basename + '_unclustered_domains.fasta')
        unclustDoms = list(SeqIO.parse(open(infile, 'rU'), 'fasta'))
        newDoms = []        # These will be added to the fasta file when we generate it at the end of the hmmer_grow function
        # Loop through the hmmer file and find its best match in the unclustered sequences file and make modifications indicated by hmmer
        with open(hmmer_results, 'r') as fileIn:        # Looping through the HMMER file line-by-line is good since we already have a dictionary structure if we want to look at other regions for the particular sequence ID on this line
            for line in fileIn:
                if line == '' or line == '\n':
                    continue        # Skip blank lines, shouldn't exist, but can't hurt
                # Parse line
                sl = line.rstrip('\n').rstrip('\r').split()
                pid = sl[0]
                dstart = sl[1]
                dend = sl[2]
                did = sl[3]
                hmmerrange = set(range(int(dstart), int(dend)+1))   # +1 to stop positions to keep everything 1-indexed
                # Find best match in the fasta file
                best = [0,0,0,0]
                secondBest = [0,0,0,0]
                positionsOverlapped = set()
                for x in range(len(unclustDoms)):
                    # Get sequence details
                    seqid_components = unclustDoms[x].id.split('_')
                    seqid = '_'.join(seqid_components[0:-3])
                    if seqid == pid:
                        start, stop = seqid_components[-1].split('-')
                        seqrange = set(range(int(start), int(stop)+1))
                        hmmerOvlAmount = 1 - (len(hmmerrange-seqrange) / len(hmmerrange))           # 1-(calc) means we're getting the amount of the HMMER hit that is overlapped (i.e., 1-[[100-40]/100] == 0.6, which means that 60% of it is overlapped)
                        seqOvlAmount = 1 - (len(seqrange-hmmerrange) / len(seqrange))               # Unsure if this is necessary yet
                        positionsOverlapped = positionsOverlapped.union(hmmerrange & seqrange)
                        #if hmmerOvlAmount > best[2] or seqOvlAmount > best[3]:
                        if hmmerOvlAmount > best[2]:
                            secondBest = best
                            best = [unclustDoms[x].id, hmmerrange, hmmerOvlAmount, seqOvlAmount]
                # Figure out if this hmmer region is novel
                if best[2] == 0:                                                                            # i.e., this sequence region does not overlap anything in the unclustered domains fasta
                    # Add unadulterated into the newDoms list for latter appending to fasta file
                    newseqname = seqnamer(args, outdir, basename, pid, dstart, dend)
                    newseq = seqgrab(args, outdir, basename, pid, dstart, dend)
                    newDoms.append([newseqname, newseq])
                    print('Added novel domain to file (' + newseqname + ')')
                    iterate = 0
                elif best[2] <= 0.2:                                                                        # i.e., this sequence only has slight overlap with other sequences. I use 0.2 since, theoretically, the most that we can trim off a hit here is 39% or so. That still means most of the sequence is left intact, so it may be a genuine domain region inbetween two other regions.
                    # Trim the mostly novel sequence and add into the newDoms list for latter appending to fasta file
                    newrange = hmmerrange-positionsOverlapped
                    newstart = str(min(newrange))
                    newend = str(max(newrange))
                    newseqname = seqnamer(args, outdir, basename, pid, newstart, newend)
                    newseq = seqgrab(args, outdir, basename, pid, newstart, newend)
                    print('Added trimmed novel domain to file (' + newseqname + ')')
                    iterate = 0
                # Figure out if there is a (nearly) guaranteed match in the fasta file. I'm allowing secondBest to be <= 0.1 since that means our main hit is still almost certainly the region that best matches the sequence that is incorporated in the HMM.
                elif secondBest[2] <= 0.1:
                    # Handle good matches
                    if best[2] >= 0.9 or best[3] >= 0.9:
                        print('Good match to HMMER hit: ' + pid + '_' + str(dstart) + '-' + str(dend) + ' : ' + best[0])
                        newseqname = pid + '_Domain_' + best[0].split('_')[-2] + '_' + str(dstart) + '-' + str(dend)
                        newseq = seqgrab(args, outdir, basename, pid, dstart, dend)
                        for x in range(len(unclustDoms)):
                            if unclustDoms[x].id == best[0]:
                                unclustDoms[x].id = newseqname
                                unclustDoms[x].seq = newseq
                    # Handle divergent matches
                    elif best[2] >= 0.5 or best[3] >= 0.5:     # note that for divergent matches we use 'and'. This is important to prevent fragmentary model matches which overlap an established sequence from overwriting the full length sequence which is part of the model.
                        print('Divergent match found: ' + pid + '_' + str(dstart) + '-' + str(dend) + ' : ' + best[0])
                        newseqname = pid + '_Domain_' + best[0].split('_')[-2] + '_' + str(dstart) + '-' + str(dend)
                        newseq = seqgrab(args, outdir, basename, pid, dstart, dend)
                        for x in range(len(unclustDoms)):
                            if unclustDoms[x].id == best[0]:
                                unclustDoms[x].id = newseqname
                                unclustDoms[x].seq = newseq
                    # Ignore poor matches
                    else:
                        print('Ignoring highly divergent overlap: ' + pid + '_' + str(dstart) + '-' + str(dend) + ' : ' + best[0])
                        print(best)
                        print(secondBest)
                # Handle ambiguity
                else:
                    # Handle the best case scenario for an ambiguous match. This still means we can be pretty sure that the best region is the correct match, but the amount it overlaps other regions is starting to be concerning.
                    if (best[2] >= 0.9 or best[3] >= 0.9) and (secondBest[2] <= 0.5 and secondBest[3] <= 0.5):
                        print('Ambiguous but good match to HMMER hit: ' + pid + '_' + str(dstart) + '-' + str(dend) + ' : ' + best[0])
                        newseqname = pid + '_Domain_' + best[0].split('_')[-2] + '_' + str(dstart) + '-' + str(dend)
                        newseq = seqgrab(args, outdir, basename, pid, dstart, dend)
                        for x in range(len(unclustDoms)):
                            if unclustDoms[x].id == best[0]:
                                unclustDoms[x].id = newseqname
                                unclustDoms[x].seq = newseq
                    # Stop handling sequences more ambiguous than this as a safety precaution
                    else:
                        print('Ambiguous overlap: ' + pid + '_' + str(dstart) + '-' + str(dend) + ' : ' + best[0])
                        print(best)
                        print(secondBest)
                        
        # Update fasta file
        if iterate != 2:        # If iterate == 2, then we're going to end the while loop. There may be some changes indicated by hmmer from this function, but we're assuming there will be very few since it's already been altered by hmmer once with no impact on the discovery of new domains.
            with open(os.path.join(os.getcwd(), outdir, basename + '_unclustered_domains.fasta'), 'w') as outfile:      # This overwrites the old file which is okay since we've already loaded the records in as a list
                for x in range(len(unclustDoms)):
                    seqid = unclustDoms[x].id
                    seq = str(unclustDoms[x].seq)
                    outfile.write('>' + seqid + '\n' + seq + '\n')     # The unclustered_domains.fasta should always end on a newline, so we don't need to preface our addition with '\n'
                for dom in newDoms:
                    outfile.write('>' + dom[0] + '\n' + dom[1] + '\n')
                
        # Return value to dictate whether we continue the looping
        return iterate

    def parse_joiner(args, outdir, basename, group_dict):
        import os, re
        from itertools import groupby
        from Bio import SeqIO
        from collections import Counter
        print('Performing domain boundary cleanup...')
        iterate2 = 'y'
        localiterate = 'y'
        ##### PARSING HMMERSEARCH
        align_out_dir = os.path.join(os.getcwd(), outdir, 'tmp_alignments')
        hmmer_file = open(os.path.join(align_out_dir, 'tmp_hmmer.results'))
        hmm_reg = re.compile(r'(Domain_\d{0,10})_align')
        grouper = lambda x: hmm_reg.search(x).group(0) if not x.startswith('#') else x[0]
        inside_grouper = lambda x: x[0]
        domdict = {}
        domlist = []
        for key, group in groupby(hmmer_file, grouper):
            if key.startswith('#'): continue
            group_hits = []
            for line in group:
                split_line = line.split()
                seqid = split_line[0]
                seqlen = split_line[2]
                dstart = int(split_line[17])
                dend = int(split_line[18])
                evalue = float(split_line[12])
                if evalue <= float(args['hmmevalnov']) and dend - dstart + 1 > int(args['cleanAA']):
                    group_hits.append([seqid, dstart, dend, seqlen])
            # Skip processing any groups that have no significant hmmer hits [should be impossible, but who knows?]
            if len(group_hits) == 0:
                continue
            # Loop through the accepted hits for this group and join any regions that require it
            updated_hits = []
            if len(group_hits) == 1:
                updated_hits = group_hits
            else:
                for seqid, regions in groupby(group_hits, inside_grouper):
                    regions = list(regions)
                    if len(regions) == 1:
                        updated_hits += regions
                        continue
                    regions.sort(key = lambda x: x[1])  # If we get to here, then we know we have multiple domain hits to the same sequence
                    while True:
                        overlap = 'n'
                        for i in range(0, len(regions)-1):
                            if regions[i+1][1] < regions[i][2] + int(args['cleanAA']):
                                regions[i+1] = [regions[i+1][0], regions[i][1], regions[i+1][2], regions[i+1][3]]
                                del regions[i]
                                overlap = 'y'
                                break
                            else:
                                #tmp.append(regions[i])
                                overlap = 'n'
                                continue
                        if overlap == 'n':
                            break
                    updated_hits += regions

            # Join overlapping domains
            for entry in updated_hits:
                seqid = entry[0]
                domlist.append([key] + entry)
                if seqid in domdict:
                    domdict[seqid].append([key] + entry[1:])
                else:
                    domdict[seqid] = [[key] + entry[1:]]
        # Sort domlist
        domlist.sort()
        #### FIND DOMAIN LENGTHS
        align_out_dir = os.path.join(os.getcwd(), outdir, 'tmp_alignments')
        msas = os.listdir(align_out_dir)
        dom_lengths = {}
        for msa in msas:
            if not msa.endswith('.fasta'):
                continue
            records = SeqIO.parse(open(os.path.join(align_out_dir, msa), 'rU'), 'fasta')
            for record in records:
                dom_lengths[msa.split('.')[0]] = len(str(record.seq))
        #### FIND CO-OCCURRING DOMAINS
        grouper = lambda x: x[0]
        domjoin = []
        vague_doms = [] # We'll make an output file from this listing domains that look like they might overlap, but don't overlap a majority of the time
        for key, group in groupby(domlist, grouper):    # key looks something like 'Domain_1_align'
            if localiterate == 'n': break
            domjoin = []
            nonejoin = []
            group = list(group)
            for val in group:       # val looks something like ['Domain_1_align', 'Seq1', 76, 154]
                seqdoms = domdict[val[1]]   # seqdoms looks something like [['Domain_1_align', 76, 154], ['Domain_2_align', 167, 203]] after sorting
                seqdoms.sort(key = lambda x: x[1])
                # SINGLE DOMAIN: Count the length of the leftover sequence and add it as a 'NONE' hit in domjoin
                if len(seqdoms) == 1:
                    leftover = int(seqdoms[0][3]) - int(seqdoms[0][2])
                    nonejoin.append(['None',leftover])
                    continue
                # MULTIDOMAIN [LAST DOMAIN]: Count the length of the leftover sequence and add it as a 'NONE' hit in domjoin [we don't add a 'continue' at the end of this loop, since it's possible this domain occurs twice in this sequence]
                if seqdoms[-1][0] == key:
                    leftover = int(seqdoms[-1][3]) - int(seqdoms[-1][2])
                    nonejoin.append(['None',leftover])
                # MULTIDOMAIN: Look to see if the current key/domain is followed by other domains
                for i in range(0, len(seqdoms)-1):
                    if seqdoms[i][0] == key:
                        if int(seqdoms[i+1][1]) < int(seqdoms[i][2]) + int(args['cleanAA']):
                            domjoin.append(seqdoms[i+1][0])
                        else:
                            nonejoin.append(['None',999])      # Make the leftover value very large to make it clear that nothing follows this domain (in this instance)
            # Tally votes, and see if any domain regions should be joined together
            if len(domjoin) == 0:
                # This domain doesn't overlap any other domains
                continue
            else:
                # Tally co-occurrences
                commonest_dom = Counter(domjoin).most_common(1)[0][0]
                commonest_count = Counter(domjoin).most_common(1)[0][1]
                percentage = commonest_count / sum(Counter(domjoin).values())
                if percentage < args['cooccur']:
                    # This domain doesn't overlap other domains frequently enough
                    continue
                # Update tally using 'None' counts
                commonest_len = dom_lengths[commonest_dom]
                nonetally = 0
                for none in nonejoin:
                    if int(none[1]) > commonest_len + args['cleanAA']:  # i.e., if the leftover sequence length is long enough to allow for a gap (cleanAA=30 by default) + the full domain length (commonest_len), then we assume it does not follow
                        nonetally += 1
                percentage = commonest_count / sum(Counter(domjoin).values()) + nonetally
                if percentage < args['cooccur']:
                    # This domain doesn't overlap other domains frequently enough
                    continue
                # If we get to here, then we assume that these two domains should overlap
                for val in group:       # We begin by re-looping through all the sequences that contain our current domain of interest (our 'key'), and joining the lengths of the 'key' + the overlapping domain ('commonest')
                    seqdoms = domdict[val[1]]
                    updated_seqdom = []
                    for i in range(len(seqdoms)):
                        if seqdoms[i][0] == commonest_dom:
                            continue
                        if seqdoms[i][0] != key:
                            updated_seqdom.append(seqdoms[i])
                        elif i != len(seqdoms)-1:   # i.e., if seqdoms[i][0] == key and this isn't the last entry, then we know that the next domain is our commonest_dom
                            joined_region = [key, seqdoms[i][1], seqdoms[i+1][2], seqdoms[i+1][3]]
                            updated_seqdom.append(joined_region)
                        else:                       # i.e., if we get to here, seqdoms[i] == key but it is also the last domain, so we know nothing follows after it
                            updated_seqdom.append(seqdoms[i])

                    domdict[val[1]] = updated_seqdom

                # Since we've joined domain groups, we now want to re-do the whole process of clustering, hmmer prediction, and joining
                print('Joined two domain regions together, reiterating the clustering process...')
                localiterate = 'n'
                break

        # Update the unclustered_domains.fasta file
        records = SeqIO.to_dict(SeqIO.parse(open(os.path.join(os.getcwd(), outdir, basename + '_cdhit.fasta'), 'rU'), 'fasta'))
        outlist = []
        for key2, value2 in domdict.items():
            seq = str(records[key2].seq)
            ongoingCount = 0
            for region in value2:
                ongoingCount += 1
                seqregion = seq[region[1]-1:region[2]]
                outlist.append('>' + key2 + '_Domain_' + str(ongoingCount) + '_' + str(region[1]) + '-' + str(region[2]) + '\n' + seqregion)
        outfile = open(os.path.join(os.getcwd(), outdir, basename + '_unclustered_domains.fasta'), 'w')
        outfile.write('\n'.join(outlist))
        outfile.close()
        
        # If we ever get to this point in the script and localiterate is still == 'y', then we know that there were no overlaps at all, and we have reached convergence
        if localiterate == 'y':
            print('No more domain regions overlap, beginning the final clustering process...')
            iterate2 = 'n'
        return iterate2
            
    def cluster_graph(clusterer, matrix):
        import seaborn as sns
        import matplotlib.pyplot as plt
        import sklearn
        asdf = clusterer.condensed_tree_.plot()
        color_palette = sns.color_palette('Paired', 20)
        cluster_colors = [color_palette[x] if x >= 0
                  else (0.5, 0.5, 0.5)
                  for x in clusterer.labels_]
        cluster_member_colors = [sns.desaturate(x, p) for x, p in
                                 zip(cluster_colors, clusterer.probabilities_)]
        eg = sklearn.manifold.TSNE().fit_transform(matrix.data)
        plt.scatter(*eg.T, s=50, linewidth=0, c=cluster_member_colors, alpha=0.25)
        plt.savefig('myfig')
