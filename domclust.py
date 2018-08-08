def alfree_matrix(fastaFile, reduceNum, alfAlgorithm):
        # Set up
        from alfpy import word_pattern, word_vector, word_distance
        from alfpy.utils import seqrecords, distmatrix
        from alfpy.utils.data import seqcontent
        # Read in unclustered domains file
        unclustDoms = open(fastaFile)
        records = seqrecords.read_fasta(unclustDoms)
        unclustDoms.close()
        # Extract details from records using alfpy-provided functions
        seqList = records.seq_list
        lengthList = records.length_list
        idList = records.id_list
        # Optional reduction of protein alphabet
        if str(reduceNum) == '15':
                murphy_15_tab = {"L":"L","V":"L","I":"L","M":"L","C":"C","A":"A","G":"G","S":"S","T":"T","P":"P","F":"F","Y":"F","W":"W","E":"E","D":"D","N":"N","Q":"Q","K":"K","R":"K","H":"H"}
                for i in range(len(seqList)):
                        newseq = ''
                        for letter in seqList[i]:
                                newseq += murphy_15_tab[letter]
                        seqList[i] = newseq
        elif str(reduceNum) == '11':
                eleven_tab = seqcontent.get_reduced_alphabet('protein')
                for i in range(len(seqList)):
                        newseq = ''
                        for letter in seqList[i]:
                                if letter in eleven_tab:
                                        newseq += eleven_tab[letter]
                                else:
                                        newseq += letter
                        seqList[i] = newseq
        elif reduceNum != 'n':
                print('I didn\'t recognise the reduceNum value provided to this function. It should be None, 11, or 15.')
                print('I\'m just going to treat this as None... if you don\'t want this behaviour, fix your input.')
        ### TEST
        # Compute distance matrix of choice
        if alfAlgorithm == 'google':
                # Calc for word_size = 2
                p = word_pattern.create(seqList, word_size=2)
                counts = word_vector.Counts(lengthList, p)
                dist = word_distance.Distance(counts, 'google')
                matrix2 = distmatrix.create(idList, dist)
                # Now 1
                p = word_pattern.create(seqList, word_size=1)
                counts = word_vector.Counts(lengthList, p)
                dist = word_distance.Distance(counts, 'google')
                matrix1 = distmatrix.create(idList, dist)
        elif alfAlgorithm == 'canberra':
                # 2
                p = word_pattern.create(seqList, word_size=2)
                weightmodel = word_vector.WeightModel(seqcontent.get_weights('protein'))
                counts = word_vector.CountsWeight(lengthList, p, weightmodel)
                dist = word_distance.Distance(counts, 'canberra')
                matrix2 = distmatrix.create(idList, dist)
                # 1
                p = word_pattern.create(seqList, word_size=1)
                weightmodel = word_vector.WeightModel(seqcontent.get_weights('protein'))
                counts = word_vector.CountsWeight(lengthList, p, weightmodel)
                dist = word_distance.Distance(counts, 'canberra')
                matrix1 = distmatrix.create(idList, dist)
        else:
                # 2
                p = word_pattern.create(seqList, word_size=2)
                counts = word_vector.Counts(lengthList, p)
                dist = word_distance.Distance(counts, 'braycurtis')
                matrix2 = distmatrix.create(idList, dist)
                # 1
                p = word_pattern.create(seqList, word_size=1)
                counts = word_vector.Counts(lengthList, p)
                dist = word_distance.Distance(counts, 'braycurtis')
                matrix1 = distmatrix.create(idList, dist)
        ### TESTING ###
        #import numpy as np
        #matrix2 = np.array(matrix2)
        #np.set_printoptions(threshold=np.inf)
        #matrix2.tofile('dataset.csv', sep=',', format="%s")
        #text = pyperclip.paste()
        #table = text.split(']\r\n [ ')
        #for i in range(len(table)):
        #        table[i] = table[i].replace('[', '')
        #        table[i] = table[i].replace(']', '')
        #for line in table:
        #        sl = line.split()
        #        out.append('\t'.join(sl))
        #pyperclip.copy('\r\n'.join(out))
        #print(idList)
        #quit()
        ### TESTING ###
        # Return value
        return matrix1, matrix2, idList

def cluster_hdb(leaf, singleclust, minsize, minsample, matrix1, matrix2, idList, rejects_list):
        # Set up
        import hdbscan
        group_dict = {}
        # Run clustering algorithm
        if leaf and singleclust:
                clusterer = hdbscan.HDBSCAN(metric='precomputed', cluster_selection_method = 'leaf', min_cluster_size = int(minsize), min_samples = int(minsample), allow_single_cluster = True)
        elif leaf and not singleclust:
                clusterer = hdbscan.HDBSCAN(metric='precomputed', cluster_selection_method = 'leaf', min_cluster_size = int(minsize), min_samples = int(minsample))
        elif not leaf and singleclust:
                clusterer = hdbscan.HDBSCAN(metric='precomputed', min_cluster_size = int(minsize), min_samples = int(minsample), allow_single_cluster = True)   
        #elif args['robust'] == 'y':
        #        #test robust single linkage
        #        clusterer = hdbscan.robust_single_linkage_.RobustSingleLinkage(metric='precomputed', gamma = min_cluster_size = int(args['minsize']))
        else:
                clusterer = hdbscan.HDBSCAN(metric='precomputed', min_cluster_size = int(minsize), min_samples = int(minsample))
        #clusterer.fit(matrix2.data)         # We look at matrix2 (or word_size == 2) first since it should, theoretically, find more 'similar' groups better than a word_size of 1
        clusterer.fit(matrix1.data) ## TESTING WORD SIZE 1 FIRST
        # Pull out domain groups
        clust_groups = clusterer.labels_
        clust_probs = clusterer.probabilities_
        print(clust_probs)
        # Sort groups
        for i in range(len(idList)):
                if clust_groups[i] != -1:
                        if clust_groups[i] not in group_dict:
                                group_dict[clust_groups[i]] = [idList[i]]
                        else:
                                group_dict[clust_groups[i]].append(idList[i])
                else:
                        rejects_list.append(idList[i])
        #### TESTING ####
        # Look through the results of the matrix1 (or word_size == 1) for any groups not found for word_size == 2
        #if args['leaf'] == 'y' and args['singleclust'] == 'y':
        #        clusterer = hdbscan.HDBSCAN(metric='precomputed', cluster_selection_method = 'leaf', min_cluster_size = int(args['minsize']), min_samples = int(args['minsample']), allow_single_cluster = True)
        #elif args['leaf'] == 'y' and args['singleclust'] == 'n':
        #        clusterer = hdbscan.HDBSCAN(metric='precomputed', cluster_selection_method = 'leaf', min_cluster_size = int(args['minsize']), min_samples = int(args['minsample']))
        #elif args['leaf'] == 'n' and args['singleclust'] == 'y':
        #        clusterer = hdbscan.HDBSCAN(metric='precomputed', min_cluster_size = int(args['minsize']), min_samples = int(args['minsample']), allow_single_cluster = True)   
        #else:
        #        clusterer = hdbscan.HDBSCAN(metric='precomputed', min_cluster_size = int(args['minsize']), min_samples = int(args['minsample']))
        #clusterer.fit(matrix2.data)         ## TESTING WORD SIZE 2 SECOND
        #new_groups = clusterer.labels_
        # Format new group into a dictionary
        #new_dict = {}
        #for i in range(len(idList)):
        #        if new_groups[i] != -1:
        #                if new_groups[i] not in new_dict:
        #                        new_dict[new_groups[i]] = [idList[i]]
        #                else:
        #                        new_dict[new_groups[i]].append(idList[i])
        # Look for groups not found in previous clustering
        #oldvals = []
        #ongoingCount = 0
        #oldlen = len(group_dict)
        #for val in group_dict.values():
        #        oldvals += val
        #skip = 'n'
        #for key, value in new_dict.items():
        #        for val in value:
        #                if val in oldvals:
        #                        skip = 'y'
        #                        break
        #        if skip == 'y':
        #                skip = 'n'
        #                continue
        #        else:
        #                group_dict[ongoingCount + oldlen] = value
        #                ongoingCount += 1
        #                for val in value:
        #                        del rejects_list[rejects_list.index(val)]           # Now that these are no longer rejects, we delete them from our rejects_list
        return group_dict, rejects_list

def tmpdir_setup(outputDir, tmpdirName):
        # Set up
        import os
        # Main function
        tmpDir = os.path.join(os.path.abspath(outputDir), tmpdirName)
        if os.path.isdir(tmpDir):
                for file in os.listdir(tmpDir):
                        filePath = os.path.join(tmpDir, file)
                        try:
                                if os.path.isfile(filePath):
                                        os.unlink(filePath)
                        except Exception as e:
                                print(e)
        else:
                os.mkdir(tmpDir)
        return tmpDir

def mafft_align(mafftdir, fastaFile, outputDir, threads, group_dict):
        # Set up
        import os, threading, math
        from Bio import SeqIO
        from Bio.Align.Applications import MafftCommandline
        # Define functions integral to this one
        def run_mafft(mafftdir, outputDir, fastaFile, startNum, endNum, thread):
                # Set up
                import platform
                for i in range(startNum, endNum):
                        # Extract sequences for this cluster and write to file
                        clustIDs = group_dict[i]           # Our group_dict value is indexed by sequential integers, so we can split work load easily by an iterating loop
                        tmpFasta = []
                        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
                        for record in records:
                                if record.id in clustIDs:
                                        tmpFasta.append('>' + record.id + '\n' + str(record.seq))
                        tmpName = os.path.join(outputDir, 'mafft_tmpfile_thread' + str(thread) + '.fasta')
                        with open(tmpName, 'w') as fileOut:
                                fileOut.write('\n'.join(tmpFasta))
                        # Run MAFFT
                        if platform.system() == 'Windows':
                                mafft_cline = MafftCommandline(os.path.join(mafftdir, 'mafft.bat'), input=tmpName)
                                stdout, stderr = mafft_cline()
                        else:
                                mafft_cline = MafftCommandline(os.path.join(mafftdir, 'mafft'), input=tmpName)
                                stdout, stderr = mafft_cline()
                        if stdout == '':
                                raise Exception('MAFFT error text below' + str(stderr))
                        # Process MAFFT output
                        stdout = stdout.split('\n')
                        while stdout[-1] == '\n' or stdout[-1] == '' or stdout[-1] == 'Terminate batch job (Y/N)?\n':   # Remove junk, sometimes MAFFT will have the 'Terminate ...' line
                                del stdout[-1]
                        stdout = '\n'.join(stdout)
                        # Create output alignment files
                        with open(os.path.join(outputDir, 'Domain_' + str(i) + '_align.fasta'), 'w') as fileOut:
                                fileOut.write(stdout)
                        # Clean up temp file
                        os.unlink(tmpName)

        # Set up threading requirements                         # This threading system is derived from chunk_fasta in (what is currently called) domfind.py
        dict_size = len(group_dict)
        rawNum = dict_size / threads                            # In cases where threads > dict_size, rawNum will be less than 1. numRoundedUp will equal the number of threads, and so we'll end up rounding these to 1. Yay!
        numRoundedUp = round((rawNum % 1) * threads, 0)         # By taking the decimal place and multiplying it by the num of threads, we can figure out how many threads need to be rounded up to process every cluster
        chunkPoints = []
        ongoingCount = 0
        for i in range(int(threads)):
                if i+1 <= numRoundedUp:                         # i.e., if two threads are being rounded up, we'll round up the first two loops of this
                        chunkPoints.append([ongoingCount, math.ceil(rawNum) + ongoingCount])    # Round up the rawNum, and also add our ongoingCount which corresponds to the number of clusters already put into a chunk
                        ongoingCount += math.ceil(rawNum)                                       # Unlike chunk_fasta, we're storing a paired value of ongoingCount and the chunk point
                else:                                                                           # Our mafft function iterates over a range, so we go up to and not including the last value; this system is compliant with that style of sorting
                        chunkPoints.append([ongoingCount, math.floor(rawNum) + ongoingCount])   # Also note that group_dict is indexed starting from 0, so if group_dict len == 10, we want to iterate over range(0,10) since the last actual index is 9
                        ongoingCount += math.floor(rawNum)
                if ongoingCount >= dict_size:                   # Without this check, if we have more threads than clusters, we can end up with "extra" numbers in the list (e.g., [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]).
                        break
        # Begin the loop
        processing_threads = []
        ongoingCount = 0
        for start, end in chunkPoints:
                build = threading.Thread(target=run_mafft, args=(mafftdir, outputDir, fastaFile, start, end, ongoingCount+1))
                processing_threads.append(build)
                build.start()
                print('Initiated thread num ' + str(ongoingCount+1) + ' for MAFFT alignment...')
                ongoingCount += 1

        # Wait for all threads to end.
        for process_thread in processing_threads:
                process_thread.join()
        print('MAFFT alignment completed')

def cluster_hmms(msaDir, hmmer3dir, concatName):
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
        # Concatenate HMMs
        concatHMM = os.path.join(msaDir, concatName)
        with open(concatHMM, 'w') as fileOut:
                for hmm in hmms:
                        fileOut.write(open(os.path.join(msaDir, hmm), 'r').read())
        # Press HMMs
        cmd = os.path.join(hmmer3dir, 'hmmpress') + ' -f "' + os.path.join(msaDir, concatHMM) + '"'
        run_hmmpress = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        hmmout, hmmerr = run_hmmpress.communicate()
        if hmmerr.decode("utf-8") != '':
                raise Exception('hmmpress error text below' + str(hmmerr.decode("utf-8")))

def hmmer_grow(domDict, unclustFasta, outputBase, group_dict, rejects_list, iterate):   # Outputbase refers to the path generated in the domain_finder script
        # Set up
        import os
        from Bio import SeqIO
        iterate += 1    # If we don't reset iterate to 0 by adding a new domain, this lets us ensure that the iterate value increases
        # Define functions integral to this one
        def seqgrab(outputBase, seqid, start, stop):
                # Pull out the protein region
                infile = os.path.join(outputBase + '_cdhit.fasta')   # We grab the hits from the cdhit.fasta file because the clean.fasta file might have a small stretch of low complexity region inside the sequence region that HMMER is hitting against
                records = SeqIO.parse(open(infile, 'r'), 'fasta')
                for record in records:
                        if record.id == seqid:
                                aaseq = str(record.seq)
                                aarecord = record
                                break
                aaseq = aaseq[int(start)-1:int(stop)]           # -1 to start to make this 0-indexed
                aarecord.seq = aaseq
                return aaseq, aarecord
        def seqnamer(unclustFasta, seqid, start, stop):
                # Figure out which domain number this region should be
                records = SeqIO.parse(open(unclustFasta, 'r'), 'fasta')
                seqnum = 1
                for record in records:
                        baseid = '_'.join(record.id.split('_')[0:-3])
                        if baseid == seqid:
                                seqnum = int(record.id.split('_')[-2])+1                # We do it this way because we modify the sequences in our record in-place, which means we might originally have Domain 1-3 for a baseid, but we might delete Domain 2. Thus, because the fasta files are ordered, the last match with the baseid should contain the highest domain number, so we just +1 to it.
                                #seqnum += 1
                # Add to fasta
                newseqname = seqid + '_Domain_' + str(seqnum) + '_' + str(start) + '-' + str(stop)
                return newseqname
        
        # Load in the unclustered domains fasta as a list of records [this lets us append items to it]
        unclustDoms = list(SeqIO.parse(open(unclustFasta, 'r'), 'fasta'))
        # Loop through the domDict object and find its best match in the unclustered sequences file and make modifications indicated by hmmer
        for key, value in domDict.items():
                for entry in value:
                        # Extract details from entry [this gives them informative value names, and also is a way of be re-using an older bit of code]
                        pid = key
                        dstart = entry[1]
                        dend = entry[2]
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
                                        hmmerOvlAmount = 1 - (len(hmmerrange-seqrange) / len(hmmerrange))                   # 1-(calc) means we're getting the amount of the HMMER hit that is overlapped (i.e., 1-[[100-40]/100] == 0.6, which means that 60% of it is overlapped)
                                        seqOvlAmount = 1 - (len(seqrange-hmmerrange) / len(seqrange))                           # Unsure if this is necessary yet
                                        positionsOverlapped = positionsOverlapped.union(hmmerrange & seqrange)
                                        #if hmmerOvlAmount > best[2] or seqOvlAmount > best[3]:
                                        if hmmerOvlAmount > best[2]:
                                                secondBest = best
                                                best = [unclustDoms[x].id, hmmerrange, hmmerOvlAmount, seqOvlAmount]
                        # Figure out if this hmmer region is novel
                        if best[2] == 0:                                                                                                                                                        # i.e., this sequence region does not overlap anything in the unclustered domains fasta
                                # Add unadulterated into the unclustDoms list
                                newseqname = seqnamer(unclustFasta, pid, dstart, dend)
                                newseq, newrec = seqgrab(outputBase, pid, dstart, dend)
                                newrec.id = newseqname
                                unclustDoms.append(newrec)
                                #newDoms.append([newseqname, newseq])                # Add this directly into the list
                                print('Added novel domain to file (' + newseqname + ')')
                                iterate = 0
                        elif best[2] <= 0.2:                                                                                                                                                # i.e., this sequence only has slight overlap with other sequences. I use 0.2 since, theoretically, the most that we can trim off a hit here is 39% or so because the secondBest hit will only overlap one side of the sequence. That still means most of the sequence is left intact, so it may be a genuine domain region inbetween two other regions.
                                # Trim the mostly novel sequence and add into the unclustDoms list
                                newrange = hmmerrange-positionsOverlapped
                                newstart = str(min(newrange))
                                newend = str(max(newrange))
                                newseqname = seqnamer(unclustFasta, pid, newstart, newend)
                                newseq, newrec = seqgrab(outputBase, pid, newstart, newend) # Add to list
                                newrec.id = newseqname
                                unclustDoms.append(newrec)
                                print('Added trimmed novel domain to file (' + newseqname + ')')
                                iterate = 0
                        # Figure out if there is a (nearly) guaranteed match in the fasta file. I'm allowing secondBest to be <= 0.1 since that means our main hit is still almost certainly the region that best matches the sequence that is incorporated in the HMM.
                        elif secondBest[2] <= 0.1:
                                # Handle good matches
                                if best[2] >= 0.9 or best[3] >= 0.9:
                                        print('Good match to HMMER hit: ' + pid + '_' + str(dstart) + '-' + str(dend) + ' : ' + best[0])
                                        newseqname = pid + '_Domain_' + best[0].split('_')[-2] + '_' + str(dstart) + '-' + str(dend)
                                        newseq, newrec = seqgrab(outputBase, pid, dstart, dend)                                 # We don't care about the newrec here since we just need to modify the existing record in place, we're not adding a new record or fusing existing ones
                                        for x in range(len(unclustDoms)):
                                                if unclustDoms[x].id == best[0]:
                                                        unclustDoms[x].id = newseqname
                                                        unclustDoms[x].seq = newseq
                                # Handle divergent matches
                                elif best[2] >= 0.5 or best[3] >= 0.5:         # note that for divergent matches we use 'and'. This is important to prevent fragmentary model matches which overlap an established sequence from overwriting the full length sequence which is part of the model.
                                        print('Divergent match found: ' + pid + '_' + str(dstart) + '-' + str(dend) + ' : ' + best[0])
                                        newseqname = pid + '_Domain_' + best[0].split('_')[-2] + '_' + str(dstart) + '-' + str(dend)
                                        newseq, newrec = seqgrab(outputBase, pid, dstart, dend)
                                        for x in range(len(unclustDoms)):
                                                if unclustDoms[x].id == best[0]:
                                                        unclustDoms[x].id = newseqname
                                                        unclustDoms[x].seq = newseq
                                # Handle poor matches
                                else:
                                        # Handle best case scenario
                                        if secondBest[2] == 0 and secondBest[3] == 0:
                                                print('Poor match found: ' + pid + '_' + str(dstart) + '-' + str(dend) + ' : ' + best[0])
                                                newseqname = pid + '_Domain_' + best[0].split('_')[-2] + '_' + str(dstart) + '-' + str(dend)
                                                newseq, newrec = seqgrab(outputBase, pid, dstart, dend)
                                                for x in range(len(unclustDoms)):
                                                        if unclustDoms[x].id == best[0]:
                                                                unclustDoms[x].id = newseqname
                                                                unclustDoms[x].seq = newseq
                                        # Ignore worse matches
                                        else:
                                                print('Ignoring highly divergent overlap: ' + pid + '_' + str(dstart) + '-' + str(dend) + ' : ' + best[0])
                                                print(best)
                                                print(secondBest)
                        # Handle ambiguity
                        else:
                                # Handle the best case scenario for an ambiguous match. This still means we can be pretty sure that the best region is the correct match.
                                if best[2] >= 0.9 and best[3] >= 0.9 and secondBest[2] <= 0.5:
                                        print('Ambiguous but good match to HMMER hit: ' + pid + '_' + str(dstart) + '-' + str(dend) + ' : ' + best[0])
                                        newseqname = pid + '_Domain_' + best[0].split('_')[-2] + '_' + str(dstart) + '-' + str(dend)
                                        newseq, newrec = seqgrab(outputBase, pid, dstart, dend)
                                        for x in range(len(unclustDoms)):
                                                if unclustDoms[x].id == best[0]:
                                                        unclustDoms[x].id = newseqname
                                                        unclustDoms[x].seq = newseq
                                # Matches that fall into this category appear to be candidates for fusing the two underlying unclustered domain sequences into a single one
                                else:
                                        print('Fusing overlap: ' + best[0] + ' : ' + secondBest[0])
                                        newseqname = pid + '_Domain_' + best[0].split('_')[-2] + '_' + str(dstart) + '-' + str(dend)
                                        newseq, newrec = seqgrab(outputBase, pid, dstart, dend)
                                        # Delete underlying sequences
                                        for i in range(0, 2):
                                                for x in range(len(unclustDoms)):
                                                        if unclustDoms[x].id == best[0] or unclustDoms[x].id == secondBest[0]:
                                                                del unclustDoms[x]
                                                                break
                                        # Append new sequence into file
                                        newrec.id = newseqname
                                        unclustDoms.append(newrec)
                                        print(newseq)
                                        #print(best)
                                        #print(secondBest)
                                        
        # Update fasta file
        if iterate != 2:                # If iterate == 2, then we're going to end the while loop. There may be some changes indicated by hmmer from this function, but we're assuming there will be very few since it's already been altered by hmmer once with no impact on the discovery of new domains.
                with open(outputBase + '_unclustered_domains.fasta', 'w') as outfile:          # This overwrites the old file which is okay since we've already loaded the records in as a list
                        for x in range(len(unclustDoms)):
                                seqid = unclustDoms[x].id
                                seq = str(unclustDoms[x].seq)
                                outfile.write('>' + seqid + '\n' + seq + '\n')         # The unclustered_domains.fasta should always end on a newline, so we don't need to preface our addition with '\n'
                        #for dom in newDoms:
                        #        outfile.write('>' + dom[0] + '\n' + dom[1] + '\n')
                        
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
        for key, group in groupby(domlist, grouper):        # key looks something like 'Domain_1_align'
                if localiterate == 'n': break
                domjoin = []
                nonejoin = []
                group = list(group)
                for val in group:           # val looks something like ['Domain_1_align', 'Seq1', 76, 154]
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
                                                nonejoin.append(['None',999])          # Make the leftover value very large to make it clear that nothing follows this domain (in this instance)
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
                        for val in group:           # We begin by re-looping through all the sequences that contain our current domain of interest (our 'key'), and joining the lengths of the 'key' + the overlapping domain ('commonest')
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
                                        else:                                           # i.e., if we get to here, seqdoms[i] == key but it is also the last domain, so we know nothing follows after it
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
