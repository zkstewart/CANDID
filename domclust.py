def protein_alphabet_reduce(proteinList, reduceNum):
        # Set up
        elevenLetter = {'E':'D','L':'I','M':'I','Q':'K','R':'K','T':'S','V':'I','W':'F','Y':'F'}
        fifteenLetter = {"L":"L","V":"L","I":"L","M":"L","C":"C","A":"A","G":"G","S":"S","T":"T","P":"P","F":"F","Y":"F","W":"W","E":"E","D":"D","N":"N","Q":"Q","K":"K","R":"K","H":"H"}
        # Ensure reduceNum is sensible and hold onto the correct reduce dict OR return our unmodified protein list
        if reduceNum == None or reduceNum == 'n' or type(reduceNum) == bool:
                return proteinList      # No modification is necessary
        elif int(reduceNum) == 11:
                reduceDict = elevenLetter
        elif int(reduceNum) == 15:
                reduceDict = fifteenLetter
        else:
                print('I didn\'t recognise the reduceNum value provided to the protein_alphabet_reduce function. It should be None, \'n\', 11, or 15.')
                print('I\'m just going to treat this as None... if you don\'t want this behaviour, fix your input.')
                return proteinList      # No modification is necessary
        # Ensure our proteinList is a list; if a single str is provided, make it a list (then return the str back from the function later)
        listAtEnd = True
        if type(proteinList) == str:
                proteinList = [proteinList]
                listAtEnd = False
        # Main function
        for i in range(len(proteinList)):
                newseq = ''
                for letter in proteinList[i]:
                        if letter in reduceDict:
                                newseq += reduceDict[letter]
                        else:
                                newseq += letter
                proteinList[i] = newseq
        # Return our modified list
        if listAtEnd == False:
                proteinList = proteinList[0]
        return proteinList

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
        seqList = protein_alphabet_reduce(seqList, reduceNum)
        # Compute distance matrix for word sizes
        matrix = []
        for num in [2, 1]:      ## TESTING: Modify this to change the ordering of word size; the first value takes priority during HDBSCAN clustering
                p = word_pattern.create(seqList, word_size=num)
                if alfAlgorithm == 'canberra':
                        weightmodel = word_vector.WeightModel(seqcontent.get_weights('protein'))
                        counts = word_vector.CountsWeight(lengthList, p, weightmodel)
                else:
                        counts = word_vector.Counts(lengthList, p)
                dist = word_distance.Distance(counts, alfAlgorithm)
                matrix.append(distmatrix.create(idList, dist))
        # Return value
        return matrix[0], matrix[1], idList

def cluster_hdb(leaf, singleClust, minSize, minSample, matrix1, matrix2, idList):
        # Set up
        import hdbscan
        groupDicts = []
        # Convert our input values into relevant parameters
        if leaf == True:
                clustSelect = 'leaf'
        else:
                clustSelect = 'eom'     # This is the default HDBSCAN clustering method
        if singleClust == True:
                allowSingle = True
        else:
                allowSingle = False
        # Run clustering algorithm for both word size matrices
        for matrix in [matrix1, matrix2]:       # TESTING WORD SIZE 1 FIRST
                clusterer = hdbscan.HDBSCAN(metric='precomputed', cluster_selection_method = clustSelect, min_cluster_size = int(minSize), min_samples = int(minSample), allow_single_cluster = allowSingle)
                clusterer.fit(matrix.data) 
                # Pull out domain groups
                clust_groups = clusterer.labels_
                #clust_probs = clusterer.probabilities_
                #print(clust_probs)
                # Sort groups
                groupDict = {}
                for i in range(len(idList)):
                        if clust_groups[i] != -1:
                                if clust_groups[i] not in groupDict:
                                        groupDict[clust_groups[i]] = [idList[i]]
                                else:
                                        groupDict[clust_groups[i]].append(idList[i])
                groupDicts.append(groupDict)
        # Merge word size matrices together
        oldVals = []
        ongoingCount = 0
        for val in groupDicts[0].values():
                oldVals += val
        for key, value in groupDicts[1].items():
                skip = False
                for val in value:
                        if val in oldVals:
                                skip = True
                                break
                if skip == True:
                        continue
                else:
                        groupDicts[0][ongoingCount + len(groupDicts[0])] = value
                        ongoingCount += 1
        return groupDicts[0]

def tmpdir_setup(tmpDir):
        # Set up
        import os
        # Main function
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
                ongoingCount += 1

        # Wait for all threads to end.
        for process_thread in processing_threads:
                process_thread.join()

def msa_trim(msaFastaIn, pctTrim, minLength, outType, msaFastaOut):
        '''msaFastaIn is the path to the aligned MSA FASTA file to be trimmed
        pctTrim refers to the minimum proportion of sequences present in a single column to demarcate the start and end of an alignment
        minLength refers to the minimum length of a MSA after trimming before we decide to not trim at all; if this value is less than 1, we assume it's a ratio, otherwise it is an absolute length.
        outType influences whether this function returns a Biopython MSA object ("obj") or an output file ("file")
        msaFastaOut is only relevant when outType == file; otherwise it will be ignored
        '''
        # Set up
        from Bio import AlignIO
        from Bio.Seq import Seq
        from Bio.Alphabet import SingleLetterAlphabet
        # Ensure outType is sensible
        if outType.lower() not in ['obj', 'file']:
                print('msa_trim: This function requires an outType to be specified with a specific format.')
                print('Your provided value "' + outType + '" should instead be "obj", to return the Biopython MSA object, or "file" to produce an output file which uses the string provided as msaFastaOut.')
                print('Format this correctly and try again.')
                quit()
        if outType.lower() == 'file' and not type(msaFastaOut) == str:
                print('msa_trim: You specified a file output but did\'nt provide a string for msaFasta out - the file output name.')
                print('Format this correctly and try again.')
                quit()
        # Process minLength and ensure it is sensible
        try:
                int(minLength)
        except:
                print('msa_trim: minLength must be an integer or capable of conversion to integer.')
                print('Format this correctly and try again.')
                quit()
        if minLength < 0:
                print('msa_trim: minLength must be greater than 0.')
                print('Format this correctly and try again.')
                quit()
        # Load in fasta file as MSA object
        msa = AlignIO.read(msaFastaIn, 'fasta')
        # Loop through aligned columns and find the first position from the 5' end that meets our pctTrim value 
        for i in range(len(msa[0].seq)):
                col = msa[:,i]
                pctBases = 1 - (col.count('-') / len(col))
                if pctBases <= pctTrim:
                        continue
                break
        # Same but for 3' end
        for x in range(len(msa[0].seq), 0, -1):
                col = msa[:,x-1]
                pctBases = 1 - (col.count('-') / len(col))
                if pctBases <= pctTrim:
                        continue
                break
        # Check our values to ensure they're sensible
        if i >= x:
                print('MSA can\'t be trimmed at this pctTrim value; no columns contain this proportion of sequences!')
                return msa      # If the user isn't expecting a returned object this should just disappear; if they want a file out, we won't modify it
        # Compare our MSA length post-trimming to our specified cut-offs to determine whether we're doing anything to this sequence or not
        seqLen = x - i          # This works out fine in 1-based notation
        if minLength < 1:
                ratio = seqLen / len(msa[0])
                if ratio < minLength:
                        return msa      # We're not going to make any changes if trimming shortens it too much
        else:
                if seqLen < minLength:
                        return msa      # As above
        # Trim our MSA object
        for y in range(len(msa)):       # Don't overwrite i from above! I made this mistake...
                msa[y].seq = Seq(str(msa[y].seq)[i:x], SingleLetterAlphabet())
        # Return results either as the MSA object or as an output file
        if outType.lower() == 'obj':
                return msa
        else:
                with open(msaFastaOut, 'w') as fileOut:
                        for row in msa:
                                fileOut.write('>' + row.description + '\n' + str(row.seq) + '\n')

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

def coord_dict_compare(currentCoordDict, prevCoordDict):
        # Set up
        import copy
        currentCoordDict = copy.deepcopy(currentCoordDict)
        prevCoordDict = copy.deepcopy(prevCoordDict)
        # Comparison 1: same keys?
        if set(currentCoordDict.keys()) != set(prevCoordDict.keys()):
                return True     # This function asks the question "Did anything change?" - if this "if" statement is True, then the answer to this question is True
        # Comparison 2: same number of values associated with keys?
        for key in currentCoordDict.keys():
                currentVals = currentCoordDict[key]
                prevVals = prevCoordDict[key]
                if len(currentVals) != len(prevVals):
                        return True
                # Comparison 3: very similar coords associated with keys?
                mergedCoords = coord_lists_merge([currentVals, prevVals], 5)    # 5bp overlap means that minor alterations to coordinates won't matter, but substantial ones will
                if len(mergedCoords) != len(currentVals):
                        return True
        # If we get to here, nothing much changed
        return False

def coord_lists_merge(coordLists, offsetNum):   # offsetNum lets us specify increased stringency for overlaps; this is useful if we don't want 1 coord overlaps to cause merging
        # Pairwise coord list merging until a single list remains
        while len(coordLists) > 1:      # This is our exit condition; we merge coord lists until we have one remaining
                # Merge two lists together
                for x in range(len(coordLists[0])-1,-1,-1):
                        pair1 = coordLists[0][x]
                        merged = False
                        for y in range(len(coordLists[1])-1,-1,-1):
                                pair2 = coordLists[1][y]
                                # Detect overlap
                                if pair1[1] >= pair2[0] + offsetNum and pair2[1] >= pair1[0] + offsetNum:
                                        # Merge coords
                                        start = min([pair1[0], pair2[0]])
                                        end = max([pair1[1], pair2[1]])
                                        coordLists[1][y] = [start, end]
                                        merged = True
                        # If we couldn't merge the coord from [0] i.e., pair1, add it to [1]
                        if merged == False:
                                coordLists[1].append(pair1)
                # Delete the first list after it has been fully merged
                del coordLists[0]
        # Process the remaining list to merge self-overlaps
        coordList = coordLists[0]       # We need to extract the last list entry out of the containing list since it's not needed anymore
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

def run_hammock(hammockDir, javaDir, outputDir, threads, inputFasta):
        import os, subprocess, shutil
        # Remove output directory if it exists
        '''Hammock will instantly die if the folder exists. Since it is possible that a Hammock run will be interrupted,
        leaving the output folder intact, we need to remove it before we start the program run. Of course, this negates
        the whole reason for Hammock's behaviour (i.e., caution) but what can we do? Ideally, the user isn't doing things
        within the hammock_out folder...'''
        if os.path.isdir(outputDir):
                shutil.rmtree(outputDir)
        # Format command
        tmpDir = os.path.join(outputDir, 'tmp')
        cmd = os.path.join(javaDir, 'java') + ' -jar ' + os.path.join(hammockDir, 'Hammock.jar') + ' full -i {} -d {} -t {} --tmp {}'.format(*[inputFasta, outputDir, threads, tmpDir])
        # Run Hammock
        run_hammock = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        hammout, hammerr = run_hammock.communicate()
        # Ensure Hammock finished successfully [it writes everything to stderr, so it's difficult for us to, at a glance, figure out if any errors actually occurred; we need to parse stderr specifically
        checks = [False, False, False]
        for line in hammerr.decode("utf-8").split('\n'):
                if 'final_clusters' in line:
                        checks[0] = True
                elif 'Final system KLD over match' in line:
                        checks[1] = True
                elif 'Final system KLD over all MSA' in line:
                        checks[2] = True
        if not checks == [True, True, True]:
                print(checks)
                raise Exception('Hammock error text below' + str(hammerr.decode("utf-8")))

def parse_hammock(hammockOutFile, originalFasta):   # originalFasta refers to the fasta file used for running Hammock; we can pair sequence IDs back up by using Hammock's original_order.tsv file
        # Set up
        from Bio import SeqIO
        groupDict = {}
        unclustDict = {0: []}
        records = SeqIO.parse(open(originalFasta, 'r'), 'fasta')
        # Read through Hammock's output file
        with open(hammockOutFile, 'r') as fileIn:
                fileIn.readline()       # Skip the header line
                for line in fileIn:
                        sl = line.split('\t')
                        # Find out the sequence ID for this sequence
                        for record in records:
                                seqid = record.description
                                break
                        # Skip unclustered sequences
                        if sl[0] == 'NA':
                                unclustDict[0].append(seqid)
                                continue
                        # Add to our groupDict
                        if sl[0] in groupDict:
                                groupDict[sl[0]].append(seqid)
                        else:
                                groupDict[sl[0]] = [seqid]
        # Cull single-entry clusters and rename clusters
        ongoingCount = 0                        # This will be our cluster number iterator
        dictKeys = list(groupDict.keys())
        originallyEmpty = True
        for key in dictKeys:
                originallyEmpty = False         # This lets us assess whether our groupDict was always empty; if we enter this loop and then groupDict == {}, we know that we deleted a cluster
                if len(groupDict[key]) == 1:    # This is useful since, if Hammock gives all NAs for clusters, it probably means it's a single cluster; if it wasn't all NA's, then there probably isn't a cluster at all
                        del groupDict[key]
                        continue
                # Find out the cluster ID for this sequence
                clustNum = ongoingCount
                ongoingCount += 1
                # Rename the dict key
                groupDict[clustNum] = groupDict.pop(key)
        # Return our clusters if identified; otherwise, return the unclustDict object [not finding any clusters means there aren't any, or there is only 1]
        if groupDict != {}:
                return groupDict
        elif originallyEmpty == True:
                return unclustDict
        else:
                return None
