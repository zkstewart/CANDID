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

def mafft_align(mafftdir, fastaFile, outputDir, prefix, suffix, threads, group_dict):
        # Set up
        import os, threading, math
        from Bio import SeqIO
        from Bio.Align.Applications import MafftCommandline
        # Define functions integral to this one
        def run_mafft(mafftdir, outputDir, prefix, suffix, fastaFile, startNum, endNum, thread):
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
                        with open(os.path.join(outputDir, prefix + str(i) + suffix), 'w') as fileOut:
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
                build = threading.Thread(target=run_mafft, args=(mafftdir, outputDir, prefix, suffix, fastaFile, start, end, ongoingCount+1))
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
        import os
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
                print('"' + os.path.basename(msaFastaIn) + '" can\'t be trimmed at this pctTrim value; no columns contain this proportion of sequences!')
                return msa      # If the user isn't expecting a returned object this should just disappear; if they want a file out, we won't modify it
        # Compare our MSA length post-trimming to our specified cut-offs to determine whether we're doing anything to this sequence or not
        seqLen = x - i          # This works out fine in 1-based notation
        if minLength < 1:
                ratio = seqLen / len(msa[0])
                if ratio < minLength:
                        print('"' + os.path.basename(msaFastaIn) + '" trimming reduces length more than minLength proportion cut-off; no trimming will be performed.')
                        return msa      # We're not going to make any changes if trimming shortens it too much
        else:
                if seqLen < minLength:
                        print('"' + os.path.basename(msaFastaIn) + '" trimming reduces length more than absolute minLength cut-off; no trimming will be performed.')
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

def msa_score(msa, scoringMethod, cutoff):
        # Set up
        import os
        from Bio import AlignIO
        import Bio.Align
        # Define functions integral to this one
        def blosum_score(pair):
                # Set up
                from Bio.SubsMat import MatrixInfo
                blosum = MatrixInfo.blosum62    # blosum62 seems to be the most commonly used matrix, I'll lean on consensus w/r/t which one is probably best here
                # Make sure pair is formatted correctly
                pair = (pair[0].upper(), pair[1].upper())
                if pair not in blosum:
                        pair = tuple(reversed(pair))
                # Return the score
                return blosum[pair]
        def sumpairs_score(pair):
                if pair[0].upper() == pair[1].upper():
                        return 1
                else:
                        return 0
        # Determine what input we are receiving
        if type(msa) == Bio.Align.MultipleSeqAlignment:
                fileType = 'obj'
        elif type(msa) == str:
                if not os.path.isfile(msa):
                        print('msa_score: You\'ve provided a string input but it isn\'t detected as a file.')
                        print('Did you type the file name correctly, or do you need to provide the full path to it? Fix your input and try again.')
                        quit()
                fileType = 'file'
        else:
                print('msa_score: You\'ve provided an input type I am not compatible with.')
                print('The input should be a Bio.Align.MultipleSeqAlignment object or a string indicating the file location. Fix your input and try again.')
                quit()
        # If we're working with a file input, load it as an obj
        if fileType == 'file':
                msa = AlignIO.read(msa, 'fasta')
        # Loop through msa columns and perform pairwise scoring
        overallScores = []
        for i in range(len(msa[0].seq)):
                col = msa[:,i]
                colScores = []
                # Column scores based on possible pairs exclusive of '-'
                for x in range(len(col)-1):
                        for y in range(x+1, len(col)):
                                pair = (col[x].upper(), col[y].upper())
                                if pair[0] != '-' and pair[1] != '-':
                                        if scoringMethod == 'blosum':
                                                colScores.append(blosum_score(pair))
                                        elif scoringMethod == 'sp':
                                                colScores.append(sumpairs_score(pair))
                                else:
                                        continue
                # Average column score with penalty for gaps
                if colScores != []:
                        colScore = (sum(colScores) / len(colScores)) * (1 - (col.count('-') / len(col)))
                else:
                        colScore = 0    # If there is only a single letter, it receives no score
                overallScores.append(colScore)
        # Compare our overall score against cutoff and determine if this alignment is "good enough"
        finalScore = sum(overallScores) / len(overallScores)
        if finalScore >= cutoff:
                return True
        else:
                return False

def testing_stuff():
        import os
        from Bio import AlignIO
        scoringMethod = 'sp'
        # Define functions integral to this one
        def blosum_score(pair):
                # Set up
                from Bio.SubsMat import MatrixInfo
                blosum = MatrixInfo.blosum62    # blosum62 seems to be the most commonly used matrix, I'll lean on consensus w/r/t which one is probably best here
                # Make sure pair is formatted correctly
                pair = (pair[0].upper(), pair[1].upper())
                if pair not in blosum:
                        pair = tuple(reversed(pair))
                # Return the score
                return blosum[pair]
        def sumpairs_score(pair):
                if pair[0].upper() == pair[1].upper():
                        return 1
                else:
                        return 0
        fastaFiles = []
        for file in os.listdir('/home/lythl/Desktop/CANDID/newtest/tmp_alignments'):
                if file.endswith('.fasta'):
                        fastaFiles.append(file)
        for fasta in fastaFiles:
                scoringMethod = 'sp'
                msa = AlignIO.read(os.path.join('/home/lythl/Desktop/CANDID/newtest/tmp_alignments', fasta), 'fasta')
                overallScores = []
                for i in range(len(msa[0].seq)):
                        col = msa[:,i]
                        colScores = []
                        # Column scores based on possible pairs exclusive of '-'
                        for x in range(len(col)-1):
                                for y in range(x+1, len(col)):
                                        pair = (col[x].upper(), col[y].upper())
                                        if pair[0] != '-' and pair[1] != '-':
                                                if scoringMethod == 'blosum':
                                                        colScores.append(blosum_score(pair))
                                                elif scoringMethod == 'sp':
                                                        colScores.append(sumpairs_score(pair))
                                        else:
                                                continue
                        # Average column score with penalty for gaps
                        if colScores != []:
                                colScore = (sum(colScores) / len(colScores)) * (1 - (col.count('-') / len(col)))
                        else:
                                colScore = 0    # If there is only a single letter, it receives no score
                        overallScores.append(colScore)
                print(sum(overallScores) / len(overallScores))

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
