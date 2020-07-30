from alfpy import word_pattern, word_vector, word_distance
from alfpy.utils import seqrecords, distmatrix
from alfpy.utils.data import seqcontent
import os, threading, math, platform, hdbscan, copy
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet
from Bio.Align import MultipleSeqAlignment

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

def alfree_matrix(fastaFile, wordSize, reduceNum, alfAlgorithm):
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
        matrices = []           # Currently I'm not returning multiple matrices, but this exists to enable multiple word sizes to be used and combined
        wordSizes = [wordSize]  # As above, not used now, but maybe in the future this will be relevant
        for num in wordSizes:
                p = word_pattern.create(seqList, word_size=wordSize)
                if alfAlgorithm == 'canberra':
                        weightmodel = word_vector.WeightModel(seqcontent.get_weights('protein'))
                        counts = word_vector.CountsWeight(lengthList, p, weightmodel)
                else:
                        counts = word_vector.Counts(lengthList, p)
                dist = word_distance.Distance(counts, alfAlgorithm)
                matrices.append(distmatrix.create(idList, dist))
        # Return value
        return matrices, idList

## Consider https://hdbscan.readthedocs.io/en/latest/api.html for detecting outliers via all_points_mutual_reachability

def cluster_hdb(leaf, singleClust, minSize, minSample, matrixList, idList):
        # Set up
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
        for matrix in matrixList:                               # This lets us possibly feed in multiple matrices; currently I'm not doing that since it doesn't seem to really help all too much
                clusterer = hdbscan.HDBSCAN(metric='precomputed', cluster_selection_method = clustSelect, min_cluster_size = int(minSize), min_samples = int(minSample), allow_single_cluster = allowSingle)
                clusterer.fit(matrix.data)
                # Pull out domain groups
                clust_groups = clusterer.labels_
                # Sort groups
                groupDict = {}
                for i in range(len(idList)):
                        if clust_groups[i] != -1:
                                if clust_groups[i] not in groupDict:
                                        groupDict[clust_groups[i]] = [idList[i]]
                                else:
                                        groupDict[clust_groups[i]].append(idList[i])
                groupDicts.append(groupDict)
        # Merge word size matrices together if relevant
        if len(groupDicts) > 1:
                oldVals = []
                origLen = len(groupDicts[0])
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
                                groupDicts[0][ongoingCount + origLen] = value
                                ongoingCount += 1
        return groupDicts[0]

def tmpdir_setup(tmpDir):
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

def mafft_align_clust_dict(mafftdir, fastaFile, outputDir, prefix, suffix, threads, group_dict, algorithm):
        # Ensure that algorithm value is sensible
        if algorithm != None:   # If this is None, we'll just use default MAFFT
                if algorithm.lower() not in ['genafpair', 'localpair', 'globalpair']:
                        print('mafft_align: algorithm option must be an option in the below list. Fix this parameter and try again.')
                        print(['genafpair', 'localpair', 'globalpair'])
                        quit()
        # Define functions integral to this one
        def run_mafft(mafftdir, outputDir, prefix, suffix, fastaFile, startNum, endNum, thread, algorithm):
                # Set up
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
                        else:
                                mafft_cline = MafftCommandline(os.path.join(mafftdir, 'mafft'), input=tmpName)
                        if algorithm != None:
                                if algorithm.lower() == 'genafpair':
                                        mafft_cline.genafpair = True
                                elif algorithm.lower() == 'localpair':
                                        mafft_cline.localpair = True
                                elif algorithm.lower() == 'globalpair':
                                        mafft_cline.globalpair = True
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
                build = threading.Thread(target=run_mafft, args=(mafftdir, outputDir, prefix, suffix, fastaFile, start, end, ongoingCount+1, algorithm))
                processing_threads.append(build)
                build.start()
                ongoingCount += 1

        # Wait for all threads to end.
        for process_thread in processing_threads:
                process_thread.join()

def mafft_align_file_list(mafftdir, outputDir, fileList, threads, algorithm):
        # Ensure that algorithm value is sensible
        if algorithm != None:   # If this is None, we'll just use default MAFFT
                if algorithm.lower() not in ['genafpair', 'localpair', 'globalpair']:
                        print('mafft_align: algorithm option must be an option in the below list. Fix this parameter and try again.')
                        print(['genafpair', 'localpair', 'globalpair'])
                        quit()
        # Define functions integral to this one
        def run_mafft(mafftdir, outputDir, fileList, startNum, endNum, thread, algorithm):
                # Set up
                for i in range(startNum, endNum):
                        # Identify file to align
                        fastaFile = fileList[i]
                        # Run MAFFT
                        if platform.system() == 'Windows':
                                mafft_cline = MafftCommandline(os.path.join(mafftdir, 'mafft.bat'), input=fastaFile)
                        else:
                                mafft_cline = MafftCommandline(os.path.join(mafftdir, 'mafft'), input=fastaFile)
                        if algorithm != None:
                                if algorithm.lower() == 'genafpair':
                                        mafft_cline.genafpair = True
                                elif algorithm.lower() == 'localpair':
                                        mafft_cline.localpair = True
                                elif algorithm.lower() == 'globalpair':
                                        mafft_cline.globalpair = True
                        stdout, stderr = mafft_cline()
                        if stdout == '':
                                raise Exception('MAFFT error text below' + str(stderr))
                        # Process MAFFT output
                        stdout = stdout.split('\n')
                        while stdout[-1] == '\n' or stdout[-1] == '' or stdout[-1] == 'Terminate batch job (Y/N)?\n':   # Remove junk, sometimes MAFFT will have the 'Terminate ...' line
                                del stdout[-1]
                        stdout = '\n'.join(stdout)
                        # Create output alignment files
                        fileOutName = os.path.basename(fastaFile).rsplit(".", maxsplit=1)[0] + "_align.fasta"
                        with open(os.path.join(outputDir, fileOutName), 'w') as fileOut:
                                fileOut.write(stdout)
                        # Clean up temp file
                        os.unlink(fastaFile)

        # Set up threading requirements                         # This threading system is derived from chunk_fasta in (what is currently called) domfind.py
        list_size = len(fileList)
        rawNum = list_size / threads                            # In cases where threads > dict_size, rawNum will be less than 1. numRoundedUp will equal the number of threads, and so we'll end up rounding these to 1. Yay!
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
                if ongoingCount >= list_size:                   # Without this check, if we have more threads than clusters, we can end up with "extra" numbers in the list (e.g., [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]).
                        break
        # Begin the loop
        processing_threads = []
        ongoingCount = 0
        for start, end in chunkPoints:
                build = threading.Thread(target=run_mafft, args=(mafftdir, outputDir, fileList, start, end, ongoingCount+1, algorithm))
                processing_threads.append(build)
                build.start()
                ongoingCount += 1

        # Wait for all threads to end.
        for process_thread in processing_threads:
                process_thread.join()

def msa_trim(msaFastaIn, pctTrim, minLength, outType, msaFastaOut, indivSeqDrop, skipOrDrop, pctTrimType='presence'):
        from collections import Counter
        '''msaFastaIn is the path to the aligned MSA FASTA file to be trimmed
        pctTrim refers to the minimum proportion of sequences present in a single column to demarcate the start and end of an alignment
        minLength refers to the minimum length of a MSA after trimming before we decide to not trim at all; if this value is less than 1,
        we assume it's a ratio, otherwise it is an absolute length. outType influences whether this function returns a Biopython MSA object
        ("obj") or an output file ("file") msaFastaOut is only relevant when outType == file; otherwise it will be ignored
        '''
        # Ensure outType is sensible
        if outType.lower() not in ['obj', 'file', 'both']:
                print('msa_trim: This function requires an outType to be specified with a specific format.')
                print('Your provided value "' + outType + '" should instead be "obj", to return the Biopython MSA object, "file" to produce an output file which uses the string provided as msaFastaOut, or "both" to do both aforementioned things.')
                print('Format this correctly and try again.')
                quit()
        if (outType.lower() == 'file' or outType.lower() == 'both') and not type(msaFastaOut) == str:
                print('msa_trim: You specified a file output but didn\'t provide a string for msaFasta out - the file output name.')
                print('Format this correctly and try again.')
                quit()
        # Ensure pctTrimType is sensible
        if pctTrimType.lower() not in ['presence', 'identical']:
                print('msa_trim: This function requires a pctTrimType to be specified with a specific format.')
                print('Your provided value "' + pctTrimType + '" should does not meet specifications')
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
        # Process indivSeqDrop and ensure it is sensible
        if indivSeqDrop != None:
                try:
                        float(indivSeqDrop)
                        if not 0 <= float(indivSeqDrop) <= 1:
                                print('msa_trim: indivSeqDrop appears to be a float, but it is not a value from 0 to 1.')
                                print('Format this correctly and try again.')
                                quit()
                        indivSeqDrop = float(indivSeqDrop)
                except:
                        print('msa_trim: indivSeqDrop was not specified as None, but is also not capable of conversion to float.')
                        print('Format this correctly and try again.')
                        quit()
        # Process skipOrDrop and ensure it is sensible
        if skipOrDrop.lower() not in ['skip', 'drop']:
                print('msa_trim: skipOrDrop must equal "skip" or "drop"; I don\'t recognise ' + skipOrDrop + '.')
                print('Format this correctly and try again.')
                quit()
        # Load in fasta file as MSA object
        msa = AlignIO.read(msaFastaIn, 'fasta')
        # Loop through aligned columns and find the first position from the 5' end that meets our pctTrim value
        if pctTrimType.lower() == 'presence':
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
        elif pctTrimType.lower() == 'identical':
                for i in range(len(msa[0].seq)):
                        col = msa[:,i]
                        colCount = Counter(col)
                        identical = 0
                        for key, value in colCount.items():
                                if key != '-':
                                        if value > identical:
                                                identical = value
                        pctBases = identical / len(col)
                        if pctBases <= pctTrim:
                                continue
                        break
                # Same but for 3' end
                for x in range(len(msa[0].seq), 0, -1):
                        col = msa[:,x-1]
                        colCount = Counter(col)
                        identical = 0
                        for key, value in colCount.items():
                                if key != '-':
                                        if value > identical:
                                                identical = value
                        pctBases = identical / len(col)
                        if pctBases <= pctTrim:
                                continue
                        break

        # Check our values to ensure they're sensible
        if i >= x:      # If i >= x, that means we'd be trimming the sequence to 1bp or a negative value; in other words, we can't trim it at this pctTrim as printed below
                if skipOrDrop.lower() == 'skip':
                        print('#"' + os.path.basename(msaFastaIn) + '" can\'t be trimmed at this pctTrim value since no columns contain this proportion of sequences; no trimming will be performed.')
                        return msa              # If the user isn't expecting a returned object this should just disappear; if they want a file out, we won't modify it
                elif skipOrDrop.lower() == 'drop':
                        print('#"' + os.path.basename(msaFastaIn) + '" can\'t be trimmed at this pctTrim value since no columns contain this proportion of sequences; msa will be dropped.')
                        return None
        # Compare our MSA length post-trimming to our specified cut-offs to determine whether we're doing anything to this sequence or not
        seqLen = x - i          # This works out fine in 1-based notation
        if minLength < 1:
                ratio = seqLen / len(msa[0])
                if ratio < minLength:
                        if skipOrDrop.lower() == 'skip':
                                print('#"' + os.path.basename(msaFastaIn) + '" trimming reduces length more than minLength proportion cut-off; no trimming will be performed.')
                                return msa      # We're not going to make any changes if trimming shortens it too much
                        elif skipOrDrop.lower() == 'drop':
                                print('#"' + os.path.basename(msaFastaIn) + '" trimming reduces length more than minLength proportion cut-off; msa will be dropped.')
                                return None
        else:
                if seqLen < minLength:
                        if skipOrDrop.lower() == 'skip':
                                print('#"' + os.path.basename(msaFastaIn) + '" trimming reduces length more than absolute minLength cut-off; no trimming will be performed.')
                                return msa      # As above
                        elif skipOrDrop.lower() == 'drop':
                                print('#"' + os.path.basename(msaFastaIn) + '" trimming reduces length more than absolute minLength cut-off; msa will be dropped.')
                                return None
        # Trim our MSA object
        origMsa = copy.deepcopy(msa)    # Since we're going to be making changes to the msa from here on but still might want to return the unedited msa, we need to create a backup
        newMsa = MultipleSeqAlignment([])
        for y in range(len(msa)):       # Don't overwrite i from above! I made this mistake...
                msa[y].seq = Seq(str(msa[y].seq)[i:x], SingleLetterAlphabet())
                # Optionally remove sequences that don't appear to "fit" the alignment [i.e., contain a lot of gaps] according to user-specified gap proportion cut-off
                if indivSeqDrop != None:
                        gapCount = str(msa[y].seq).count('-')
                        gapProp = gapCount / len(msa[y].seq)
                        if gapProp > indivSeqDrop:
                                continue
                        newMsa.append(msa[y])
        # If we dropped sequences, make sure we don't have any blank columns now
        if indivSeqDrop != None and len(newMsa) > 1:
                for a in range(len(newMsa[0].seq), 0, -1):
                        col = newMsa[:,a-1]
                        if set(col) == {'-'}:
                                for b in range(len(newMsa)):
                                        newMsa[b].seq = Seq(str(newMsa[b].seq)[0:a-1] + str(newMsa[b].seq)[a:], SingleLetterAlphabet())
        # If we dropped sequences, ensure that our newMsa still has more than one entry in it
        if indivSeqDrop != None:
                if len(newMsa) < 2:
                        if skipOrDrop.lower() == 'skip':
                                print('#"' + os.path.basename(msaFastaIn) + '" removing gappy sequences according to indivSeqDrop cut-off means we do not have >= 2 sequences in this msa; no trimming will be performed.')
                                return origMsa
                        elif skipOrDrop.lower() == 'drop':
                                print('#"' + os.path.basename(msaFastaIn) + '" removing gappy sequences according to indivSeqDrop cut-off means we do not have >= 2 sequences in this msa; msa will be dropped.')
                                return None
                msa = newMsa
        # Return results either as the MSA object, as an output file, or as both
        if outType.lower() == 'file' or outType.lower() == 'both':
                with open(msaFastaOut, 'w') as fileOut:
                        fileOut.write(msa.format('fasta'))
        if outType.lower() == 'obj' or outType.lower() == 'both':
                return msa

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

def odseq_outlier_detect(msaFileNameList, rScriptDir, tmpDir, threshold, distMetric, bootStraps):
        # Set up
        import os, subprocess, pathlib
        ## Ensure input parameters are sensible
        # rScriptDir
        if rScriptDir == None:
                rScriptDir = ''         # We'll assume Rscript is locatable in PATH if unspecified
        elif rScriptDir != '' and not (os.path.isfile(os.path.join(rScriptDir, 'Rscript.exe')) or os.path.isfile(os.path.join(rScriptDir, 'Rscript'))):
                print('odseq_outlier_detect: rScriptDir does not appear to contain the Rscript file.')
                print('Fix your input and try again.')
                quit()
        # tmpDir
        if tmpDir == None:
                tmpDir = '.'            # We'll just put the file in the current directory if it isn't specified
        elif tmpDir != '' and not os.path.isdir(tmpDir):
                print('odseq_outlier_detect: tmpDir is not an existing directory or not able to be located.')
                print('Fix your input and try again.')
                quit()
        # threshold
        try:
                threshold = float(threshold)
        except:
                print('odseq_outlier_detect: threshold needs to be a float or capable of conversion to float.')
                print('Fix your input and try again.')
                quit()
        # distMetric
        if distMetric.lower() not in ['linear', 'affine']:
                print('odseq_outlier_detect: distMetric needs to be an option available in the list below. Fix your input and try again.')
                print(['linear', 'affine'])
                quit()
        # bootStraps
        try:
                bootStraps = int(bootStraps)
        except:
                print('odseq_outlier_detect: bootStraps needs to be an integer or capable of conversion to integer.')
                print('Fix your input and try again.')
                quit()
        # msaFileNameList
        if type(msaFileNameList) == str:
                msaFileNameList = [msaFileNameList]
        elif type(msaFileNameList) != list:
                print('odseq_outlier_detect: msaFileNameList type is not recognisable. It should be a list, but instead it is ' + str(type(msaFileNameList)) + '.')
                print('Fix your input and try again.')
                quit()
        if msaFileNameList == []:
                print('odseq_outlier_detect: msaFileNameList is empty. I don\'t know what to do in this situation since it shouldn\'t happen.')
                print('Code your call to this function properly to skip it.')
                quit()
        # Ensure that the msaFileNames are locatable
        ongoingCount = 0
        for fileName in msaFileNameList:
                if not os.path.isfile(fileName):
                        print('odseq_outlier_detect: index ' + str(ongoingCount) + ' in msaFileNameList is not able to be located. You might need to specify the full path to the file.')
                        print('Fix your input and try again.')
                        quit()
                msaFileNameList[ongoingCount] = os.path.abspath(fileName)
                ongoingCount += 1
        # Create script file
        scriptText = ['library("msa")', 'library("odseq")']
        for fileName in msaFileNameList:
                fileName = pathlib.Path(fileName).as_posix()
                scriptText.append('filename = "' + fileName + '"')
                scriptText.append('alig <- readAAMultipleAlignment(filename)')
                scriptText.append('y <- odseq(alig, threshold = {}, distance_metric = "{}", B = {})'.format(threshold, distMetric.lower(), bootStraps))
                scriptText.append('print(filename)')
                scriptText.append('print(y)')
        scriptFile = file_name_gen(os.path.join(tmpDir, 'tmp_odseq_script'), '.R')
        with open(scriptFile, 'w') as fileOut:
                fileOut.write('\n'.join(scriptText))
        # Format cmd
        cmd = '"' + os.path.join(rScriptDir, 'Rscript') + '" ' + scriptFile
        run_odseq = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        odseqout, odseqerr = run_odseq.communicate()
        odseqout = odseqout.decode("utf-8")
        if 'FALSE' not in odseqout and 'TRUE' not in odseqout:  # Rscript writes module loading details to stderr so we need to check that it worked properly in other ways
                raise Exception('Rscript ODseq error text below' + str(odseqerr.decode("utf-8")))
        # Parse ODseq results
        odseqTable = odseqout.split('[1]')
        if odseqTable[0] == '':
                del odseqTable[0]
        odseqDict = {}
        for i in range(len(odseqTable)):
                # Extract the chunk of information from this ODseq result
                chunk = odseqTable[i].split('\n')
                del chunk[0]    # We don't need the filename line anymore; printing it was useful enough to give us something to split by
                if chunk[-1] == '':
                        del chunk[-1]
                # Assemble the boolean results in a list [Note that this is an ordered list which corresponds to the input MSA, hence why there should be no need to parse file names for reassociation - we already have msaFileNameList]
                odseqResults = []
                for x in range(len(chunk)):
                        if x % 2 == 1:                  # i.e., if we're looking at an odd number in the list, it should be a boolean result line rather than sequence ID line
                                sl = chunk[x].split()   # It's probably overly cautious, but it does let us use sequences named 'TRUE' and 'FALSE' (why would you do this...) which a simple if == statement would not
                                odseqResults += sl
                # Convert 'TRUE' to True, 'FALSE' to False
                for x in range(len(odseqResults)):
                        if odseqResults[x] == 'TRUE':
                                odseqResults[x] = True
                        elif odseqResults[x] == 'FALSE':
                                odseqResults[x] = False
                        else:
                                print('odseq_outlier_detect: unrecognised output in odseqResults (' + str(odseqResults[x]) + '). What\'s going on?')
                                quit()
                # Add to our dictionary using index as key [Refer to comment in brackets above, this should be an easy data structure to work with since msaFileNameList[0]'s result will be odseqDict[0]]
                odseqDict[i] = odseqResults
        # Ensure everything worked fine
        if len(odseqDict) != len(msaFileNameList):
                print('odseq_outlier_detect: length of odseqDict != length of msaFileNameList. Inspect the below stderr report to see what went wrong and try to fix it.')
                print(odseqerr.decode("utf-8"))
                quit()
        # Clean up tmp file
        os.unlink(scriptFile)
        # Return results
        return odseqDict

def msa_outlier_detect(msaFileNameList, statsSave, removeIdentical):
        # Set up
        import hdbscan, statistics
        from Bio import AlignIO
        from Bio.Align import MultipleSeqAlignment
        import numpy as np
        # Set default HDBSCAN parameters for detecting outliers in MSA
        allowSingle = True      # Making this user tunable is probably incorrect.
        clustSelect = 'eom'     # This method is pretty imprecise as it stands,
        minSize = 2             # but it seems like these settings are the only
        minSample = 2           # ones that really "work"
        # Define functions integral to this one
        def sumpairs_score(pair):
                if pair[0].upper() == pair[1].upper():
                        return 1
                else:
                        return 0
        def pairwise_sumpairs_matrix(msa):
                spScore = []
                for a in range(len(msa)):
                        spScore.append([])
                        for b in range(len(msa)):
                                colScores = []
                                gapCount = 0
                                if a == b:
                                        spScore[-1].append(0.00)
                                        continue
                                for i in range(len(msa[a].seq)):
                                        pair = (msa[a][i], msa[b][i])
                                        if pair[0] != '-' and pair[1] != '-':
                                                colScores.append(sumpairs_score(pair))
                                        elif pair[0] == '-' and pair[1] == '-':         # Don't penalise an alignment that has a gap induced by another sequence
                                                continue
                                        elif set(msa[a][i:]) == {'-'}:                  # Don't penalise an alignment that has ended
                                                continue
                                        else:                                           # If this induces a gap in another alignment or is gapped itself, penalise it
                                                gapCount += 1
                                # Average column score with penalty for gaps
                                if colScores != []:
                                        colScore = 1 - (sum(colScores) / len(colScores)) * (1 - (gapCount / len(msa[a])))       # 1 - (..) since we want this to be a measure of dissimilarity in line with what alfpy does
                                        #colScore = 1 - (sum(colScores) / len(colScores))
                                else:
                                        colScore = 1
                                spScore[-1].append(colScore)
                return spScore
        
        # Ensure msaFileNameList is correctly formatted
        if type(msaFileNameList) == str:
                msaFileNameList = [msaFileNameList]
        elif type(msaFileNameList) != list:
                print('msa_outlier_detect: msaFileNameList type is not recognisable. It should be a list, but instead it is ' + str(type(msaFileNameList)) + '.')
                print('Fix your input and try again.')
                quit()
        if msaFileNameList == []:
                print('msa_outlier_detect: msaFileNameList is empty. I don\'t know what to do in this situation since it shouldn\'t happen.')
                print('Code your call to this function properly to skip it.')
                quit()
        # Ensure statsSave is correctly formatted
        if statsSave != True and statsSave != False:
                print('msa_outlier_detect: statsSave should equal True or False, not ' + str(statsSave) + '.')
                print('Fix your input and try again.')
                quit()
        # Ensure removeIdentical is correctly formatted
        if removeIdentical != True and removeIdentical != False:
                print('msa_outlier_detect: removeIdentical should equal True or False, not ' + str(removeIdentical) + '.')
                print('Fix your input and try again.')
                quit()
        # Loop through MSA files and try to identify outliers
        outlierDict = {}
        ongoingCount = 0
        for fileName in msaFileNameList:
                msa = AlignIO.read(fileName, 'fasta')
                # Remove identical sequences if relevant [these can skew our pairwise scoring matrix if highly similar sequence bits slip through]
                if removeIdentical == True:
                        madeChanges = False
                        newMsa = MultipleSeqAlignment([])
                        identicalPairs = []
                        for y in range(len(msa)):
                                for z in range(len(msa)):
                                        if y == z:
                                                continue
                                        elif msa[y].seq == msa[z].seq:
                                                identicalPairs.append([y, z])
                        if identicalPairs != []:
                                # Define our full identical group(s)
                                removeGroups = [identicalPairs[0]]      # Seed removeGroups with our first pair
                                for pair in identicalPairs:
                                        for n in range(len(identicalPairs)):
                                                if pair[0] in identicalPairs[n] or pair[1] in identicalPairs[n]:
                                                        if pair[0] not in identicalPairs[n]:
                                                                identicalPairs[n].append(pair[0])
                                                        elif pair[1] not in identicalPairs[n]:
                                                                identicalPairs[n].append(pair[1])
                                                else:
                                                        removeGroups.append(pair)
                                # Remove all but one entry from each identical group from our newMsa
                                for y in range(len(msa)):
                                        found = False
                                        for group in removeGroups:
                                                if y in group[1:]:      # This means we can allow the first number in each removeGroup
                                                        found = True
                                                        break
                                        if found == False:
                                                newMsa.append(msa[y])
                                if len(newMsa) > 2:     # We don't want to work with a structure with less than 3 sequences since our pairwise scoring becomes less impactful. It's a bit arbitrary in some respects,
                                        madeChanges = True
                                        origMsaLen = len(msa)
                                        msa = newMsa    # but this should help a domain model to not be overwhelmed by identical sequences while still letting us meaningfully use means and stdev for outlier detection.
                # Perform pairwise scoring
                spScore = pairwise_sumpairs_matrix(msa)
                spdm = np.array(spScore)
                # Cluster with HDBSCAN
                clusterer = hdbscan.HDBSCAN(metric='precomputed', cluster_selection_method = clustSelect, min_cluster_size = int(minSize), min_samples = int(minSample), allow_single_cluster = allowSingle)
                clusterer.fit(spdm)
                # Calculate basic statistics from pairwise scoring excluding HDBSCAN detected outliers
                if statsSave == True:
                        spMeanList = []
                        spAllMeansList = []
                        for i in range(len(spScore)):
                                spMean = statistics.mean(spScore[i][:i] + spScore[i][i+1:])     # Exclude self-match
                                spAllMeansList.append(spMean)
                                if clusterer.labels_[i] == -1:
                                        continue
                                spMeanList.append(spMean)
                        spMeansMean = statistics.mean(spMeanList)
                        spMeansPsdt = statistics.stdev(spMeanList)     # If len(spMeanList) == 1, stdev == 0. This makes it more likely we remove a legitimate sequence. TESTING: See if I should make this 10% of mean or something like that?
                # Convert HDBSCAN groups into boolean list of True == outlier, False == not outlier
                outlierList = []
                for i in range(len(clusterer.labels_)):
                        label = clusterer.labels_[i]
                        if label == -1 or (-1 not in clusterer.labels_ and 1 in clusterer.labels_):     # If the second condition is True, HDBSCAN didn't find any outliers but it did find at least 2 separate clusters. In this case, it seems
                                if statsSave == True:                                                   # appropriate to treat every sequence as a potential outlier, and use our statistical distribution to pick out huge outliers
                                        '''Note: This serves as a _really_ rough heuristic check to justify HDBSCAN's clustering decision. It involves 
                                        a basic comparison using pstdev to see if this row's mean SP score is anomalous compared to others. At the start
                                        is an additional hard cut-off check to make sure we don't remove something "a little bit" different to a group of 
                                        sequences that are otherwise very similar. I added this cut-off check as a result of manual inspection of results
                                        to prevent a mistake from happening to a specific cluster I was testing. The test scenario was this [0.26951219512195124,
                                        0.33048780487804874, 0.2, 0.22073170731707314, 0.2182926829268293, 0.21707317073170732]. 0.33 was detected as an outlier,
                                        but 0.33 is still really similar for a SP score. 0.2 * 2 == 0.4, and this helps to rescue our example. 0.5 seems to
                                        be a point where the cluster is no longer highly homogenous.
                                        The second hard cut-off check was derived from [0.5294117647058824, 0.47500000000000003, 0.46029411764705885, 0.5042016806722689,
                                        0.4613445378151261]. It exceeds our 0.5 cut-off so we want to be less lenient with it, but the first sequence is still
                                        quite similar to the others from manual inspection. 0.1 seems to be a good point where, even if the minimum mean is 0.5,
                                        a sequence with distance 0.6 to the others still looks "normal" in a MSA. 0.7 as the max for this cut-off is a bit
                                        arbitrary and it shouldn't really happen, but it's just to prevent any weirdness from happening (e.g., a cluster
                                        of 0.9 distances should not group with a 1.0 distance [even though a 0.9 cluster shouldn't exist])
                                        '''
                                        if spAllMeansList[i] < 0.5 and spAllMeansList[i] < min(spAllMeansList) * 2:
                                                outlierList.append(False)
                                        elif spAllMeansList[i] < 0.6 and spAllMeansList[i] < min(spAllMeansList) + 0.15:        # This is another hard cut-off case derived from real data
                                                outlierList.append(False)                                                       # tldr; 0.45 and 0.6 are compatible in a cluster
                                        elif spAllMeansList[i] < 0.7 and spAllMeansList[i] < min(spAllMeansList) + 0.1:
                                                outlierList.append(False)
                                        elif spAllMeansList[i] > spMeansMean + (1.5*spMeansPsdt):
                                                outlierList.append(True)
                                        else:
                                                outlierList.append(False)
                                else:
                                        outlierList.append(True)
                        else:
                                outlierList.append(False)
                # Add in previously deleted identical values if removeIdentical is specified and we made changes
                if removeIdentical == True:
                        if madeChanges == True:
                                for x in range(origMsaLen):
                                        found = False
                                        for group in removeGroups:
                                                if x in group[1:]:
                                                        found = group
                                        if found != False:
                                                outlierList.insert(x, outlierList[found[0]])    # This took me a bit of mental effort to devise, then a bit more to understand why it worked, but it's quite simple
                outlierDict[ongoingCount] = outlierList                                         # When we find an index that was removed, we just insert an identical copy of its remaining sequence result, the index of which is given by found[0]
                ongoingCount += 1
        return outlierDict

def outlier_dict_merge(dict1, dict2):
        # Ensure dict inputs are compatible
        if set(dict1.keys()) != set(dict2.keys()):
                print('outlier_dict_merge: dict1 and dict2 don\'t have identical keys. This shouldn\'t be true if they were produced using the same MSA file name list.')
                print('Fix your input and try again.')
                quit()
        # Merge into a new output dict
        mergedDict = {}
        for key in dict1.keys():
                value1 = dict1[key]
                value2 = dict2[key]
                if len(value1) != len(value2):
                        print('outlier_dict_merge: values within dict1 and dict2 don\'t have identical length. This shouldn\'t be true if they were produced using the same MSA file name list.')
                        print('Fix your input and try again. A bit of debug info is below.')
                        print('Key = ' + str(key) + '. Value1 = ' + str(value1) + '. Value2 = ' + str(value2) + '.')
                        quit()
                mergedList = []
                for i in range(len(value1)):
                        if value1[i] == True and value2[i] == True:     # If both outlier detection methods agree that it is an outlier, we mark it as outlier
                                mergedList.append(True)
                        else:                                           # If they disagree whether it's an outlier or both agree it is not an outlier, we marked it as not outlier
                                mergedList.append(False)
                mergedDict[key] = mergedList
        return mergedDict

def curate_msa_from_outlier_dict(outlierDict, msaFileNameList):
        # Set up
        from Bio import AlignIO
        from Bio.Align import MultipleSeqAlignment
        from Bio.Seq import Seq
        from Bio.Alphabet import SingleLetterAlphabet
        # Ensure msaFileNameList is correctly formatted
        if type(msaFileNameList) == str:
                msaFileNameList = [msaFileNameList]
        elif type(msaFileNameList) != list:
                print('curate_msa_from_outlier_dict: msaFileNameList type is not recognisable. It should be a list, but instead it is ' + str(type(msaFileNameList)) + '.')
                print('Fix your input and try again.')
                quit()
        if msaFileNameList == []:
                print('curate_msa_from_outlier_dict: msaFileNameList is empty. I don\'t know what to do in this situation since it shouldn\'t happen.')
                print('Code your call to this function properly to skip it.')
                quit()
        # Ensure outlierDict and msaFileNameList are compatible
        if len(outlierDict) != len(msaFileNameList):
                print('curate_msa_from_outlier_dict: msaFileNameList length is not identical to outlierDict length. This shouldn\'t be true if they were produced using the same MSA file name list.')
                print('Fix your input and try again. A bit of debug info is below.')
                print('Len outlierDict = ' + str(len(outlierDict)) + '. Len msaFileNameList = ' + str(len(msaFileNameList)) + '.')
                quit()
        # Loop through msaFileNameList and make modifications to MSAs in place
        for i in range(len(msaFileNameList)):
                # Parse MSA and make sure it corresponds to outlierDict correctly
                msa = AlignIO.read(msaFileNameList[i], 'fasta')
                outlierList = outlierDict[i]
                if len(msa) != len(outlierList):
                        print('curate_msa_from_outlier_dict: MSA file "' + msaFileNameList[i] + '" length is not the same as outlierDict entry. This shouldn\'t be true if they were produced using the same MSA file name list.')
                        print('Fix your input and try again.')
                        quit()
                newMsa = MultipleSeqAlignment([])
                for y in range(len(msa)):
                        if outlierList[y] == False:
                                newMsa.append(msa[y])
                # If we dropped sequences, ensure that our newMsa still has more than one entry in it [Note: this should technically never happen]
                if len(newMsa) < 2:
                        print('curate_msa_from_outlier_dict: MSA file "' + msaFileNameList[i] + '" after outlier removal has less than two entries. This shouldn\'t be possible, so something is wrong with the code.')
                        print('A bit of debug info is below.')
                        print('Len msa = ' + str(len(msa)) + '. Len newMsa = ' + str(len(newMsa)) + '. outlistList = ' + str(outlierList) + '.')
                        quit()
                # If we dropped sequences, make sure we don't have any blank columns now
                if len(newMsa) < len(msa):
                        for a in range(len(newMsa[0].seq), 0, -1):
                                col = newMsa[:,a-1]
                                if set(col) == {'-'}:
                                        for b in range(len(newMsa)):
                                                newMsa[b].seq = Seq(str(newMsa[b].seq)[0:a-1] + str(newMsa[b].seq)[a:], SingleLetterAlphabet())
                # Produce our updated MSA fasta file
                with open(msaFileNameList[i], 'w') as fileOut:
                        fileOut.write(newMsa.format('fasta'))

def file_name_gen(prefix, suffix):
        import os
        ongoingCount = 2
        while True:
                if not os.path.isfile(prefix + '1' + suffix):
                        return prefix + '1' + suffix
                elif os.path.isfile(prefix + str(ongoingCount) + suffix):
                        ongoingCount += 1
                else:
                        return prefix + str(ongoingCount) + suffix

def cluster_hmms(msaFileList, outputDir, hmmer3dir, concatName):
        # Set up
        import os, subprocess
        # Build HMMs from MSA directory
        hmms = []
        for msa in msaFileList:
                msa = os.path.abspath(msa)
                outputFileName = msa.rsplit('.', maxsplit=1)[0] + '.hmm'
                hmms.append(outputFileName)
                # Format cmd
                cmd = os.path.join(hmmer3dir, 'hmmbuild') + ' "' + os.path.join(outputDir, outputFileName) + '" "' + msa + '"'
                # Run hmmbuild
                run_hmmbuild = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                hmmout, hmmerr = run_hmmbuild.communicate()
                if hmmerr.decode("utf-8") != '':
                        raise Exception('hmmbuild error text below\n' + str(hmmerr.decode("utf-8")) + '\nMake sure that you define the -h3dir argument if this directory is not in your PATH')
        # Concatenate HMMs
        concatHMM = os.path.join(outputDir, concatName)
        with open(concatHMM, 'w') as fileOut:
                for hmm in hmms:
                        fileOut.write(open(hmm, 'r').read())
        # Press HMMs
        cmd = os.path.join(hmmer3dir, 'hmmpress') + ' -f "' + concatHMM + '"'
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

def coord_lists_overlap_cull(coordDict1, coordDict2, ovlPropCutoff):    # In this function, coordDict1 is the static value used for comparison. coordList2 is what we are culling overlaps from using ovlPropCutoff as our criteria.
        # Ensure ovlPropCutoff is sensible
        try:
                ovlPropCutoff = float(ovlPropCutoff)
        except:
                print('coord_lists_overlap_cull: ovlPropCutoff needs to be a float or capable of conversion to float')
                print('Fix your input and try again.')
                quit()
        if not 0 <= ovlPropCutoff <= 1:
                print('coord_lists_overlap_cull: ovlPropCutoff needs to be a float from 0 -> 1.')
                print('Fix your formatting and try again.')
                quit()
        # Loop through coordDict2 and find overlaps
        for key, value2 in coordDict2.items():
                if key not in coordDict1:
                        continue
                value1 = coordDict1[key] # I'm trying to keep consistent names; value1 comes from coordDict1
                # Compare all pairs from the lists of coords
                i = 0
                while True:
                        dropped = False
                        if i >= len(value2):
                                break
                        for x in range(len(value1)):
                                pair1 = value1[x]
                                pair2 = value2[i]
                                if pair1[1] >= pair2[0] and pair2[1] >= pair1[0]:
                                        minTail = min(pair1[1], pair2[1])
                                        maxStart = max(pair1[0], pair2[0])
                                        ovlLen = minTail - maxStart + 1     # When we know two ranges overlap, minTail (end portion of range) - maxStart gives our overlap length; +1 to offset 1-based stuff since 100-1 == length of 100
                                        ovlProp1 = ovlLen / (pair1[1] - pair1[0] + 1)
                                        ovlProp2 = ovlLen / (pair2[1] - pair2[0] + 1)
                                        if ovlProp1 >= ovlPropCutoff or ovlProp2 >= ovlPropCutoff:
                                                del value2[i]
                                                dropped = True
                                                break
                        if dropped == False:
                                i += 1
        # Delete empty dict keys now
        for key in list(coordDict2.keys()):
                if coordDict2[key] == []:
                        del coordDict2[key]
        return coordDict2
                

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

def parse_mms2tab_to_clusters(mms2Table):
        import math
        '''Note: Function currently assumes
        file is sorted by E-value.
        '''
        clustDict = {}
        newClustNumber = 0
        idPointerDict = {}
        evalueDict = {}
        EVALUE_REDUCTION = 2
        EVALUE_SAFETY_NET = -15 # REALLY arbitrary, might help clusters a bit?
        def update_cluster_evalues(clustDict, evalueDict, clustNumber, newEvalue):
                clustList = clustDict[clustNumber]
                for sid in clustList:
                        evalueDict[sid] = max([evalueDict[sid], newEvalue])
                
        # Load in file as a groupby iterator
        with open(mms2Table, 'r') as fileIn:
                for line in fileIn:
                        if line == "\n":
                                continue
                        sl = line.rstrip("\r\n").split("\t")
                        qid = sl[0]
                        tid = sl[1]
                        evalue = float(sl[10])
                        # Skip self hits
                        if qid == tid:
                                continue
                        # New qid hits
                        if qid not in idPointerDict:
                                # New qid + new tid hits
                                if tid not in idPointerDict:
                                        # Handle qid
                                        clustDict[newClustNumber] = [qid, tid]
                                        idPointerDict[qid] = newClustNumber
                                        evalueDict[qid] = evalue
                                        # Handle tid
                                        idPointerDict[tid] = newClustNumber
                                        evalueDict[tid] = evalue
                                        newClustNumber += 1
                                # New qid + old tid hits
                                else:
                                        passableEvalue = max([math.log10(evalueDict[tid]) / EVALUE_REDUCTION, EVALUE_SAFETY_NET])
                                        if math.log10(evalue) <= passableEvalue:
                                                # Obtain clustNumber
                                                clustNumber = idPointerDict[tid]
                                                # Handle qid
                                                evalueDict[qid] = evalue
                                                idPointerDict[qid] = clustNumber
                                                # Update tid cluster
                                                clustDict[clustNumber].append(qid)
                                                update_cluster_evalues(clustDict, evalueDict, clustNumber, evalue)
                        # Old qid hits
                        elif qid in idPointerDict:
                                # Old qid + new tid hits
                                if tid not in idPointerDict:
                                        passableEvalue = max([math.log10(evalueDict[qid]) / EVALUE_REDUCTION, EVALUE_SAFETY_NET])
                                        if math.log10(evalue) <= passableEvalue:
                                                # Obtain clustNumber
                                                clustNumber = idPointerDict[qid]
                                                # Handle tid
                                                evalueDict[tid] = evalue
                                                idPointerDict[tid] = clustNumber
                                                # Update qid cluster
                                                clustDict[clustNumber].append(tid)
                                                update_cluster_evalues(clustDict, evalueDict, clustNumber, evalue)
                                # Old qid + old tid hits
                                else:
                                        # Check if E-value is compatible
                                        passableEvalue1 = max([math.log10(evalueDict[qid]) / EVALUE_REDUCTION, EVALUE_SAFETY_NET])
                                        passableEvalue2 = max([math.log10(evalueDict[tid]) / EVALUE_REDUCTION, EVALUE_SAFETY_NET])
                                        if math.log10(evalue) <= passableEvalue1 and math.log10(evalue) <= passableEvalue2:
                                                # Obtain clustNumbers
                                                clustNumber1 = idPointerDict[qid]
                                                clustNumber2 = idPointerDict[tid]
                                                # Skip inverse cluster scenario
                                                if clustNumber1 == clustNumber2:
                                                        continue
                                                # Merge clusters
                                                clustDict[clustNumber1] += clustDict[clustNumber2]
                                                clustDict[clustNumber1] = list(set(clustDict[clustNumber1]))
                                                # Update pointers
                                                for sid in clustDict[clustNumber2]:
                                                        idPointerDict[sid] = clustNumber1
                                                del clustDict[clustNumber2]
                                                # Update cluster evalues
                                                update_cluster_evalues(clustDict, evalueDict, clustNumber1, evalue)
        # Update cluster numbers for uniformity
        finalClustDict = {}
        ongoingCount = 0
        for key, value in clustDict.items():
                finalClustDict[ongoingCount] = value
                ongoingCount += 1
        return finalClustDict
