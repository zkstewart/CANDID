class benchparse:
    def benchparse(args, outdir, basename):
        import argparse, os, re
        from itertools import groupby
        from collections import Counter

        def hmmersort(args, outdir, basename):
            # Script modified from one written by Andrzej Zielezinski (http://www.staff.amu.edu.pl/~andrzejz/)
            import os
            print('Parsing hmmsearch output...')
            # Load in file
            fh = open(os.path.join(os.getcwd(), outdir, basename + '_hmmer.results'), 'r')
            #fh = open(args.file, 'r')
            domdict = {}
            # Loop through file
            for line in fh:
                if line.startswith('#'):
                    continue  # Skip header and footer of a domtblout file.
                if line == '' or line == '\n':
                    continue     # Skip blank lines, shouldn't exist, but can't hurt
                sl = line.split()
                evalue = float(sl[12])
                if evalue > float(args['hmmeval']):
                    continue
                pid = sl[0]
                if sl[3].startswith('cath'):
                    did = sl[3]
                else:
                    did = os.path.basename(sl[3])
                # Start of removal argument handling
                #if did in args['removal']:
                #    print('Skipping ' + did + ': ' + pid)
                #    continue
                # End of removal argument handling
                dstart = int(sl[17])
                dend = int(sl[18])
                # Add into domain dictionary    # We need to use a dictionary for later sorting since hmmsearch does not produce output that is ordered in the way we want to work with. Hmmscan does, but it is SIGNIFICANTLY slower.
                if pid not in domdict:          # This unfortunately serves as the only aspect of this script that should pose potential problems for memory use. To my knowledge there is no way to re-order hmmsearch output into the same format as hmmscan that is both fast and memory efficient.
                    domdict[pid] = [[did, dstart, dend, evalue]]
                else:
                    domdict[pid].append([did, dstart, dend, evalue])
            fh.close()
            # Sort
            for key, value in domdict.items():
                value.sort(key = lambda x: x[1])
                domdict[key] = value
            return domdict

        # Sort hmmersearch output
        domdict = hmmersort(args, outdir, basename)

        # Loop through domdict
        for key, value in domdict.items():
            # If more than one domain prediction is present, look for overlaps
            if len(value) == 1:
                continue
            oldval = value[:]
            done = 'n'
            while True:
                if done == 'y' or len(value) == 1:
                    break
                for i in range(0, len(value)-1):
                    if value[i][2] > value[i+1][1]: # If the domain predictions overlap...
                        if value[i][3] < value[i+1][3]: # If the first domain prediction has a better E-value...
                            del value[i+1]
                            break
                        elif value[i][3] == value[i+1][3]:  # If the two domain predictions have equal E-values...
                            val1range = value[i][2] - value[i][1]
                            val2range = value[i+1][2] - value[i+1][1]
                            if val1range == val2range:          # If the two domain predictions also have equal length, we delete the one that (presumably) has a later starting point...
                                del value[i+1]
                                break
                            else:
                                longer_range = max(val1range, val2range)
                                if longer_range == val1range:       # If the two domain predictions differ in length, we choose the longer domain prediction...
                                    del value[i+1]
                                    break
                                else:
                                    del value[i]
                                    break
                        else:                                   # If the second domain prediction has a better E-value...
                            del value[i]
                            break
                    elif i == len(value)-2:
                        done = 'y'
                        break
                    else:
                        continue
            # Figure out what sections have been cut out, and make sure that we re-mask these regions
            domregions = set()
            for val in value:
                domregions = domregions.union(set(range(val[1], val[2]+1)))
            maskregions = set()
            for val in oldval:
                if val not in value:
                    maskregions = maskregions.union(set(range(val[1], val[2]+1)))
            interregions = domregions & maskregions
            masked = list(sorted(list(maskregions - interregions)))
            if masked == []:
                donothing = 0       # masked == [] means the best E-value region is also the biggest region in general, and overlaps all other regions
            else:
                # Reorganise the sections to be masked
                mask_list = []
                firstnum = masked.pop(0)
                lastnum = firstnum
                for num in masked:
                    if num == lastnum + 1:
                        lastnum = num
                        continue
                    else:
                        mask_list.append([firstnum, lastnum])
                        lastnum = num
                        firstnum = num
                mask_list.append([firstnum, num])
                # Add masking regions back into the dictionary for file output
                for region in mask_list:
                    value.append(['mask', region[0], region[1], 0])
            value.sort(key=lambda x:x[1])
            # Remove any values that should be allowed to exit for benchmarking
            for i in range(len(value)-1, -1, -1):
                if value[i][0] in args['removal']:
                    print('Skipping ' + value[i][0] + ': ' + key)
                    del value[i]
            domdict[key] = value

        # Produce output
        ongoingCount = 0
        outname = os.path.join(os.getcwd(), outdir, basename + '_hmmerParsed.results')
        outputText = []
        dom_list = []
        for key, value in domdict.items():
            value.sort(key = lambda x: x[1])
            #print(value)
            #quit()
            for domregion in value:
                outputText.append(key + '\t' + str(domregion[1]) + '\t' + str(domregion[2]) + '\t' + str(domregion[0]))
                dom_list.append(domregion[0])
            ongoingCount += 1
            # Save to reduce memory consumption
            if ongoingCount%10000 == 0 and os.path.isfile(outname) == False:
                output = open(outname, 'w')
                output.write('\n'.join(outputText))
                output.close()
                outputText = []
            elif ongoingCount%10000 == 0 and os.path.isfile(outname) == True:
                output = open(outname, 'a')
                output.write('\n'.join(outputText))
                output.close()
                outputText = []
        # Dump the last few results after this process has finished
        if os.path.isfile(outname) == False:
                output = open(outname, 'w')
                output.write('\n'.join(outputText))
                output.close()
        elif os.path.isfile(outname) == True:
                output = open(outname, 'a')
                output.write('\n'.join(outputText))
                output.close()
        # Produce counter output so user can determine what to remove for benching
        counted = Counter(dom_list)
        outname = os.path.join(os.getcwd(), outdir, basename + '_domCounts.txt')
        outputText = ['key\toccurrence']
        for pair in counted.most_common():
            outputText.append(str(pair[0]) + '\t' + str(pair[1]))
        output = open(outname, 'w')
        output.write('\n'.join(outputText))
        output.close()

    def reject_novelty(args, outdir, basename):
        import os
        from Bio import SeqIO
        # Get reject regions to mask
        records = SeqIO.parse(open(os.path.join(os.getcwd(), outdir, args['rejects']), 'rU'), 'fasta')
        rejects = {}
        for record in records:
            seqname = '_'.join(record.id.split('_')[0:-3])
            reject_range = record.id.split('_')[-1].split('-')
            reject_range = [int(reject_range[0]), int(reject_range[1])]
            if seqname not in rejects:
                rejects[seqname] = [reject_range]
            else:
                rejects[seqname].append(reject_range)
        # Mask reject regions from segcoils file
        coils_name = os.path.join(os.getcwd(), outdir, basename + '_segcoils.fasta')
        records = SeqIO.to_dict(SeqIO.parse(open(coils_name, 'rU'), 'fasta'))
        outlist = []
        for key, value in records.items():
            if key in rejects:
                maskrange = rejects[key]
                seq = str(value.seq)
                for val in maskrange:
                    seq = seq[0:val[0]-1] + ('x' * (val[1] - val[0]+1)) + seq[val[1]:]
                outlist.append('>' + key + '\n' + seq)
            else:
                outlist.append('>' + key + '\n' + str(value.seq))
        # Make new outfile
        outfile = open(coils_name, 'w')
        outfile.write('\n'.join(outlist))
        outfile.close()

