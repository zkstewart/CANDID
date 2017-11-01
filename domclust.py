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
            print(newseq)
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
        if args['leaf'] == 'y':
            clusterer = hdbscan.HDBSCAN(metric='precomputed', cluster_selection_method = 'leaf', min_cluster_size = int(args['minsize']), min_samples = int(args['minsample']))
        else:
            clusterer = hdbscan.HDBSCAN(metric='precomputed', min_cluster_size = int(args['minsize']), min_samples = int(args['minsample']))
        clusterer.fit(matrix1.data)
        new_groups = clusterer.labels_
        # Relabel group
        #for x in range(len(new_groups)):
        #    if new_groups[x] != -1:
        #        new_groups[x] = new_groups[x] + len(group_dict)
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
            shutil.rmtree(align_out_dir, ignore_errors=False, onerror=None)
            os.mkdir(align_out_dir)
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

    def clust_update(args, outdir, basename, group_dict, rejects_list):
        import os, re
        from Bio import SeqIO
        print('Updating clusters...')
        # Read in hmmer output file
        align_out_dir = os.path.join(os.getcwd(), outdir, 'tmp_alignments')
        hmmer_results = os.path.join(align_out_dir, 'tmp_hmmer.results')
        hmmer_file = open(hmmer_results, 'r')
        # Parse contents and compare against current domain groups
        #hmm_reg = re.compile(r'Domain_(\d{0,10})_align')
        ongoingCount = 0
        iterate = 'n'
        for line in hmmer_file:
            if line.startswith('#') or line == '' or line == '\n': continue
            split_line = line.split()
            seqid = split_line[0]
            seqrange = set(range(int(split_line[17])-1, int(split_line[18])))
            evalue = float(split_line[12])
            if evalue > float(args['hmmevalnov']):
                continue
            # Check if this is already part of a group
            overlap = 'n'
            seqnum = 1      # Every time we find another domain region with this same seqid, we add +1 here so we can properly label it if we add it into our unclusted_domains.fasta file
            for key, value in group_dict.items():
                if overlap == 'y':
                    break
                for val in value:
                    val_split = val.split('_')
                    valid = '_'.join(val_split[0:-3])           # We split off the range value here (e.g. 1-72) and just get the sequence ID
                    if valid != seqid:
                        continue
                    seqnum += 1
                    # Check for overlap [if we get here, then the seqid is in another group]
                    valrange = val_split[-1].split('-')
                    valrange = set(range(int(valrange[0])-1, int(valrange[1])))
                    if len(seqrange-valrange) / len(seqrange) < 0.9:    # i.e., if more than 10% of seqrange (new sequence) overlaps valrange (established sequence), we discard it since we want to minimise overlapping sequences. Decisions like these make this script designed to find the locations of new domain regions, but it isn't necessarily capable of finding the exact borders of domains.
                        overlap = 'y'
                        break
            if overlap == 'y':
                continue
            # Check if this region hasn't previously been recognised and rejected [if we get here, then the sequence does not significantly overlap any sequences that are currently clustered]
            overlap = 'n'
            for reject in rejects_list:
                reject_split = reject.split('_')
                rejectid = '_'.join(reject_split[0:-3])
                if rejectid != seqid:
                    continue
                seqnum += 1
                rejectrange = reject_split[-1].split('-')               # If we get here, this sequence has previously been recognised and rejected - we just need to see if it is the same range
                rejectrange = set(range(int(rejectrange[0])-1, int(rejectrange[1])))
                if len(seqrange-rejectrange) / len(seqrange) < 0.9:
                    overlap = 'y'
                    break
            if overlap == 'y':
                continue
            # Check if this sequence passes our length cutoff criteria
            infile = os.path.join(os.getcwd(), outdir, basename + '_cdhit.fasta')
            records = SeqIO.parse(open(infile, 'rU'), 'fasta')
            for record in records:
                if record.id == seqid:
                    aaseq = str(record.seq)
                    break
            aaseq = aaseq[int(split_line[17])-1:int(split_line[18])]
            if not len(aaseq) >= args['cleanAA']:
                continue
            # Add any sequence that gets this deep in our loop into our unclustered domains group for a further iteration
            iterate = 'y'
            newseqname = '>' + seqid + '_Domain_' + str(seqnum) + '_' + split_line[17] + '-' + split_line[18]
            outfile = open(os.path.join(os.getcwd(), outdir, basename + '_unclustered_domains.fasta'), 'a')
            outfile.write(newseqname + '\n' + aaseq + '\n')     # The unclustered_domains.fasta should always end on a newline, so we don't need to preface our addition with '\n', but we do finish it with one
            outfile.close()
            # Add the new sequence into the group_dict for any further loops
            ongoingCount += 1
            group_dict['irrelevant_' + str(ongoingCount)] = [newseqname[1:]]

        # Let user know if this script is going to go through another iterative round
        if iterate == 'n':
            print('No more domain regions were identified with iterative HMMER searching.')
        else:
            print('New potential domain region(s) were identified with iterative HMMER searching...')
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
