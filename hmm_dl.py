class hmm_dl:
    def fasta_dl(args, filename):
        # Function adapted from https://stackoverflow.com/questions/32288113/python-3-how-to-create-a-text-progress-bar-for-downloading-files
        import urllib.request, os, traceback
        print('Downloading CDD ' + filename + ' file...')
        url = args.url
        u = urllib.request.urlopen(url)
        meta = u.info()
        #print(str(meta).split())
        metaInfo = str(meta).split()
        #print(len(metaInfo))
        print ("Download size: " + str(round(int(metaInfo[3])/1000000, 1)) + " megabytes")
        fileTotalbytes=int(metaInfo[3])

        data_blocks = []
        total=0
        # Put rest of function in try:except clause to delete the file if the download fails
        try:
            # Open output file
            out_name = os.path.join(args.outdir, filename)
            f = open(out_name, "wb")
            # Download loop
            while True:
                block = u.read(1024)
                data_blocks.append(block)
                total += len(block)
                hash = ((60*total)//fileTotalbytes)
                print("[{}{}] {}%".format('#' * hash, ' ' * (60-hash), int(total/fileTotalbytes*100)), end="\r")

                # Save to file as we go to reduce memory usage
                if len(data_blocks) == 10000:   # Saves every 10mb
                    f.write(b''.join(data_blocks))
                    data_blocks = []

                # Break out of dl loop when complete
                if not len(block):
                    break
            
            # Final save
            f.write(b''.join(data_blocks))
            f.close()
            u.close()
        except:
            var = traceback.format_exc()
            print("Unexpected error:")
            print(var)
            f.close()
            u.close()
            os.remove(out_name)
            quit()

    def pfam_dl(args, pfam_filename):
        import urllib.request, os, traceback
        print('Downloading PFAM ' + filename + ' file...')
        url = args.pfamurl
        u = urllib.request.urlopen(url)
        meta = u.info()
        #print(str(meta).split())
        metaInfo = str(meta).split()
        #print(len(metaInfo))
        print ("Download size: " + str(round(int(metaInfo[3])/1000000, 1)) + " megabytes")
        fileTotalbytes=int(metaInfo[3])

        data_blocks = []
        total=0
        # Put rest of function in try:except clause to delete the file if the download fails
        try:
            # Open output file
            out_name = os.path.join(args.outdir, filename)
            f = open(out_name, "wb")
            # Download loop
            while True:
                block = u.read(1024)
                data_blocks.append(block)
                total += len(block)
                hash = ((60*total)//fileTotalbytes)
                print("[{}{}] {}%".format('#' * hash, ' ' * (60-hash), int(total/fileTotalbytes*100)), end="\r")

                # Save to file as we go to reduce memory usage
                if len(data_blocks) == 10000:   # Saves every 10mb
                    f.write(b''.join(data_blocks))
                    data_blocks = []

                # Break out of dl loop when complete
                if not len(block):
                    break
            
            # Final save
            f.write(b''.join(data_blocks))
            f.close()
            u.close()
        except:
            var = traceback.format_exc()
            print("Unexpected error:")
            print(var)
            f.close()
            u.close()
            os.remove(out_name)
            quit()
        
    def untar(args, filename, extractdir):
        import os, tarfile, gzip
        print('Extracting ' + filename + '...')
        # Get file details
        file = os.path.join(args.outdir, filename)
        extract_dir = os.path.join(args.outdir, extractdir)
        # Jump into output directory
        if file.endswith(".tar.gz"):
            tar = tarfile.open(file, "r:gz")
            tar.extractall(extract_dir)
            tar.close()
        elif file.endswith(".tar"):
            tar = tarfile.open(file, "r:")
            tar.extractall(extract_dir)
            tar.close()
        elif file.endswith(".tar.bz2"):
            tar = tarfile.open(file, "r:bz2")
            tar.extractall(extract_dir)
            tar.close()
        elif file.endswith(".gz"):
            tar = gzip.open(file, "rb")
            outfile = open(os.path.join(extract_dir, filename[0:-3]), 'wb')
            for line in tar:
                outfile.write(line)
            tar.close()
            outfile.close()
        else:
            print('I don\'t recognise the file format for ' + filename)
            print('This script assume it is compressed with .tar, .tar.gz, or .tar.bz2 extension')
            print('To proceed further, I recommend you uncompress this file yourself into a directory titled')
            print(extract_dir)
            print('and then re-run this script, as it will automatically detect this directory.')

    def hmmbuild(args, filename):
        # Define thread function for later use
        def build_loop(loop_args, startNum, endNum):
            import subprocess, os
            for i in range(startNum, endNum):
            # Loop through all individual msas to produce hmms
                file_prefix = '.'.join(loop_args[2][i].split('.')[0:-1])
                # Continue from previous run if cancelled
                hmmout_name = os.path.join(loop_args[1], file_prefix + '.hmm')
                if os.path.isfile(hmmout_name):
                    if os.path.getsize(hmmout_name) != 0:
                        continue
                cmd = loop_args[3] + ' ' + hmmout_name + ' ' + os.path.join(loop_args[0], loop_args[2][i])
                run_hmmbuild = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                hmmout, hmmerr = run_hmmbuild.communicate()
                # Error handler (inside loop)
                if hmmerr.decode("utf-8") != '':
                    raise Exception('hmmbuild error text below\n' + str(hmmerr.decode("utf-8")) + '\nMake sure that you define the -hmmer3dir argument if this directory is not in your PATH')
                
        # Continue defining the thread function-using function
        import os, threading
        print('Building individual HMMs from ' + filename + ' MSAs...')
        # Get file details
        extract_dir = os.path.join(args.outdir, 'cdd_extraction')
        ind_hmm_dir = os.path.join(args.outdir, 'cdd_individual_models')
        extracted_msas = os.listdir(extract_dir)
        cmd_prefix = os.path.join(args.hmmer3dir, 'hmmbuild')
        dirsize = len(extracted_msas)
        # Set up threading requirements
        threads = args.threads
        chunk_size = int(dirsize / threads)
        loop_args = [extract_dir, ind_hmm_dir, extracted_msas, cmd_prefix, dirsize]
        processing_threads = []
        # Begin the loop
        for i in range(threads):
            start = chunk_size * i
            if i+1 != int(args['threads']):
                end = chunk_size * (i+1)
            else:
                end = dirsize
            build = threading.Thread(target=build_loop, args=(loop_args, start, end))
            processing_threads.append(build)
            build.start()
            print('Initiated thread num ' + str(i+1) + ' for individual hmm building...')

        # Wait for all threads to end.
        for process_thread in processing_threads:
            process_thread.join()
        print('Individual HMM building completed')

    def concat_hmms(args):
        import os
        print('Concatenating individual HMMs into a single .hmm...')
        # Get file details
        ind_hmm_dir = os.path.join(args.outdir, 'cdd_individual_models')
        ind_hmms = os.listdir(ind_hmm_dir)
        outhmm_name = os.path.join(args.outdir, 'cdd_db', 'CDD.hmm')
        outhmm_file = open(outhmm_name, 'w')
        # Concatenate files
        for ind_hmm in ind_hmms:
            outhmm_file.write(open(os.path.join(ind_hmm_dir, ind_hmm), 'r').read())
        outhmm_file.close()

    def convert(args, extractdir):
        import subprocess, os
        print('Converting SUPERFAMILY database to HMMER 3.1 version...')
        # Get file details
        if args.superfamily.endswith('.gz'):
            supfam_name = args.superfamily[0:-3]        # If, for whatever reason, SUPERFAMILY either is compressed with something other than just gz, or if it has already been decompressed and the user is providing its uncompressed name as an argument, we need to handle its file name
        else:
            supfam_name = args.superfamily              # I assume that, if the file name doesn't end with .gz, it was already decompressed. This may cause issues if SUPERFAMILY is distributed with another type of compression in the future.
        file_in = os.path.join(args.outdir, extractdir, supfam_name)
        out_name = file_in + '_3.1'
        # Format command
        cmd_prefix = os.path.join(args.hmmer3dir, 'hmmconvert')
        cmd = os.path.join(args.hmmer3dir, 'hmmconvert') + ' ' + file_in + ' > ' + out_name
        # Run
        run_hmmconvert = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        hmmout, hmmerr = run_hmmconvert.communicate()
        # Error handler (inside loop)
        if hmmerr.decode("utf-8") != '':
            raise Exception('hmmconvert error text below\n' + str(hmmerr.decode("utf-8")) + '\nMake sure that you define the -hmmer3dir argument if this directory is not in your PATH')

    def concat_additional(args):
        import os
        print('Concatenating additional HMM databases into a single .hmm...')
        # Open CDD hmm for append
        cddhmm_name = os.path.join(args.outdir, 'cdd_db', 'CDD.hmm')
        cddhmm_file = open(cddhmm_name, 'a')
        # Concatenate SUPERFAMILY file if arguments were provided
        if args.superfamily != 'n':
            print('Concatenating SUPERFAMILY...')
            if args.superfamily.endswith('.gz'):
                supfam_name = args.superfamily[0:-3]
            else:
                supfam_name = args.superfamily
            if supfam_name.endswith('_3.1'):
                donothing = 0
            else:
                supfam_name += '_3.1'
            supfam_hmm_name = os.path.join(args.outdir, 'superfamily_extraction', supfam_name)
            supfam_hmm_file = open(supfam_hmm_name, 'r')
            for line in supfam_hmm_file:
                cddhmm_file.write(line)
            supfam_hmm_file.close()
        # Concatenate CATH file if arguments were provided
        if args.cath != 'n':
            print('Concatenating CATH...')
            if args.cath.endswith('.gz'):
                cath_name = args.cath[0:-3]
            else:
                cath_name = args.cath
            cath_hmm_name = os.path.join(args.outdir, 'cath_extraction', cath_name)
            cath_hmm_file = open(cath_hmm_name, 'r')
            for line in cath_hmm_file:
                cddhmm_file.write(line)
            cath_hmm_file.close()
        cddhmm_file.close()
            
        # Leave a note saying that the CDD file also includes SUPERFAMILY/CATH
        note_name = os.path.join(args.outdir, 'cdd_db', 'note.txt')
        note_file = open(note_name, 'w')
        if args.superfamily != 'n':
            note_file.write('This .hmm file contains SUPERFAMILY models as well\n')
        if args.superfamily != 'n':
            note_file.write('This .hmm file contains CATH models as well\n')
        note_file.close()
