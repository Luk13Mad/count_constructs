import click
import sys,logging,os,gzip,gc,time
import itertools as itt
from Bio import SeqIO
import ray
from .counter_2D import Global_counts_2D,Specific_counts_2D
from .matching_2D import find_wMM
from . import IO_2D



@click.command()
@click.option("--forward","-f",type = click.Path(), required = True, help = "Path to FASTQ.gz file with forward reads.")
@click.option("--reverse","-r",type = click.Path(), required = True, help = "Path to FASTQ.gz file with reverse reads.")
@click.option("--expected","-e",type = click.Path(), required = True, help = "Path to TSV file with expected constructs.")
@click.option("--samplefolder",type = click.Path(exists = True), required = True, help = "Absolute path to your sample folder.")
@click.option("--startspacer1","-s1",type = click.INT, required = True, help = "Start default extraction window first spacer.")
@click.option("--startspacer2","-s2",type = click.INT, required = True, help = "Start default extraction window second spacer.")
@click.option("--endspacer1","-e1",type = click.INT, required = True, help = "End default extraction window first spacer.")
@click.option("--endspacer2","-e2",type = click.INT, required = True, help = "End default extraction window second spacer.")
@click.option("--resfile", default = "count_results.tsv", help = "Name of results file. STDOUT prints results to standard out.")
@click.option("--workers","-w",type = click.INT,default = 1,help = "Max number of workers to be used.")
@click.option("--batchsize",type = click.INT, default = 200000, help = "Batchsize for processing FASTQ pairs.")
@click.option("--rayobjectstore",type = click.INT, default = 10, help = "Memory size in GB for ray object store.")
def setup_2D(expected,samplefolder,forward,reverse,batchsize,workers,resfile,startspacer1,startspacer2,endspacer1,endspacer2,rayobjectstore):

    #check if samplefolder is absolute path
    if not os.path.isabs(samplefolder):
        logging.error("Sample folder must be supplied as absolute path.")
        sys.exit(1)

    #check if samplefolder path ends with separation
    if not samplefolder.endswith(os.sep):
        samplefolder = os.path.join(samplefolder, "")

    os.chdir(samplefolder)

    #init logger
    logging.basicConfig(filename=samplefolder + 'count_constructs.log', filemode='w', encoding='utf-8', level=logging.INFO, format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    
    logging.info(f"Initializing ray with number of CPU's: {workers}")
    #setting number of ray workers
    ray.init(num_cpus = (workers - 3),  # subtract 3 to reserve those cores for the counters and writer
        object_store_memory = (rayobjectstore * (10**9)))

    logging.info(f"Forward FASTQ from : {forward}")
    logging.info(f"Reverse FASTQ from : {reverse}")
    logging.info(f"Expected constructs from : {expected}")
    logging.info(f"Batchsize : {batchsize}")
    logging.info(f"Object store memory size : {rayobjectstore} GB")

    #put some constants in ray object store
    startspacer1_ref,endspacer1_ref,startspacer2_ref,endspacer2_ref = ray.put(startspacer1),ray.put(endspacer1),ray.put(startspacer2),ray.put(endspacer2)

    #load expected constructs from file
    logging.info("Started reading expected constructs.")
    expected_constructs = []

    with open(expected,"r") as expected_handle:
        try:
            for line in expected_handle:
                gene1,gene2,gene1_sgrnaid,gene2_sgrnaid,seq1,seq2 = line.strip().split("\t")
                expected_constructs.append(IO_2D.expected_construct_2D(gene1,gene2,gene1_sgrnaid,gene2_sgrnaid,seq1,seq2))
        except Exception as e:
            logging.error(f"Encountered error while reading expected construct: {e}")
            sys.exit(1) 

    expected_constructs_ref = ray.put(expected_constructs)

    logging.info(f"Finished reading expected. {len(expected_constructs)} constructs read.")

    #dict with seq1 as key for faster lookup
    lookup_seq1 = {}
    for cnstr in expected_constructs:
        if cnstr.seq1 in lookup_seq1:
            lookup_seq1[cnstr.seq1].append(cnstr.seq2)
        else:
            lookup_seq1.update({cnstr.seq1 : [cnstr.seq2] })
    lookup_seq1_ref = ray.put(lookup_seq1)
    del lookup_seq1

    #removing discarded_reads output file if it already exists from previous run
    discarded_reads_outfile_fwd = samplefolder + "discarded_reads_forward.fastq.gz"
    discarded_reads_outfile_rev = samplefolder + "discarded_reads_reverse.fastq.gz"
    if os.path.isfile(discarded_reads_outfile_fwd):
        logging.warning(f"File already exists: {discarded_reads_outfile_fwd}")
        try:
            logging.warning(f"Try deleting file: {discarded_reads_outfile_fwd}")
            os.remove(discarded_reads_outfile_fwd)
        except:
            logging.error(f"Failed deleting: {discarded_reads_outfile_fwd}")
            sys.exit(1)

    if os.path.isfile(discarded_reads_outfile_rev):
        logging.warning(f"File already exists: {discarded_reads_outfile_rev}")
        try:
            logging.warning(f"Try deleting file: {discarded_reads_outfile_rev}")
            os.remove(discarded_reads_outfile_rev)
        except:
            logging.error(f"Failed deleting: {discarded_reads_outfile_rev}")
            sys.exit(1)

    #init counters and writer as ray actors
    global_counter_ref = Global_counts_2D.remote()
    specific_counter_ref = Specific_counts_2D.remote(expected_constructs_ref)
    writer_ref = IO_2D.Writer.remote(samplefolder = samplefolder,resfile = resfile)

    logging.info(f"Global counter status: {ray.get(global_counter_ref.ping.remote())}")
    logging.info(f"Specific counter status: {ray.get(specific_counter_ref.ping.remote())}")
    logging.info(f"Writer status: {ray.get(writer_ref.ping.remote())}")
    logging.info("Start processing....")

    #list which will hold ray object IDs
    ray_tasks = list()

    #loop distributing tasks with back pressure
    with gzip.open(forward,"rt") as fwd_handle,gzip.open(reverse,"rt") as rev_handle:
        fwd_seqreader = SeqIO.parse(fwd_handle, "fastq")
        rev_seqreader = SeqIO.parse(rev_handle, "fastq")
        for idx,read_batch in enumerate(batched(iterable = zip(fwd_seqreader,rev_seqreader), n = batchsize)):
            if len(ray_tasks) > (workers):
                logging.info("Max number of processes reached, waiting for some to finish.")
                done_ids,ray_tasks = ray.wait(ray_tasks)
                logging.info(f"Task(s) finished: {len(done_ids)}.")
                del done_ids

            logging.info(f"Distributed batch number: {idx}")
            read_batch_ref = ray.put(read_batch)
            #read_batch will be list with length n, containing tuples (forward_read,reverse_read)
            ray_tasks.append(main_loop_2D.remote(read_batch = read_batch_ref,lookup_seq1 = lookup_seq1_ref,
            global_counter = global_counter_ref,specific_counter = specific_counter_ref, writer = writer_ref,
            start_spacer1 = startspacer1_ref, end_spacer1 = endspacer1_ref,
            start_spacer2 = startspacer2_ref, end_spacer2 = endspacer2_ref))

            del read_batch_ref,read_batch
            gc.collect()

    logging.info("All batches distributed, waiting for results.")


    #wait for all tasks to finish before further processing
    # write discarded reads to file
    while len(ray_tasks):
        done_ids,ray_tasks = ray.wait(ray_tasks)
        logging.info(f"Task finished.")
        del done_ids
    gc.collect()

    logging.info("Processing finished.")

    #print global counter to file
    while ray.get(global_counter_ref.get_request_queue_len.remote()) != 0:
        time.sleep(60)
    
    logging.info("Waiting on global count results.")
    global_counter_res = ray.get(global_counter_ref.get_counter.remote())
    logging.info(f'success_0MM : {global_counter_res["success_0MM"]}, success_1MM: {global_counter_res["success_1MM"]}, discarded_count: {global_counter_res["discarded_count"]},total_reads_processed: {global_counter_res["total_reads_processed"]}')

    #print specific counter to file
    while ray.get(specific_counter_ref.get_request_queue_len.remote()) != 0:
        time.sleep(60)

    logging.info("Writing specific count results to file.")
    ray.get(specific_counter_ref.send_to_writer.remote(writer_ref,[expected_constructs_ref]))

    while ray.get(writer_ref.get_request_queue_len.remote()) != 0:
        time.sleep(60)

    logging.info("Shutting down ray.")
    ray.shutdown()

    logging.info("Exiting.")
    sys.exit(0)

        

@ray.remote
def main_loop_2D(lookup_seq1,global_counter,specific_counter,read_batch,start_spacer1,end_spacer1,start_spacer2,end_spacer2,writer):

    read_pairs = [IO_2D.read_pair_2D(id1 = f.id, id2 = r.id,
    desc1 = f.description, desc2 = r.description,
    read1 = str(f.seq),read2 = str(r.seq),
    start_spacer1 = start_spacer1, end_spacer1 = end_spacer1,
    start_spacer2 = start_spacer2, end_spacer2 = end_spacer2) for f,r in read_batch]

    #help variables
    seq1 = set(lookup_seq1.keys())

    #read pairs which will be discarded
    read_pairs_to_discard = []

    #calls to counter are better if batched

    #global counter local
    gcounter_local = {"success_0MM": 0 , "success_1MM" : 0, "discarded_count" : 0}

    #specific counter local
    scounter_local = {}

    for idx,current_readpair in enumerate(read_pairs):
        if current_readpair.pot_guide1 in lookup_seq1 : #guide1 was found with 0MM
            if current_readpair.pot_guide2 in lookup_seq1[current_readpair.pot_guide1]:
                #both guides found with 0MM
                key = current_readpair.pot_guide1 + ";" + current_readpair.pot_guide2
                try:
                    scounter_local[key]["count_0MM_0MM"] += 1
                except:
                    scounter_local.update({key : {"count_0MM_0MM" : 1,"count_0MM_1MM" :0, "count_1MM_0MM" :0, "count_1MM_1MM":0}})
                gcounter_local["success_0MM"] += 1
            else:
                # search pot_guide2 with 1MM
                res_guide2_1MM = list(itt.filterfalse(lambda x: find_wMM(current_readpair.pot_guide2,x,1) == -1,lookup_seq1[current_readpair.pot_guide1]))
                if len(res_guide2_1MM) == 1:
                    #found guide1_0MM and guide2_1MM
                    tmp_guide2_1MM = res_guide2_1MM[0]
                    key = current_readpair.pot_guide1 + ";" + tmp_guide2_1MM
                    try:
                        scounter_local[key]["count_0MM_1MM"] += 1
                    except:
                        scounter_local.update({key : {"count_0MM_0MM" : 0,"count_0MM_1MM" :1, "count_1MM_0MM" :0, "count_1MM_1MM":0}})
                    gcounter_local["success_1MM"] += 1
                elif len(res_guide2_1MM) > 1: #if multiple guides 2 found discard
                    gcounter_local["discarded_count"] += 1
                    read_pairs_to_discard.append(idx)
                else: #did not find pot_guide2 with 1MM
                    #try full search for guide2
                    res_guide2_full = list(itt.filterfalse(lambda x: find_wMM(current_readpair.read2,x,1) == -1,lookup_seq1[current_readpair.pot_guide1]))
                    if len(res_guide2_full) == 1:
                        tmp_guide2_full = res_guide2_full[0]
                        key = current_readpair.pot_guide1 + ";" + tmp_guide2_full
                        try:
                            scounter_local[key]["count_0MM_1MM"] += 1
                        except:
                            scounter_local.update({key : {"count_0MM_0MM" : 0,"count_0MM_1MM" :1, "count_1MM_0MM" :0, "count_1MM_1MM":0}})
                        gcounter_local["success_1MM"] += 1
                    else:
                        gcounter_local["discarded_count"] += 1
                        read_pairs_to_discard.append(idx)
        else:
            #try finding current_readpair.pot_guide1 with one mismatch
            res_guide1_1MM = list(itt.filterfalse(lambda x: find_wMM(current_readpair.pot_guide1,x,1) == -1,seq1))
            if len(res_guide1_1MM) == 1: #pot_guide1 was found with 1MM
                tmp_guide1_1MM = res_guide1_1MM[0]
                if current_readpair.pot_guide2 in lookup_seq1[tmp_guide1_1MM]:
                    #guide1 found with 1MM guide2 found with 0MM
                    key = tmp_guide1_1MM + ";" + current_readpair.pot_guide2
                    try:
                        scounter_local[key]["count_1MM_0MM"] += 1
                    except:
                        scounter_local.update({key : {"count_0MM_0MM" : 0,"count_0MM_1MM" :0, "count_1MM_0MM" :1, "count_1MM_1MM":0}})
                    gcounter_local["success_1MM"] += 1
                else:
                    # search pot_guide2 with 1MM
                    res_guide2_1MM = list(itt.filterfalse(lambda x: find_wMM(current_readpair.pot_guide2,x,1) == -1,lookup_seq1[tmp_guide1_1MM]))
                    if len(res_guide2_1MM) == 1:
                        #found guide1_1MM and guide2_1MM
                        tmp_guide2_1MM = res_guide2_1MM[0]
                        key = tmp_guide1_1MM + ";" + tmp_guide2_1MM
                        try:
                            scounter_local[key]["count_1MM_1MM"] += 1
                        except:
                            scounter_local.update({key : {"count_0MM_0MM" : 0,"count_0MM_1MM" :0, "count_1MM_0MM" :0, "count_1MM_1MM":1}})
                        gcounter_local["success_1MM"] += 1
                    elif len(res_guide2_1MM) > 1: #discard if multiple matches for guide2
                        gcounter_local["discarded_count"] += 1
                        read_pairs_to_discard.append(idx)
                    else: #did not find guide in pot_guide2, do full search
                        res_guide2_full = list(itt.filterfalse(lambda x: find_wMM(current_readpair.read2,x,1) == -1,lookup_seq1[tmp_guide1_1MM]))
                        if len(res_guide2_full) == 1:
                            tmp_guide2_full = res_guide2_full[0]
                            key = tmp_guide1_1MM + ";" + tmp_guide2_full
                            try:
                                scounter_local[key]["count_1MM_1MM"] += 1
                            except:
                                scounter_local.update({key : {"count_0MM_0MM" : 0,"count_0MM_1MM" :0, "count_1MM_0MM" :0, "count_1MM_1MM":1}})
                            gcounter_local["success_1MM"] += 1
                        else: #discard if full search failed
                            gcounter_local["discarded_count"] += 1
                            read_pairs_to_discard.append(idx)
            elif len(res_guide1_1MM) > 1: #discard if multiple guide1 found
                gcounter_local["discarded_count"] += 1
                read_pairs_to_discard.append(idx)
            else: #pot_guide1 not found with 1MM so do full search
                res_guide1_full = list(itt.filterfalse(lambda x: find_wMM(current_readpair.read1,x,1) == -1,seq1))
                if len(res_guide1_full) == 1: #if only one guide1 was found with full search
                    tmp_guide1_full = res_guide1_full[0]
                    if current_readpair.pot_guide2 in lookup_seq1[tmp_guide1_full]:
                    #guide1 found with 1MM guide2 found with 0MM
                        key = tmp_guide1_full + ";" + current_readpair.pot_guide2
                        try:
                            scounter_local[key]["count_1MM_0MM"] += 1
                        except:
                            scounter_local.update({key : {"count_0MM_0MM" : 0,"count_0MM_1MM" :0, "count_1MM_0MM" :1, "count_1MM_1MM":0}})
                        gcounter_local["success_1MM"] += 1
                    else:
                    # search pot_guide2 with 1MM
                        res_guide2_1MM = list(itt.filterfalse(lambda x: find_wMM(current_readpair.pot_guide2,x,1) == -1,lookup_seq1[tmp_guide1_full]))
                        if len(res_guide2_1MM) == 1:
                            #found guide1_1MM and guide2_1MM
                            tmp_guide2_1MM = res_guide2_1MM[0]
                            key = tmp_guide1_full + ";" + tmp_guide2_1MM
                            try:
                                scounter_local[key]["count_1MM_1MM"] += 1
                            except:
                                scounter_local.update({key : {"count_0MM_0MM" : 0,"count_0MM_1MM" :0, "count_1MM_0MM" :0, "count_1MM_1MM":1}})
                            gcounter_local["success_1MM"] += 1
                        elif len(res_guide2_1MM) > 1: #discard if multiple matches for guide2
                            gcounter_local["discarded_count"] += 1
                            read_pairs_to_discard.append(idx)
                        else: #did not find pot_guide2 or found multiple
                            res_guide2_full = list(itt.filterfalse(lambda x: find_wMM(current_readpair.read2,x,1) == -1,lookup_seq1[tmp_guide1_full]))
                            if len(res_guide2_full) == 1:
                                tmp_guide2_full = res_guide2_full[0]
                                key = tmp_guide1_full + ";" + tmp_guide2_full
                                try:
                                    scounter_local[key]["count_1MM_1MM"] += 1
                                except:
                                    scounter_local.update({key : {"count_0MM_0MM" : 0,"count_0MM_1MM" :0, "count_1MM_0MM" :0, "count_1MM_1MM":1}})
                                gcounter_local["success_1MM"] += 1
                            else: #discard if full search failed
                                gcounter_local["discarded_count"] += 1
                                read_pairs_to_discard.append(idx)
                else:
                    gcounter_local["discarded_count"] += 1
                    read_pairs_to_discard.append(idx)

    #batched call to increase counter
    #putting local counter into ray object storage
    gcounter_local_ref = ray.put(gcounter_local)
    scounter_local_ref = ray.put(scounter_local)

    #forcing ray to wait for counters to process call
    global_counter.batch_increase.remote(gcounter_local_ref)
    specific_counter.batch_increase.remote(scounter_local_ref)

    discarded_reads_for_writing = [read_batch[x] for x in read_pairs_to_discard]
    discarded_reads_for_writing_ref = ray.put(discarded_reads_for_writing)
    ray.get(writer.write_discarded_reads_to_file.remote(discarded_reads_for_writing_ref))

    del gcounter_local,scounter_local,read_pairs,seq1,gcounter_local_ref,scounter_local_ref,read_batch,read_pairs_to_discard,discarded_reads_for_writing


    

def batched(iterable, n):
    "Batch data into lists of length n. The last batch may be shorter. From itertools recipes."
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError('n must be at least one')
    while (batch := list(itt.islice(iterable, n))):
        yield batch
