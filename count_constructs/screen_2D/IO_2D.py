from dataclasses import dataclass,field
from .matching_2D import find_wMM
from Bio import SeqIO
import sys,gzip
import gc
import itertools as itt
import ray


@dataclass
class expected_construct_2D:
    gene1 : str
    gene2 : str
    gene1_sgrnaid : str
    gene2_sgrnaid : str
    seq1 : str
    seq2 : str

    def __post_init__(self):
        self.seq2 = revcomp(self.seq2)

@dataclass
class read_pair_2D:
    id1 : str
    id2 : str
    desc1 : str
    desc2 : str
    pot_guide1 : str = field(init = False) #potential guide1
    pot_guide2 : str = field(init = False) #potential guide2
    read1 : str
    read2 : str #notation should be 5' -> 3'
    start_spacer1 : int
    end_spacer1 : int
    start_spacer2 : int
    end_spacer2 :int

    def __post_init__(self):
        self.pot_guide1 = find_spacer1(self.read1,self.start_spacer1,self.end_spacer1)
        self.pot_guide2 = find_spacer2(self.read2,self.start_spacer2,self.end_spacer2)


def find_spacer1(read : str,start : int, end : int) -> str:
    '''Attempts to search for spacer1, defaults to window around expected position.'''
    constant_part = "TACCCCTACCAACTGGTCGGGGTTTGAAAC"
    constant_part_len = len(constant_part)
    search_res = read.find(constant_part)
    read_len = len(read)

    if search_res != -1:
        if (search_res + constant_part_len + 23) <= read_len:
            return read[(search_res + constant_part_len):(search_res + constant_part_len + 23)]
        else:
            return read[start:end]
    else:
        search_res_wMM = find_wMM(template = read, query = constant_part,mismatches_allowed = 1) #search with mismatch 
        if search_res_wMM != -1:
            if (search_res_wMM + constant_part_len + 23) <= read_len:
                return read[(search_res_wMM + 30):(search_res_wMM + constant_part_len + 23)]
            else:
                return read[start:end]
        else:
            return read[start:end]


def find_spacer2(read : str,start : int, end : int) -> str:
    '''Attempts to search for spacer2, defaults to window around expected position.'''
    constant_part = "ATGAAAAAAA"
    constant_part_len = len(constant_part)
    search_res = read.find(constant_part)
    read_len = len(read)

    if search_res != -1:
        if (search_res + 10 + 23) <= read_len:
            return read[(search_res + constant_part_len):(search_res + constant_part_len + 23)]
        else:
            return read[start:end]
    else:
        search_res_wMM = find_wMM(template = read, query = constant_part,mismatches_allowed = 1) #search with mismatch 
        if search_res_wMM != -1:
            if (search_res_wMM + constant_part_len + 23) <= read_len:
                return read[(search_res_wMM + constant_part_len):(search_res_wMM + constant_part_len + 23)]
            else:
                return read[start:end]
        else:
            return read[start:end]

def revcomp(seq : str) -> str:
    '''Return reverse complement.'''
    complementary = { 'A':'T', 'T':'A', 'G':'C','C':'G' }
    return ''.join(reversed([complementary[i] for i in seq]))

@ray.remote(num_cpus=1)
class Writer:

    def __init__(self,samplefolder,resfile):
        self.samplefolder = samplefolder
        self.resfile = resfile
        self._num_pending_requests = 0

    def ping(self) -> str:
        return "Running."
    
    def get_request_queue_len(self) -> int:
        return self._num_pending_requests

    def write_discarded_reads_to_file(self,discarded_reads):
        self._num_pending_requests += 1

        outfile_fwd = self.samplefolder + "discarded_reads_forward.fastq.gz"
        outfile_rev = self.samplefolder + "discarded_reads_reverse.fastq.gz"

        with gzip.open(outfile_fwd,"at") as forward_handle,gzip.open(outfile_rev,"at") as reverse_handle:
            for curr_fwd,curr_rev in discarded_reads:
                SeqIO.write(curr_fwd, forward_handle, "fastq")
                SeqIO.write(curr_rev, reverse_handle, "fastq")

        del discarded_reads
        gc.collect()

        self._num_pending_requests -= 1

    def write_count_results_to_file(self,specific_counts_res,expected_constructs):
        self._num_pending_requests += 1

        if self.resfile == "STDOUT":
            print("gene1\tgene2\tgene1_sgrnaid\tgene2_sgrnaid\tseq1\tseq2\tcount_0MM_0MM\tcount_0MM_1MM\tcount_1MM_0MM\tcount_1MM_1MM", file = sys.stdout)
            for cnstr in expected_constructs:
                key = cnstr.seq1 + ";" + cnstr.seq2

                out_first_part = f"{cnstr.gene1}\t{cnstr.gene2}\t{cnstr.gene1_sgrnaid}\t{cnstr.gene2_sgrnaid}\t{cnstr.seq1}\t{revcomp(cnstr.seq2)}\t"
                out_second_part = f'{specific_counts_res[key]["count_0MM_0MM"]}\t{specific_counts_res[key]["count_0MM_1MM"]}\t{specific_counts_res[key]["count_1MM_0MM"]}\t{specific_counts_res[key]["count_1MM_1MM"]}'

                print(out_first_part + out_second_part, file = sys.stdout)
        else:
            with open(self.samplefolder + self.resfile,"w+") as handle:
                handle.write("gene1\tgene2\tgene1_sgrnaid\tgene2_sgrnaid\tseq1\tseq2\tcount_0MM_0MM\tcount_0MM_1MM\tcount_1MM_0MM\tcount_1MM_1MM\n")
                for cnstr in expected_constructs:
                    key = cnstr.seq1 + ";" + cnstr.seq2

                    out_first_part = f"{cnstr.gene1}\t{cnstr.gene2}\t{cnstr.gene1_sgrnaid}\t{cnstr.gene2_sgrnaid}\t{cnstr.seq1}\t{revcomp(cnstr.seq2)}\t"
                    try:
                        out_second_part = f'{specific_counts_res[key]["count_0MM_0MM"]}\t{specific_counts_res[key]["count_0MM_1MM"]}\t{specific_counts_res[key]["count_1MM_0MM"]}\t{specific_counts_res[key]["count_1MM_1MM"]}\n'
                    except:
                        out_second_part = '0\t0\t0\t0\n'

                    handle.write(out_first_part + out_second_part)
        
        self._num_pending_requests -= 1
