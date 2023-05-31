from dataclasses import dataclass
import ray
import gc

@ray.remote(num_cpus=1,max_restarts=0)
@dataclass
class Global_counts_2D:
    success_0MM : int = 0
    success_1MM : int = 0
    discarded_count : int = 0
    total_reads_processed : int = 0
    _num_pending_requests : int = 0

    def batch_increase(self,cdict : dict[str,int] ):
        self._num_pending_requests += 1
        self._increase_success_0MM(by = cdict["success_0MM"])
        self._increase_success_1MM(by = cdict["success_1MM"])
        self._increase_discarded_count(by = cdict["discarded_count"])

        del cdict
        gc.collect() #manual call to GC
        self._num_pending_requests -= 1

    def _increase_success_0MM(self,by : int):
        self.success_0MM += by
        self._increase_total_reads_processed(2 * by)

    def _increase_success_1MM(self,by : int):
        self.success_1MM += by
        self._increase_total_reads_processed(2 * by)

    def _increase_discarded_count(self, by :int):
        self.discarded_count += by
        self._increase_total_reads_processed(2 * by)

    def _increase_total_reads_processed(self, by :int):
        self.total_reads_processed += by

    def ping(self) -> str:
        return "Running."
    
    def get_request_queue_len(self) -> int:
        return self._num_pending_requests

    def get_counter(self):
        return {"success_0MM": self.success_0MM , "success_1MM" : self.success_1MM, "discarded_count" : self.discarded_count, "total_reads_processed" : self.total_reads_processed}


@ray.remote(num_cpus=1,max_restarts=0)
class Specific_counts_2D:

    def __init__(self,expected_constructs):
        self._num_pending_requests = 0
        self.lookup = dict()
        for o in expected_constructs:
            self.lookup.update({f"{o.seq1};{o.seq2}":{"count_0MM_0MM" : 0,"count_0MM_1MM" :0, "count_1MM_0MM" :0, "count_1MM_1MM":0}})

            
    def ping(self) -> str:
        return "Running."
    
    def batch_increase(self, cdict : dict[str,dict[str,int]]):
        self._num_pending_requests += 1
        for k1 in cdict.keys():
            for k2 in cdict[k1].keys():
                self.lookup[k1][k2] += cdict[k1][k2]

        del cdict 
        gc.collect() #manual call to GC
        self._num_pending_requests -= 1

    def get_counter(self):
        return self.lookup
    
    def get_request_queue_len(self) -> int:
        return self._num_pending_requests
    
    def send_to_writer(self,writer,expected_constructs_ref):
        self._num_pending_requests += 1
        expected_constructs_ref = expected_constructs_ref[0]
        ray.get(writer.write_count_results_to_file.remote(specific_counts_res = self.get_counter(),expected_constructs = expected_constructs_ref))
        self._num_pending_requests -= 1

