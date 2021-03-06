import Queue
import ROOT
from copy import deepcopy
import re
import os
from time import sleep
import logging
logger = logging.getLogger("main")

from messages import HypoTestChunk


def run_cls_jobs(job_queue, model_filename, poi, n_jobs=20, combined_path="combined", depth=0):
    """Submit n_jobs copies of the job to job_queue"""
    m1, m2 = map(int, re.search("_(\d+)_(\d+)_", model_filename).groups())
    merged_fn = "{0}/{1}_{2}_{3}.root".format(combined_path, m1, m2, round(poi, 3))
    if not os.path.isfile(merged_fn):
        # if the file is already there, use it
        res_queue = Queue.Queue()
        logger.info("Submitting jobs for {0}, {1} at {2}".format(m1, m2, poi))
        for i in xrange(n_jobs):
            job_queue.put(HypoTestChunk(res_queue, model_filename, poi))
        out_files = []
        for i in xrange(n_jobs):
            while True:
                try:
                    output_file = res_queue.get(True, 0.1)
                except Queue.Empty:
                    sleep(1)
                    continue
                break
            if output_file != "FAIL":
                out_files.append(output_file)
        merge_results(out_files, merged_fn)
    else:
        logger.info("CLs for {0}, {1} at {2} already exists".format(m1, m2, poi))
    try:
        result = get_cls(merged_fn)
        return result
    except AttributeError:
        # delete the offending file
        os.remove(merged_fn)
        # try again
        if depth < 4:
            return run_cls_jobs(job_queue, model_filename, poi, n_jobs, combined_path, depth+1)
        else:
            raise



def merge_results(files, output_file):
    """Merge a bunch of hypothesis tests into a single file"""
    logger.debug("merging {0} to {1}".format(files, output_file))
    rf_out = ROOT.TFile(output_file, "RECREATE")
    rf = ROOT.TFile(files[0], "READ")
    res = deepcopy(rf.Get("result_sig_strength"))
    rf.Close()
    for fn in files[1:]:
        rf = ROOT.TFile(fn, "READ")
        res.Add(rf.Get("result_sig_strength"))
        rf.Close()
    rf_out.cd()
    res.Write()
    rf_out.Close()
    for fn in files:
        os.remove(fn)
    return output_file


def get_cls(merged_fn):
    """Get the CLS value from the hypothesis test in merged_fn"""
    rf = ROOT.TFile(merged_fn)
    res = rf.Get("result_sig_strength")
    try:
        obs_cls = res.CLs(0)
        exp_cls_dist = res.GetExpectedPValueDist(0)
        exp_cls = exp_cls_dist.InverseCDF(0.5)
        expp1_cls = exp_cls_dist.InverseCDF(0.84)
        expm1_cls = exp_cls_dist.InverseCDF(0.16)
    except AttributeError:
        raise
    rf.Close()
    logger.debug("Obs: {0}, Expected: {1}+{2}-{3}".format(obs_cls, exp_cls, expm1_cls, expp1_cls))
    return obs_cls, exp_cls, expm1_cls, expp1_cls

if __name__ == '__main__':
    import chunk
    chunker = chunk.ChunkSubmitter(100)
    chunker.start()
    try:
        result = run_cls_jobs(chunker.job_queue,
                              "../limits/slep_templates_200_0_combined_meas_model.root",
                              0.5)
        print result
    finally:
        chunker.stop()
