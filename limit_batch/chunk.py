import threading
import Queue
import logging
logger = logging.getLogger("main")

from sh import bjobs, bsub, bkill, ErrorReturnCode
import re
import os
import time
import ROOT


class ChunkSubmitter(threading.Thread):
    """Handle submission of chunks to LSF"""
    def __init__(self, lsf_queue, max_jobs, n_tries=3):
        super(ChunkSubmitter, self).__init__()
        self.job_queue = Queue.Queue()
        self.lsf_queue = lsf_queue
        self.max_jobs = max_jobs
        self.n_tries = n_tries
        self.running_jobs = {}  # key is lsf id, value is HypoTestChunk object
        self.lsf_status_re = re.compile("^(\d*?)\s", re.MULTILINE)
        self.lsf_submit_re = re.compile("^Job <(\d*)>")
        self.stop_flag = False

    def run(self):
        while True:
            if self.stop_flag:
                break
            n_jobs_running = self.process_finished_jobs()
            if n_jobs_running < self.max_jobs:
                self.submit_jobs(self.max_jobs - n_jobs_running)
            for i in xrange(60):
                if self.stop_flag:
                    break
                time.sleep(1)

    def stop(self):
        self.stop_flag = True
        for job in self.running_jobs.keys():
            logger.info("Killing job {0}".format(job))
            try:
                bkill(job)
            except:
                continue

    def process_finished_jobs(self):
        """Process the jobs that have finished. Return the number of jobs currently running"""
        lsf_running_jobs = self.get_lsf_running_jobs()
        finished_job_ids = set(self.running_jobs.keys()) - lsf_running_jobs
        for job_id in finished_job_ids:
            job = self.running_jobs.pop(job_id)
            if self.job_successful(job_id):
                logger.debug("Job Successful: "+self.get_job_result_filename(job_id))
                # send finished message to return_queue
                job.return_queue.put(self.get_job_result_filename(job_id))
            else:
                job.n_attempts += 1
                if job.n_attempts > self.n_tries:
                    logger.error("Job failed {0} times, abandoning\n".format(job.n_attempts) +
                                  "Model file: {0}, Signal Strength: {1}".format(job.model_filename, job.sig_strength))
                    job.return_queue.put("FAIL")
                else:
                    logger.warn("Job failed {0} times, resubmitting\n".format(job.n_attempts) +
                                 "Model file: {0}, Signal Strength: {1}".format(job.model_filename, job.sig_strength))
                    # increment number of attempts and return to job queue
                    self.job_queue.put(job)
        return len(lsf_running_jobs)

    def get_lsf_running_jobs(self, depth=0):
        """Query LSF for status of all jobs. Returns Set of job IDs"""
        try:
            status_str = bjobs().stdout
        except ErrorReturnCode, e:
            if depth < 4:
                logger.warn("bjobs failed with message {0}".format(str(e)))
                return self.get_lsf_running_jobs(depth+1)
            else:
                raise
        logger.debug(status_str)
        return set(self.lsf_status_re.findall(status_str))

    def job_successful(self, job_id):
        """Check to see if output file for job_id exists and is valid"""
        # does the file exist?
        if not os.path.isfile(self.get_job_result_filename(job_id)):
            logger.warn("File {0} does not exist".format(self.get_job_result_filename(job_id)))
            return False
        # now check if it's a valid ROOT file
        rf = ROOT.TFile(self.get_job_result_filename(job_id))
        if rf.IsZombie() or rf.TestBit(ROOT.TFile.kRecovered):
            logger.warn("File {0} is corrupt".format(self.get_job_result_filename(job_id)))
            return False
        rf.Close()
        return True

    def get_job_result_filename(self, job_id):
        return "tmp/{0}.root".format(job_id)  # todo

    def submit_jobs(self, n):
        """
        Submit the next n jobs in the job_queue. If there aren'try:
        enough jobs, submit what's there and return
        """
        for i in xrange(n):
            try:
                job = self.job_queue.get_nowait()
            except Queue.Empty:
                return
            try:
                submitted = self.submit_job(job)
                self.running_jobs[submitted] = job
            except ErrorReturnCode, e:
                logger.warning("Failed to submit job {0} {1}\n{2}".format(job.model_filename, job.sig_strength, str(e)))
                continue

    def submit_job(self, job):
        """Submit job to LSF and return the job id"""
        res = bsub(
            _in="limit_wrapper.sh {0} {1}".format(job.model_filename, job.sig_strength),
            q=self.lsf_queue, o="/dev/null")
        logger.debug("submitting {0}".format(res.call_args['in']))
        return self.lsf_submit_re.match(res.stdout).groups()[0]


if __name__ == '__main__':
    import messages
    jobs = Queue.Queue()
    r = Queue.Queue()
    jobs.put(messages.HypoTestChunk(r, "../limits/slep_templates_200_0_combined_meas_model.root", 0.5))
    jobs.put(messages.HypoTestChunk(r, "../limits/slep_templates_200_0_combined_meas_model.root", 0.5))
    jobs.put(messages.HypoTestChunk(r, "../limits/slep_templates_200_0_combined_meas_model.root", 1.0))
    jobs.put(messages.HypoTestChunk(r, "../limits/slep_templates_300_0_combined_meas_model.root", 1.0))
    jobs.put(messages.HypoTestChunk(r, "../limits/blahblah.root", 1.0))

    chunk_handler = ChunkSubmitter(jobs, 30)
    chunk_handler.start()
    try:
        for i in xrange(4):
            while True:
                try:
                    print r.get(True, 1)
                except Queue.Empty:
                    continue
                break
    except KeyboardInterrupt:
        chunk_handler.stop()
