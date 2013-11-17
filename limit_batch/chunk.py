import threading
import Queue
import logging
logging.basicConfig(level=logging.INFO)

from sh import bjobs, bsub, bkill
import re
import os
import time
import ROOT


class ChunkSubmitter(threading.Thread):
    """Handle submission of chunks to LSF"""
    def __init__(self, max_jobs, n_tries=3):
        super(ChunkSubmitter, self).__init__()
        self.job_queue = Queue.Queue()
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
            time.sleep(5*60)

    def stop(self):
        self.stop_flag = True
        for job in self.running_jobs:
            logging.info("Killing job {0}".format(job))
            bkill(job)

    def process_finished_jobs(self):
        """Process the jobs that have finished. Return the number of jobs currently running"""
        lsf_running_jobs = self.get_lsf_running_jobs()
        finished_job_ids = set(self.running_jobs.keys()) - lsf_running_jobs
        for job_id in finished_job_ids:
            job = self.running_jobs.pop(job_id)
            if self.job_successful(job_id):
                logging.info("Job Successful: "+self.get_job_result_filename(job_id))
                # send finished message to return_queue
                job.return_queue.put(self.get_job_result_filename(job_id))
            else:
                job.n_attempts += 1
                if job.n_attempts > self.n_tries:
                    logging.error("Job failed {0} times, abandoning\n".format(job.n_attempts) +
                                  "Model file: {0}, Signal Strength: {1}".format(job.model_filename, job.sig_strength))
                    job.return_queue.put("FAIL")
                else:
                    logging.warn("Job failed {0} times, resubmitting\n".format(job.n_attempts) +
                                 "Model file: {0}, Signal Strength: {1}".format(job.model_filename, job.sig_strength))
                    # increment number of attempts and return to job queue
                    self.job_queue.put(job)
        return len(lsf_running_jobs)

    def get_lsf_running_jobs(self):
        """Query LSF for status of all jobs. Returns Set of job IDs"""
        status_str = bjobs().stdout
        logging.debug(status_str)
        return set(self.lsf_status_re.findall(status_str))

    def job_successful(self, job_id):
        """Check to see if output file for job_id exists and is valid"""
        # does the file exist?
        if not os.path.isfile(self.get_job_result_filename(job_id)):
            logging.warn("File {0} does not exist".format(self.get_job_result_filename(job_id)))
            return False
        # now check if it's a valid ROOT file
        rf = ROOT.TFile(self.get_job_result_filename(job_id))
        if rf.IsZombie() or rf.TestBit(ROOT.TFile.kRecovered):
            logging.warn("File {0} is corrupt".format(self.get_job_result_filename(job_id)))
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
            self.running_jobs[self.submit_job(job)] = job

    def submit_job(self, job):
        """Submit job to LSF and return the job id"""
        res = bsub(
            _in="limit_wrapper.sh {0} {1}".format(job.model_filename, job.sig_strength),
            q="1nh", o="/dev/null")
        logging.info("submitting {0}".format(res.call_args['in']))
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
