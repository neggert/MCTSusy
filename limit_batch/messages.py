"""
Messages that get passed between threads
"""

class HypoTestChunk(object):
    """Represents one piece of a hypothesis test"""
    def __init__(self, return_queue, model_filename, sig_strength):
        super(HypoTestChunk, self).__init__()
        self.return_queue = return_queue # where to send the message when the job is done
        self.model_filename = model_filename
        self.sig_strength = sig_strength
        self.n_attempts = 0 # number of attempts to run this chunk

