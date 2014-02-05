import threading
import cls_finder
import Queue
import ROOT
from glob import glob
import re
import signal
import logging
logger = logging.getLogger('main.get_limit')

class BracketValues(object):
    def __init__(self):
        self.low_poi = 0.
        self.low_poi_cls = 1.
        self.high_poi = 1000.
        self.high_poi_cls = 0.

    def new_val(self, poi, cls):
        if cls > 0.05 and poi < self.high_poi:
            self.new_low(poi, cls)
        elif poi > self.low_poi:
            self.new_high(poi, cls)
        logger.debug("Low: {0}, High: {1}".format((self.low_poi, self.low_poi_cls), (self.high_poi, self.high_poi_cls)))
        assert(self.low_poi < self.high_poi)

    def new_low(self, poi, cls):
        if self.low_poi_cls < 0.05 and poi < self.high_poi:
            # we're still bracketing
            self.low_poi = poi
            self.low_poi_cls = cls
        elif cls < self.low_poi_cls and poi < self.high_poi:
            self.low_poi = poi
            self.low_poi_cls = cls
        elif cls == self.low_poi_cls:
            self.low_poi = max(poi, self.low_poi)

    def new_high(self, poi, cls):
        if self.high_poi_cls > 0.05 and poi > self.low_poi:
            # we're still bracketing
            self.high_poi = poi
            self.high_poi_cls = cls
        elif cls > self.high_poi_cls and poi > self.low_poi:
            self.high_poi = poi
            self.high_poi_cls = cls
        elif cls == self.high_poi_cls:
            self.high_poi = min(poi, self.high_poi)

    def suggest_poi(self):
        if self.low_poi_cls < 0.05:
            return self.low_poi / 2
        elif self.high_poi_cls > 0.05:
            return self.high_poi * 2 if self.high_poi < 1000. else -1.
        else:
            return (self.low_poi + self.high_poi) / 2


class LimitFinder(threading.Thread):
    def __init__(self, model_file, job_queue, result_queue, combined_dir):
        super(LimitFinder, self).__init__()
        self.model_file = model_file
        self.result_queue = result_queue
        self.job_queue = job_queue
        self.stop_flag = False
        self.combined_dir = combined_dir
        self.bracket_values = [BracketValues() for i in xrange(4)]

    def run(self):
        # check to see if we've already processed some points
        m1, m2 = map(int, re.search("_(\d+)_(\d+)_", self.model_file).groups())
        fn_expr = "{0}/{1}_{2}_*.root".format(self.combined_dir, m1, m2)
        files = glob(fn_expr)
        if len(files) >= 2:
            # pick up where we left off
            logger.info("Restoring progress from files: {0}".format(files))
            for f in files:
                rfile = ROOT.TFile(f)
                try:
                    poi = rfile.Get("result_sig_strength").GetLastXValue()
                except AttributeError:
                    continue
                rfile.Close()
                all_cls = cls_finder.get_cls(f)
                for j in xrange(4):
                    self.bracket_values[j].new_val(poi, all_cls[j])

        else:
            # initial values
            low_poi, high_poi = self.get_initial(self.model_file)
            logger.info("Initial poi: {0} {1}".format(low_poi, high_poi))
            self.get_cls(low_poi, 0)
            self.get_cls(high_poi, 0)

        last_poi = 0.
        for i in [2, 3, 1, 0]:
            end = False
            while not end:
                end = self.bracket_values[i].high_poi_cls - self.bracket_values[i].low_poi_cls > 0.001 or\
                      self.bracket_values[i].high_poi - self.bracket_values[i].low_poi < 0.01
                new_poi = self.bracket_values[i].suggest_poi()
                if new_poi == last_poi:
                    logger.warn("Stuck at poi {0}".format(new_poi))
                    break
                elif new_poi < 0:
                    logger.warn("Couldn't suggest poi. Probably suggested poi would be > 2000. Aborting.")
                    self.result_queue.put("FAIL")
                    return
                self.get_cls(new_poi, i)
                last_poi = new_poi

        self.result_queue.put((self.model_file,
                               self.get_obs_limit(),
                               self.get_exp_limit(),
                               self.get_expm1_limit(),
                               self.get_expp1_limit()))


    def get_cls(self, poi, i):
        logger.info("Getting CLS for {0} at {1}".format(self.model_file, poi))
        all_cls = cls_finder.run_cls_jobs(self.job_queue, self.model_file, poi, combined_path=self.combined_dir)
        for j in xrange(4):
            self.bracket_values[j].new_val(poi, all_cls[j])
        return all_cls[i]

    def get_initial(self, model_file):

        logger.info("Getting initial parameters for "+model_file)
        rfile = ROOT.TFile(model_file)

        ws = rfile.Get("combined")

        data = ws.data("obsData")

        sbmodel = ws.obj("ModelConfig")

        constr = sbmodel.GetNuisanceParameters()
        ROOT.RooStats.RemoveConstantParameters(constr)

        # run the initial fit
        ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
        ROOT.RooMsgService.instance().setSilentMode(True)
        sbmodel.GetPdf().fitTo(data, ROOT.RooFit.Constrain(constr), ROOT.RooFit.PrintEvalErrors(-1))

        poi_hat = sbmodel.GetParametersOfInterest().first().getVal()
        poi_hat_err = sbmodel.GetParametersOfInterest().first().getError()

        low = max(poi_hat - poi_hat_err, 0.001)
        high = poi_hat + poi_hat_err

        while low > high:
            high *= 2

        return low, high

  #   def get_initial_from_existing(self, files):
  #       logger.info("Getting initial from existing files"+str(files))
  #       low = 0
  #       high = 1000.
  #       for f in files:
  #           rfile = ROOT.TFile(f)
  #           x = rfile.Get("result_sig_strength").GetLastXValue()
     #    cls = rfile.Get("result_sig_strength").GetLastYValue() 
     #    if cls > 0.05:
        # low = max(low, x)
     #    else:
        # high = min(high, x)
  #           rfile.Close()
  #       assert(low < high)
  #       return low, high

    def interpolate(self, bk):
        return (bk.low_poi * (bk.low_poi_cls - 0.05) + bk.high_poi * (0.05 - bk.high_poi_cls)) / (bk.low_poi_cls - bk.high_poi_cls)

    def get_obs_limit(self):
        return self.interpolate( self.bracket_values[0] )

    def get_exp_limit(self):
        return self.interpolate( self.bracket_values[1] )

    def get_expm1_limit(self):
        return self.interpolate( self.bracket_values[2] )

    def get_expp1_limit(self):
        return self.interpolate( self.bracket_values[3] )

if __name__ == '__main__':
    import chunk
    chunker = chunk.ChunkSubmitter(100)
    chunker.start()
    results = Queue.Queue()
    lfinder = LimitFinder("../limits/slep_templates_200_0_combined_meas_model.root",
                          chunker.job_queue,
                          results)
    lfinder.start()
    try:
        while True:
            try:
                res = results.get(True, 10)
            except Queue.Empty:
                continue
            break
        print res
    finally:
        chunker.stop()
