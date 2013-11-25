import threading
import cls_finder
import Queue
import ROOT
import logging
logging.basicConfig(level=logging.DEBUG)


class BracketValues(object):
    def __init__(self):
        self.low_poi = 0.
        self.low_poi_cls = 1.
        self.high_poi = 1000.
        self.high_poi_cls = 0.

    def new_val(self, poi, cls):
        if cls > 0.05:
            self.new_low(poi, cls)
        else:
            self.new_high(poi, cls)
        logging.debug("Low: {0}, High: {1}".format((self.low_poi, self.low_poi_cls), (self.high_poi, self.high_poi_cls)))

    def new_low(self, poi, cls):
        if self.low_poi_cls < 0.05:
            # we're still bracketing
            self.low_poi = poi
            self.low_poi_cls = cls
        elif cls < self.low_poi_cls:
            self.low_poi = poi
            self.low_poi_cls = cls
        elif cls == self.low_poi_cls:
            self.low_poi = max(poi, self.low_poi)

    def new_high(self, poi, cls):
        if self.high_poi_cls > 0.05:
            # we're still bracketing
            self.high_poi = poi
            self.high_poi_cls = cls
        elif cls > self.high_poi_cls:
            self.high_poi = poi
            self.high_poi_cls = cls
        elif cls == self.high_poi_cls:
            self.high_poi = min(poi, self.high_poi)

    def suggest_poi(self):
        if self.low_poi_cls < 0.05:
            return self.low_poi / 2
        elif self.high_poi_cls > 0.05:
            return self.high_poi * 2
        else:
            return (self.low_poi + self.high_poi) / 2


class LimitFinder(threading.Thread):
    def __init__(self, model_file, job_queue, result_queue):
        super(LimitFinder, self).__init__()
        self.model_file = model_file
        self.result_queue = result_queue
        self.job_queue = job_queue
        self.stop_flag = False
        self.bracket_values = [BracketValues() for i in xrange(4)]

    def run(self):
        # initial values
        low_poi, high_poi = self.get_initial(self.model_file)
        self.get_cls(low_poi, 0)
        self.get_cls(high_poi, 0)

        last_poi = 0.
        for i in [2, 3, 1, 0]:
            end = self.bracket_values[i].high_poi_cls - self.bracket_values[i].low_poi_cls > 0.001 or\
                  self.bracket_values[i].high_poi - self.bracket_values[i].low_poi < 0.01
            while not end:
                new_poi = self.bracket_values[i].suggest_poi()
                if new_poi == last_poi:
                    logging.warn("Stuck at poi {0}".format(new_poi))
                    break
                self.get_cls(new_poi, i)
                last_poi = new_poi

        self.result_queue.put((self.model_file,
                               self.get_obs_limit(),
                               self.get_exp_limit(),
                               self.get_expm1_limit(),
                               self.get_expp1_limit()))

    def get_cls(self, poi, i):
        logging.debug("Getting CLS for {0} at {1}".format(self.model_file, poi))
        all_cls = cls_finder.run_cls_jobs(self.job_queue, self.model_file, poi)
        for j in xrange(4):
            self.bracket_values[j].new_val(poi, all_cls[j])
        return all_cls[i]

    def get_initial(self, model_file):
        rfile = ROOT.TFile(model_file)

        ws = rfile.Get("combined")

        data = ws.data("obsData")

        sbmodel = ws.obj("ModelConfig")

        constr = sbmodel.GetNuisanceParameters()
        ROOT.RooStats.RemoveConstantParameters(constr)

        # run the initial fit
        sbmodel.GetPdf().fitTo(data, ROOT.RooFit.Constrain(constr))

        poi_hat = sbmodel.GetParametersOfInterest().first().getVal()
        poi_hat_err = sbmodel.GetParametersOfInterest().first().getError()

        return max(poi_hat - 2 * poi_hat_err, 0.001), poi_hat + 2 * poi_hat_err

    def get_obs_limit(self):
        return self.bracket_values[0].suggest_poi()

    def get_exp_limit(self):
        return self.bracket_values[1].suggest_poi()

    def get_expm1_limit(self):
        return self.bracket_values[2].suggest_poi()

    def get_expp1_limit(self):
        return self.bracket_values[3].suggest_poi()

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
