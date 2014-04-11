import get_limits
import chunk
import Queue
import re
import sys
import logging
from time import sleep


def main(files, out_file, combined_dir, lsf_queue):
    chunker = chunk.ChunkSubmitter(lsf_queue, 10000)
    chunker.start()
    try:
        results = Queue.Queue()
        lfinders = []
        for f in files:
            logger.info("Starting LimitFinder for "+f)
            lfinder = get_limits.LimitFinder(f,
                                  chunker.job_queue,
                                  results,
                                  combined_dir)
            lfinder.daemon = True
            lfinder.start()
            lfinders.append(lfinder)
        with open(out_file, 'w') as f:
            i = 0
            while i < len(files):
                while True:
                    try:
                        res = results.get(True, 0.1)
                    except Queue.Empty:
                        sleep(10)
                        continue

                    i += 1
                    break
                if res == "FAIL":
                    continue
                model, obs, exp, expm1, expp1 = map(str, res)
                m1, m2 = re.search("_(\d*?)_(\d*?)_", model).groups()[:2]
                f.write("\t".join((m1, m2, obs, exp, expm1, expp1)))
                f.write("\n")
                logger.info("Finished point {0}, {1}".format(m1, m2))
                f.flush()
    finally:
        chunker.stop()

if __name__ == '__main__':
    logger = logging.getLogger('main')
    logger.addHandler(logging.StreamHandler(sys.stdout))
    logger.setLevel(logging.INFO)

    limit_logger = logging.getLogger("main.get_limit")
    limit_logger.setLevel(logging.DEBUG)
    # limit_logger.addHandler(logging.StreamHandler(sys.stdout))

    out_file = sys.argv[1]
    files = sys.argv[4:]
    combined_dir = sys.argv[2]
    lsf_queue = sys.argv[3]
    sys.argv = [sys.argv[0], '-b']
    logger.info("Go!")
    main(files, out_file, combined_dir, lsf_queue)
