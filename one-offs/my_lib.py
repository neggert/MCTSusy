from collections import defaultdict

def results_to_dict(r):
    results = defaultdict(lambda: defaultdict(dict))
    results['sf']['sum']['low'] = (r.sf.sum.low.val, r.sf.sum.low.error)
    results['sf']['sum']['high'] = (r.sf.sum.high.val, r.sf.sum.high.error)
    results['sf']['vv']['low'] = (r.sf.vv.low.val, r.sf.vv.low.error)
    results['sf']['vv']['high'] = (r.sf.vv.high.val, r.sf.vv.high.error)
    results['sf']['top']['low'] = (r.sf.top.low.val, r.sf.top.low.error)
    results['sf']['top']['high'] = (r.sf.top.high.val, r.sf.top.high.error)
    results['sf']['z']['low'] = (r.sf.z.low.val, r.sf.z.low.error)
    results['sf']['z']['high'] = (r.sf.z.high.val, r.sf.z.high.error)
    results['sf']['fake']['low'] = (r.sf.fake.low.val, r.sf.fake.low.error)
    results['sf']['fake']['high'] = (r.sf.fake.high.val, r.sf.fake.high.error)
    results['of']['sum']['low'] = (r.of.sum.low.val, r.of.sum.low.error)
    results['of']['sum']['high'] = (r.of.sum.high.val, r.of.sum.high.error)
    results['of']['vv']['low'] = (r.of.vv.low.val, r.of.vv.low.error)
    results['of']['vv']['high'] = (r.of.vv.high.val, r.of.vv.high.error)
    results['of']['top']['low'] = (r.of.top.low.val, r.of.top.low.error)
    results['of']['top']['high'] = (r.of.top.high.val, r.of.top.high.error)
    results['of']['z']['low'] = (r.of.z.low.val, r.of.z.low.error)
    results['of']['z']['high'] = (r.of.z.high.val, r.of.z.high.error)
    results['of']['fake']['low'] = (r.of.fake.low.val, r.of.fake.low.error)
    results['of']['fake']['high'] = (r.of.fake.high.val, r.of.fake.high.error)

    results['of']['generated_sum']['high'] = r.of.generated_sum.val
    results['sf']['generated_sum']['high'] = r.sf.generated_sum.val

    return results
