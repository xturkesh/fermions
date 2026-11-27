import sys, numpy as np
sys.path.insert(0, '../temporal_ipr')
import distribution as dist                      # reuses record_probs (prints its own table on import? guarded below)
d = 2; out = {}
for D in [1, 4, 16, 64]:
    vals = []
    for _ in range(30):
        pr = dist.record_probs(d, D, 12); vals.append(np.log((d**12)*pr[pr > 0]))
    out['H_D%d' % D] = np.concatenate(vals)
rows = []
for D in [4, 16, 64]:
    for n in [2, 4, 6, 8, 10, 12]:
        v = np.concatenate([np.log((d**n)*dist.record_probs(d, D, n)) for _ in range(20)])
        rows.append((D, n, v.var()))
out['V'] = np.array(rows)
np.savez('data_ipr_mc.npz', **out); print("ok", [(k, out[k].shape) for k in out])
