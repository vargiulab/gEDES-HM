import time
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, "common"))
from adk_geom import iter_models, coords
from mapping import keymap_for, ZEN

t0 = time.time()
out = open('model_coordinates.tsv', 'w')
out.write("target\tmodel\tdLID\tdNMP\ttheta1\ttheta2\tcov_dist_LID\tcov_ang_LID\tcov_ang_CORE\tcov_ang_NMP\tcov_ang_HINGE\n")
n = 0
for i in range(1, 23):
    t = str(i)
    km, _ = keymap_for(t)
    path = f"{ZEN}/{t}/concat_selected_clean.pdb.gz"
    cnt = 0
    for mno, res in iter_models(path):
        c = coords(res, keymap=km)
        out.write(f"{i}\t{mno}\t{c['dLID']:.4f}\t{c['dNMP']:.4f}\t{c['theta1']:.4f}\t{c['theta2']:.4f}\t"
                  f"{c['cov_dist_LID']:.4f}\t{c['cov_ang_LID']:.4f}\t{c['cov_ang_CORE']:.4f}\t"
                  f"{c['cov_ang_NMP']:.4f}\t{c['cov_ang_HINGE']:.4f}\n")
        cnt += 1; n += 1
    print(f"  target {i:>2}: {cnt} models  (LIDcov {c['cov_dist_LID']:.2f})", flush=True)
out.close()
print(f"wrote model_coordinates.tsv, {n} models, {time.time()-t0:.0f} s")
