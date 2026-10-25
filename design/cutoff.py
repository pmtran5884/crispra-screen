"""Fine read-depth sweep where the main sweep showed ONT and NGS diverge."""
import itertools, os, sys
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sweep import LIB, GUIDES, POOL, job
import sweep as W

W.CELLS=[1000,3000,10000,30000]
W.READS=[3_000,10_000,30_000,50_000,100_000,300_000,1_000_000]
W.ENRICH=[30,100,300]
W.REPS=300

if __name__=="__main__":
    grid=[(i,ns,nr,E) for i,(ns,nr,E) in enumerate(itertools.product(W.CELLS,W.READS,W.ENRICH))]
    print(f"{len(grid)} conditions x {W.REPS} reps = {len(grid)*W.REPS:,} screens",flush=True)
    with Pool(10) as p: res=p.map(W.job,grid,chunksize=1)
    pd.DataFrame(res).to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),"cutoff_sweep.csv"),index=False)
    print("done",flush=True)
