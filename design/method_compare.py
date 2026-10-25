"""MAGeCK (RRA) vs the production rank hypergeometric, on identical screens.

Both methods see exactly the same simulated count table: the shallow sorted
sample as 'treat' and the 12.4M-read Mock as 'ctrl'. MAGeCK is run with its
conventional defaults (median normalisation). Note there is one sample per
condition, which is how these screens are actually run -- so MAGeCK has no
replicates from which to estimate negative-binomial dispersion, while a
rank statistic needs no variance model at all. That asymmetry is part of
what the comparison measures.
"""
import itertools, os, shutil, subprocess, sys, tempfile
import numpy as np, pandas as pd
from multiprocessing import Pool
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import simulate as S
from sweep import LIB, GUIDES, POOL

GENE=np.array([LIB.gene_names[g] for g in LIB.gene_id])
SG=np.array([f"sg{i:05d}" for i in range(len(LIB.share))])

def mageck_once(counts, hit_name, tmp):
    k=os.path.join(tmp,"c.txt")
    pd.DataFrame({"sgRNA":SG,"Gene":GENE,"ctrl":LIB.mock_counts,"treat":counts}).to_csv(
        k,sep="\t",index=False)
    r=subprocess.run(["mageck","test","-k",k,"-t","treat","-c","ctrl","-n",
                      os.path.join(tmp,"m"),"--norm-method","median"],
                     capture_output=True,text=True,cwd=tmp)
    f=os.path.join(tmp,"m.gene_summary.txt")
    if not os.path.exists(f): return None
    g=pd.read_csv(f,sep="\t").set_index("id")
    if hit_name not in g.index: return None
    row=g.loc[hit_name]
    if isinstance(row,pd.DataFrame): row=row.iloc[0]
    real=g[~g.index.str.startswith("CTRL_pseudo_")]
    return dict(fdr=float(row["pos|fdr"]), rank=int((real["pos|score"]<row["pos|score"]).sum())+1,
                lfc=float(row["pos|lfc"]),
                n_sig=int((real["pos|fdr"]<0.1).sum()),
                n_ctrl_sig=int((g[g.index.str.startswith("CTRL_pseudo_")]["pos|fdr"]<0.1).sum()))

def job(args):
    ci,ns,nr,E,reps=args
    rng=np.random.default_rng(500000+ci); out=[]
    tmp=tempfile.mkdtemp(prefix="mgk_")
    try:
        for _ in range(reps):
            b=int(POOL[rng.integers(len(POOL))]); hit=LIB.block_names[b]
            c=S.simulate_counts(LIB,GUIDES[b],E,ns,nr,rng)
            sc,pv,fdr,di=S.fast_scores(LIB,c)
            real=LIB.keep&(~LIB.block_is_control); ctrl=LIB.keep&LIB.block_is_control
            hy=dict(fdr=float(fdr[b]),rank=int((sc[real]>sc[b]).sum())+1,
                    det=bool(real[b] and di[b]==0 and fdr[b]<0.1),
                    n_sig=int(((fdr<0.1)&real&(di==0)).sum()),
                    n_ctrl_sig=int(((fdr<0.1)&ctrl&(di==0)).sum()))
            mg=mageck_once(c,hit,tmp)
            out.append(dict(hyper_det=hy["det"],hyper_rank=hy["rank"],hyper_fdr=hy["fdr"],
                            hyper_nsig=hy["n_sig"],hyper_ctrlsig=hy["n_ctrl_sig"],
                            mg_det=bool(mg and mg["fdr"]<0.1),
                            mg_rank=mg["rank"] if mg else np.nan,
                            mg_fdr=mg["fdr"] if mg else np.nan,
                            mg_nsig=mg["n_sig"] if mg else np.nan,
                            mg_ctrlsig=mg["n_ctrl_sig"] if mg else np.nan))
    finally:
        shutil.rmtree(tmp,ignore_errors=True)
    d=pd.DataFrame(out)
    return dict(n_sorted=ns,n_reads=nr,E=E,reps=reps,
                hyper_detected=d.hyper_det.mean(),mageck_detected=d.mg_det.mean(),
                hyper_rank1=(d.hyper_rank==1).mean(),mageck_rank1=(d.mg_rank==1).mean(),
                hyper_rank_med=d.hyper_rank.median(),mageck_rank_med=d.mg_rank.median(),
                hyper_fdr_med=d.hyper_fdr.median(),mageck_fdr_med=d.mg_fdr.median(),
                hyper_nsig_med=d.hyper_nsig.median(),mageck_nsig_med=d.mg_nsig.median(),
                hyper_ctrl_fp=(d.hyper_ctrlsig>0).mean(),mageck_ctrl_fp=(d.mg_ctrlsig>0).mean())

if __name__=="__main__":
    grid=[(i,ns,nr,E,50) for i,(ns,nr,E) in enumerate(
        itertools.product([300,1000,3000],[30_000],[100,300,1000]))]
    print(f"{len(grid)} conditions x 50 reps = {len(grid)*50} paired screens",flush=True)
    with Pool(9) as p: res=p.map(job,grid,chunksize=1)
    pd.DataFrame(res).to_csv(os.path.join(HERE,"method_comparison.csv"),index=False)
    print("done",flush=True)
