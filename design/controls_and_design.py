"""Two focused sweeps: control-guide count, and library shape at a fixed cell budget."""
import itertools, os, sys
import numpy as np, pandas as pd
from multiprocessing import Pool
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import simulate as S

AB=pd.read_csv(os.path.join(HERE,"mock_abundance.csv"))
POOL_AB=AB.mock_counts.values
INPUT_CELLS=12_500_000

def _prep(lib):
    guides={b:np.flatnonzero(lib.gene_id==g) for b,g in enumerate(lib.block_gene)}
    elig=np.array([b for b in range(len(lib.block_gene))
                   if lib.sizes[b]==lib.sizes.max() and not lib.block_is_control[b]
                   and (lib.mock_counts[guides[b]]>=1).all()])
    return guides,elig

# ---- control-guide count -------------------------------------------------
def job_ctrl(args):
    ci,n_ctrl,ns,nr,E,reps=args
    rng=np.random.default_rng(300000+ci)
    lib=S.build_synthetic(5891,4,n_ctrl,POOL_AB,np.random.default_rng(11))
    guides,elig=_prep(lib)
    out=[]
    for _ in range(reps):
        b=int(elig[rng.integers(len(elig))])
        c=S.simulate_counts_pooled(lib,guides[b],E,INPUT_CELLS,ns,nr,rng)
        score,pval,fdr,direction=S.fast_scores(lib,c)
        real=lib.keep&(~lib.block_is_control); ctrl=lib.keep&lib.block_is_control
        det=bool(real[b] and direction[b]==0 and fdr[b]<0.1)
        out.append(dict(detected=det,
                        rank=int((score[real]>score[b]).sum())+1 if real[b] else -1,
                        n_ctrl_sig=int(((fdr<0.1)&ctrl&(direction==0)).sum()),
                        n_ctrl_blocks=int(ctrl.sum()),
                        n_real_sig=int(((fdr<0.1)&real&(direction==0)).sum()),
                        null_sd=float(np.std(score[ctrl])) if ctrl.sum()>2 else np.nan))
    d=pd.DataFrame(out)
    return dict(n_controls=n_ctrl,n_sorted=ns,n_reads=nr,E=E,reps=reps,
                detected=d.detected.mean(),rank_med=d['rank'].replace(-1,np.nan).median(),
                ctrl_fp_rate=(d.n_ctrl_sig>0).mean(),
                ctrl_sig_per_screen=d.n_ctrl_sig.mean(),
                ctrl_blocks=d.n_ctrl_blocks.iloc[0],
                real_sig_med=d.n_real_sig.median(),null_sd=d.null_sd.mean())

# ---- library shape -------------------------------------------------------
def job_design(args):
    ci,ng,gpg,ns,nr,E,reps=args
    rng=np.random.default_rng(400000+ci)
    lib=S.build_synthetic(ng,gpg,1000,POOL_AB,np.random.default_rng(12))
    guides,elig=_prep(lib)
    out=[]
    for _ in range(reps):
        b=int(elig[rng.integers(len(elig))])
        c=S.simulate_counts_pooled(lib,guides[b],E,INPUT_CELLS,ns,nr,rng)
        score,pval,fdr,direction=S.fast_scores(lib,c)
        real=lib.keep&(~lib.block_is_control)
        out.append(dict(detected=bool(real[b] and direction[b]==0 and fdr[b]<0.1),
                        rank1=bool(real[b] and (score[real]>score[b]).sum()==0),
                        guides_seen=int((c[guides[b]]>0).sum())))
    d=pd.DataFrame(out)
    return dict(n_genes=ng,guides_per_gene=gpg,n_guides=len(lib.share),
                coverage=round(INPUT_CELLS/len(lib.share)),
                n_sorted=ns,n_reads=nr,E=E,reps=reps,
                detected=d.detected.mean(),rank1=d.rank1.mean(),
                guides_seen_med=d.guides_seen.median())

if __name__=="__main__":
    which=sys.argv[1]
    if which=="ctrl":
        grid=[(i,nc,ns,nr,E,400) for i,(nc,ns,nr,E) in enumerate(
            itertools.product([0,250,500,1000,2000],[300,1000,3000],[30000],[100,300]))]
        fn,fun="control_count.csv",job_ctrl
    else:
        grid=[(i,ng,gpg,ns,nr,E,300) for i,(ng,gpg,ns,nr,E) in enumerate(
            itertools.product([1500,3000,6000],[4,6,8],[300,1000,3000,10000],[30000],[100,300]))]
        fn,fun="library_design.csv",job_design
    print(f"{len(grid)} conditions",flush=True)
    with Pool(10) as p: res=p.map(fun,grid,chunksize=1)
    pd.DataFrame(res).to_csv(os.path.join(HERE,fn),index=False); print("done",flush=True)
