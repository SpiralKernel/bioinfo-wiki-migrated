#!/usr/bin/env python3
"""Standard-library arithmetic lessons; no platform or professional-method validation.

Run from the project root. Each case writes JSON to a NEW/EMPTY output directory.
No external input, network, package installation or biological data are required.
"""
import argparse
import json
import math
from pathlib import Path
import platform
import statistics as st

PROJECT = Path(__file__).resolve().parents[1]


def close(a, b):
    assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-9), (a, b)


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def rca():
    # RCA copies are repeats of ONE circle, not independent original molecules.
    originals, copies = 4, 100
    r = dict(original_circles=originals, repeats_per_DNB=copies,
             total_template_repeats=originals*copies,
             original_allele_counts=dict(A=3, G=1), allele_fraction_G=1/4,
             after_library_bias=dict(A=6, G=1), biased_fraction_G=1/7)
    close(r['allele_fraction_G'], .25)
    return r


def illumina():
    f, cycles = .01, 100
    r = dict(phasing_fraction=f, cycles=cycles,
             synchronous_fraction=(1-f)**cycles,
             error_Q20=10**(-20/10), error_Q30=10**(-30/10),
             expected_errors_150bp_Q30=150*10**(-30/10))
    close(r['expected_errors_150bp_Q30'], .15)
    return r


def nanopore():
    # Artificial overlapping 2-mer means. This is not a nanopore basecaller.
    model={'AA': 10., 'AC': 15., 'CG': 20.}
    observed=[10.2, 14.8, 20.1]
    paths=['AACG','AAAC']
    score={s:sum((observed[i]-model[s[i:i+2]])**2 for i in range(3)) for s in paths}
    assert min(score, key=score.get)=='AACG'
    return dict(observed=observed, candidate_squared_error=score, chosen='AACG',
                note='Toy 2-mer regression; no real current or trained ONT model')


def pcr():
    r=dict(start=10, cycles=10, ideal=10*2**10,
           efficiency_0_8=10*1.8**10,
           relative_bias_after_20_cycles=(1.9/1.8)**20)
    close(r['ideal'], 10240)
    return r


def hifi():
    p=.1
    majority=sum(math.comb(5,k)*p**k*(1-p)**(5-k) for k in range(3,6))
    close(majority, .00856)
    return dict(independent_binary_error=p, passes=5,
                majority_error=majority, original_molecules=1,
                note='Binary independent-error illustration, not CCS algorithm or HiFi accuracy')


def sanger():
    ref='ACGTTACG'
    deletion='ACGTACG'
    mixture=[{'position':i+1,'ref':a,'other':b,'mixed':a!=b}
             for i,(a,b) in enumerate(zip(ref,deletion))]
    return dict(reference=ref, deletion_allele=deletion, aligned_by_read_index=mixture,
                reverse_complement=ref.translate(str.maketrans('ACGT','TGCA'))[::-1],
                note='Base-index mixture illustration, no AB1 or quantitative peak model')


def ion():
    signal=2.6
    likelihood={str(n):math.exp(-.5*((signal-n)/.4)**2) for n in [2,3]}
    total=sum(likelihood.values())
    posterior={n:v/total for n,v in likelihood.items()}
    assert posterior['3']>posterior['2']
    return dict(artificial_flow_signal=signal, sigma=.4,
                normalized_two_candidate_likelihood=posterior,
                note='Linear Gaussian two-candidate toy; not instrument calibration')


def coverage():
    n,l,g=1000,100,100000
    return dict(reads=n, read_length=l, genome_length=g,
                mean_depth=n*l/g, poisson_zero_probability=math.exp(-n*l/g),
                paired_150_fragment_250_unique_span=250,
                paired_150_sequenced_bases=300)


def interaction():
    # Four independent Poisson cells: an illustration of log-linear interaction.
    y=[[10,20],[20,80]]
    rows=[sum(r) for r in y];cols=[sum(y[i][j] for i in range(2)) for j in range(2)]
    total=sum(rows);expected=[[rows[i]*cols[j]/total for j in range(2)] for i in range(2)]
    deviance=2*sum(y[i][j]*math.log(y[i][j]/expected[i][j]) for i in range(2) for j in range(2))
    effect=math.log2(80/20)-math.log2(20/10)
    close(effect,1.)
    return dict(cell_counts=y, ratio_of_ratios=2., interaction_log2=effect,
                reduced_expected=expected, poisson_deviance=deviance,
                note='Poisson illustration without dispersion/replicates; no RNA P value')


def rna():
    # Known library factors, not DESeq2 size-factor estimation or apeglm.
    counts=[20,40,80,160];size=[1,2,1,2]
    normalized=[c/s for c,s in zip(counts,size)]
    beta,se,tau=3.,2.,1.
    close(math.log2(st.mean(normalized[2:])/st.mean(normalized[:2])),2.)
    return dict(counts=counts, known_size_factors=size, normalized=normalized,
                normalized_log2_fold_change=2., NB_mean=100, NB_alpha=.1,
                NB_variance=100+.1*100**2,
                normal_prior_toy=dict(beta=beta,se=se,tau=tau,posterior_mean=beta*tau*tau/(tau*tau+se*se)))


def metabolomics():
    features=[[1000,2000,3000],[1,2,3]]
    raw_variance=[st.variance(x) for x in features]
    standardized=[[(v-st.mean(x))/st.stdev(x) for v in x] for x in features]
    return dict(feature_sample_values=features, sample_variances=raw_variance,
                standardized_sample_variances=[st.variance(x) for x in standardized],
                fold_change_raw=200/100, log2_fold_change=math.log2(200/100),
                note='Scaling arithmetic only; no metabolite identification/PCA/PLS fit')


def pseudobulk():
    # Donor means and cell counts deliberately unequal; not a count-model fit.
    donors=[('A1','A',100,1.),('A2','A',10,9.),('B1','B',100,4.),('B2','B',10,6.)]
    r={}
    for group in ['A','B']:
        rows=[d for d in donors if d[1]==group]
        r[group]=dict(donors=len(rows),cells=sum(d[2] for d in rows),
                      donor_equal_mean=st.mean(d[3] for d in rows),
                      cell_weighted_mean=sum(d[2]*d[3] for d in rows)/sum(d[2] for d in rows))
    close(r['A']['donor_equal_mean'],5.);close(r['B']['donor_equal_mean'],5.)
    return dict(groups=r,note='Unequal weighting example, not pseudobulk normalization or DE test')


def neighbors():
    x=[0., .1, 3.,3.1]
    def graph(k):
        edges=set()
        for i in range(len(x)):
            js=sorted((j for j in range(len(x)) if i!=j),key=lambda j:(abs(x[i]-x[j]),j))[:k]
            edges.update(tuple(sorted((i,j))) for j in js)
        return sorted([list(e) for e in edges])
    assert len(graph(1))==2 and len(graph(2))>2
    return dict(coordinates=x, k1_edges=graph(1),k2_edges=graph(2),
                note='Symmetrized unweighted kNN toy, not Scanpy fuzzy graph/Leiden/UMAP')


def multiomics():
    z=[-1.,0.,1.];w=[2.,-3.]
    y=[[a*b for b in w] for a in z]
    flipped=[[(-a)*(-b) for b in w] for a in z]
    assert y==flipped
    return dict(factor=z, loadings=w, reconstructed=y,
                sign_flip_same_reconstruction=True,
                block_trace_before_weighting=dict(RNA_100_features=100,protein_2_features=2),
                note='Exact rank-one illustration, no MOFA/DIABLO fit')


def enrichment():
    N,K,n,k=20,5,4,3
    p=sum(math.comb(K,i)*math.comb(N-K,n-i) for i in range(k,min(K,n)+1))/math.comb(N,n)
    hits={1,2,6};running=[];s=0.
    for i in range(1,9):
        s+=1/len(hits) if i in hits else -1/(8-len(hits));running.append(s)
    close(p,155/4845);close(max(running),2/3);close(running[-1],0.)
    return dict(ORA=dict(N=N,K=K,n=n,k=k,upper_tail_p=p),
                unweighted_running_score=running, positive_ES=max(running),
                note='ORA exact arithmetic/unweighted running score; no GSEA NES/FDR')


def epigenetics():
    counts, factors = [20,40], [1,2]
    methylation=[]
    for methylated,total in [(2,4),(50,100)]:
        fraction=methylated/total
        methylation.append(dict(methylated=methylated,total=total,fraction=fraction,
                                binomial_SE=math.sqrt(fraction*(1-fraction)/total)))
    return dict(peak_counts=counts,library_factors=factors,
                normalized_counts=[c/s for c,s in zip(counts,factors)],methylation=methylation,
                note='Independent binomial SE illustration; no biological overdispersion or DMR fit')


def microbiome():
    absolute_a=[10,10];absolute_b=[20,10]
    a=[x/sum(absolute_a) for x in absolute_a];b=[x/sum(absolute_b) for x in absolute_b]
    clr=lambda x:[math.log(v)-st.mean(math.log(t) for t in x) for v in x]
    return dict(absolute_A=absolute_a,absolute_B=absolute_b,relative_A=a,relative_B=b,
                taxon2_relative_fold_change=b[1]/a[1],clr_A=clr(a),clr_B=clr(b),
                Bray_Curtis_relative=sum(abs(x-y) for x,y in zip(a,b))/sum(a+b),
                Shannon_A=-sum(x*math.log(x) for x in a),
                note='No absolute-abundance inference from sequencing counts alone')


def batch():
    # Both groups present in both batches; subtract known additive batch shift.
    rows=[dict(batch='b1',group='A',value=1),dict(batch='b1',group='B',value=3),
          dict(batch='b2',group='A',value=11),dict(batch='b2',group='B',value=13)]
    corrected=[r['value']-(10 if r['batch']=='b2' else 0) for r in rows]
    assert corrected==[1,3,1,3]
    confounded_rows=[('A','b1'),('A','b1'),('B','b2'),('B','b2')]
    group_column=[int(g=='B') for g,b in confounded_rows]
    batch_column=[int(b=='b2') for g,b in confounded_rows]
    return dict(crossed_rows=rows,corrected=corrected,
                group_difference_preserved=st.mean(corrected[1::2])-st.mean(corrected[::2]),
                confounded_design_columns_equal=group_column==batch_column,
                note='Known additive shift; no Harmony/scVI or learned batch effect')


def proteomics():
    # A simplified equal-target/decoy null competition convention, not q-values.
    p=[.001,.01,.04,.2];m=len(p)
    bh=[0.]*m;current=1.
    for i in reversed(range(m)):
        current=min(current,p[i]*m/(i+1));bh[i]=current
    ratios=[math.log2(b/a) for a,b in [(100,200),(200,400)]]
    return dict(target_hits=100,decoy_hits=2,illustrative_decoy_over_target=2/100,
                raw_p=p,BH=bh,peptide_log2_ratios=ratios,protein_median_log2_ratio=st.median(ratios),
                note='Illustrative D/T and BH only, not a proteomics search-engine FDR implementation')


def cnv():
    depth=[100,100,50,50,150,150]
    ratios=[math.log2(x/100) for x in depth]
    cn=[2*2**x for x in ratios]
    purity=.5;observed_ratio=(purity*1+(1-purity)*2)/2
    assert cn==[2.,2.,1.,1.,3.,3.]
    return dict(depth=depth,baseline=100,log2_ratios=ratios,ideal_CN=cn,
                purity_half_CN1_observed_log2=math.log2(observed_ratio),
                note='Known baseline, no GC correction/segmentation/absolute CN fitting')


def intervals():
    a=(10,20);b=(20,30);c=(15,25)
    overlap=lambda u,v:max(0,min(u[1],v[1])-max(u[0],v[0]))
    assert overlap(a,b)==0 and overlap(a,c)==5
    query=(40,50);genes={'left_gene':(10,30),'right_gene':(60,80)}
    gap=lambda u,v:max(0,v[0]-u[1],u[0]-v[1])
    return dict(half_open_A=a,B=b,C=c,overlap_AB=overlap(a,b),overlap_AC=overlap(a,c),
                position_1based_20_BED=[19,20],
                tie_query=query,tie_genes=genes,
                nearest_boundary_gap_tie={name:gap(query,interval) for name,interval in genes.items()},
                note='Zero-based geometric gaps, not bedtools closest -d distance convention')


def rare():
    # Orthogonal centered genotype columns and centered phenotype residuals.
    g1=[1.,-1.,0.,0.];g2=[0.,0.,1.,-1.];residual=[1.,-1.,-1.,1.]
    scores=[dot(g1,residual),dot(g2,residual)]
    burden=[a+b for a,b in zip(g1,g2)]
    close(dot(burden,residual),0.);close(sum(u*u for u in scores),8.)
    return dict(centered_columns=[g1,g2],residual=residual,variant_scores=scores,
                burden_score=dot(burden,residual),unit_weight_kernel_Q=sum(u*u for u in scores),
                note='Raw score cancellation illustration, no SKAT null distribution or P value')


CASES={f.__name__:f for f in [rca,illumina,nanopore,pcr,hifi,sanger,ion,coverage,
       interaction,rna,metabolomics,pseudobulk,neighbors,multiomics,enrichment,
       epigenetics,microbiome,batch,proteomics,cnv,intervals,rare]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=['all',*CASES])
    parser.add_argument('--list',action='store_true',help='List arithmetic cases and exit without writing output')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.list:
        print('\n'.join(CASES))
        return
    if args.case is None:
        parser.error('--case is required unless --list is used')
    output=args.output or PROJECT/'examples/output/methods-omics'/args.case
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.exit(1,'Use a new/empty output directory; previous results are preserved\n')
    output.mkdir(parents=True,exist_ok=True)
    results={name:fn() for name,fn in CASES.items() if args.case in ['all',name]}
    report=dict(python=platform.python_version(),scope='Artificial standard-library arithmetic only',cases=results)
    (output/'summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False))
    print('Output:',output/'summary.json')


if __name__=='__main__':
    main()
