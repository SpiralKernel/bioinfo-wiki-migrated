#!/usr/bin/env python3
"""Standard-library arithmetic companions to method articles, not production callers.

Without --output this prints JSON only. An output directory must be new or empty.
All examples are artificial; no REML optimizer, FarmCPU, ADMIXTURE, MCMC,
demographic inference, SV caller or Truvari implementation is provided here.
"""
import argparse
import json
import math
from pathlib import Path


def solve(a, b):
    """Pivoted elimination for tiny nonsingular teaching systems only."""
    n = len(b)
    if len(a) != n or any(len(row) != n for row in a):
        raise ValueError('Expected a square system')
    rows = [[float(x) for x in row] + [float(y)] for row, y in zip(a, b)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(rows[i][k]))
        rows[k], rows[pivot] = rows[pivot], rows[k]
        if abs(rows[k][k]) < 1e-12:
            raise ValueError('Singular teaching system')
        scale = rows[k][k]
        rows[k] = [x / scale for x in rows[k]]
        for i in range(n):
            if i != k:
                scale = rows[i][k]
                rows[i] = [x - scale * y for x, y in zip(rows[i], rows[k])]
    return [row[-1] for row in rows]


def ols(columns, y):
    return solve([[sum(x*z for x, z in zip(c, d)) for d in columns] for c in columns],
                 [sum(x*z for x, z in zip(c, y)) for c in columns])


def fst():
    p, q, n = .25, .75, 8  # n counts haploid copies, not diploid individuals
    pooled = (p + q) / 2
    ht = 2 * pooled * (1 - pooled)
    hs = p*(1-p) + q*(1-q)
    numerator = (p-q)**2 - p*(1-p)/(n-1) - q*(1-q)/(n-1)
    denominator = p*(1-q) + q*(1-p)
    answer = dict(descriptive=(ht-hs)/ht, hudson_one_site=numerator/denominator,
                  ratio_of_sums=.21/.9, mean_of_ratios=(.2/.8+.01/.1)/2)
    assert math.isclose(answer['descriptive'], .25)
    assert math.isclose(answer['hudson_one_site'], 11/35)
    return answer


def background():
    n = 100
    bic = [n*math.log(rss/n)+k*math.log(n) for rss, k in [(100, 1), (80, 2)]]
    assert bic[1] < bic[0]
    return dict(n=n, rss=[100, 80], regression_parameters=[1, 2], bic=bic,
                note='BIC arithmetic only; not FarmCPU or BLINK')


def gblup():
    z = [[-1, -.5], [1, .5]]
    s = .875
    k = [[sum(x*y for x, y in zip(a, b))/s for b in z] for a in z]
    y = [-1, 1]
    alpha = solve([[k[i][j] + (i == j) for j in range(2)] for i in range(2)], y)
    individual = [sum(a*b for a, b in zip(row, alpha)) for row in k]
    zt = list(zip(*z))
    effects = solve([[sum(x*y for x, y in zip(a, b)) + s*(i == j)
                      for j, b in enumerate(zt)] for i, a in enumerate(zt)],
                    [sum(a*b for a, b in zip(col, y)) for col in zt])
    marker = [sum(a*b for a, b in zip(row, effects)) for row in z]
    assert max(abs(a-b) for a, b in zip(individual, marker)) < 1e-12
    assert math.isclose(individual[1], 20/27)
    return dict(scale=s, lambda_G=1, lambda_marker=s, individual=individual,
                marker=marker, effects=effects, fixed_intercept=0)


def reml():
    y = [1, 2, 3, 4]
    ss = sum((v-sum(y)/len(y))**2 for v in y)
    assert ss == 5
    return dict(n=4, fixed_rank=1, residual_SS=ss, ML_variance=ss/4,
                REML_variance=ss/3, note='Intercept-only Gaussian model; no two-component REML fit')


def ld():
    pairs = [(0, 0)]*3 + [(0, 1)] + [(1, 0)] + [(1, 1)]*3
    p = sum(x for x, _ in pairs)/len(pairs)
    q = sum(y for _, y in pairs)/len(pairs)
    d = sum(x*y for x, y in pairs)/len(pairs)-p*q
    r2 = d*d/(p*(1-p)*q*(1-q))
    assert math.isclose(r2, .25)
    return dict(p=p, q=q, D=d, r2=r2, monomorphic_r2=None,
                one_cM_in_Morgan=.01)


def roh():
    intervals = [(1, 4), (3, 5), (8, 9)]  # half-open, lengths in Mb
    union = []
    for start, end in sorted(intervals):
        if union and start <= union[-1][1]:
            union[-1] = (union[-1][0], max(end, union[-1][1]))
        else:
            union.append((start, end))
    length = sum(b-a for a, b in union)
    assert length == 5
    return dict(intervals_Mb=intervals, union_Mb=union, FROH=length/10,
                naive_double_counted=sum(b-a for a, b in intervals)/10)


def ancestry():
    q, f1, f2 = .25, .8, .2
    p = q*f1+(1-q)*f2
    assert math.isclose(p, .35)
    return dict(ancestry_weight=q, ancestral_frequencies=[f1, f2], mixed_p=p,
                expected_diploid_dosage=2*p, P_heterozygote=2*p*(1-p),
                swapped_labels_same_p=(1-q)*f2+q*f1)


def pedigree():
    a = [[1, 0, .5], [0, 1, .5], [.5, .5, 1]]
    inv_column = solve(a, [0, 0, 1])
    assert math.isclose(inv_column[-1], 2)
    return dict(A=a, inverse_of_A22=1, lower_right_of_A_inverse=inv_column[-1],
                G22=1.2, H_inverse_correction=1/1.2-1,
                F_offspring_of_uninbred_full_sibs=.25,
                A_diagonal_of_that_offspring=1.25)


def lmm():
    icc = 4/(4+1)
    de = 1+(10-1)*icc
    return dict(between_variance=4, residual_variance=1, ICC=icc,
                records=100, cluster_size=10, design_effect=de,
                variance_equivalent_n=100/de,
                note='Exchangeable balanced clusters; not a denominator-df calculation')


def history():
    theta, mu, t, g = .0004, 1e-8, .5, 25
    n0 = theta/(4*mu)
    generations = 2*n0*t
    assert math.isclose(n0, 10000)
    return dict(theta_per_site=theta, mu_per_site_per_generation=mu, N0=n0,
                scaled_time=t, generations=generations, years=generations*g,
                years_if_generation_length_doubled=generations*2*g,
                Ne_if_mu_doubled=n0/2)


def structure():
    g = [0, 0, 1, 1, 1, 1, 2, 2]
    group = [0]*4+[1]*4
    y = [0]*4+[4]*4
    naive = ols([[1]*8, g], y)
    adjusted = ols([[1]*8, g, group], y)
    assert math.isclose(naive[1], 2) and abs(adjusted[1]) < 1e-12
    return dict(naive_intercept_slope=naive, adjusted_intercept_slope_group=adjusted,
                note='Deterministic confounding arithmetic, no P-value')


def tajima(n, counts):
    s = len(counts)
    a1 = sum(1/i for i in range(1, n))
    a2 = sum(1/i**2 for i in range(1, n))
    pi = sum(2*k*(n-k)/(n*(n-1)) for k in counts)
    b1 = (n+1)/(3*(n-1))
    b2 = 2*(n*n+n+3)/(9*n*(n-1))
    e1 = (b1-1/a1)/a1
    e2 = (b2-(n+2)/(a1*n)+a2/a1**2)/(a1*a1+a2)
    variance = e1*s+e2*s*(s-1)
    return dict(pi_total=pi, theta_W_total=s/a1,
                D=(pi-s/a1)/math.sqrt(variance) if variance > 0 else None)


def selection():
    intermediate = tajima(8, [5, 5, 4, 3, 4])
    rare = tajima(8, [1]*5)
    assert math.isclose(intermediate['D'], 1.932104504549585)
    assert rare['D'] < 0
    return dict(n_chromosomes=8, intermediate_counts=intermediate, singletons=rare,
                note='D values only, no calibrated null distribution or selection test')


def bayes():
    return dict(observation=2, prior_variance=1, noise_variance=1,
                posterior_mean=1, posterior_variance=.5,
                gaussian_kernel_at_distance_one={str(h):math.exp(-1/(2*h*h)) for h in [.1, 1, 10]})


def genotype():
    eps = .01
    logs = [10*math.log(.99)+10*math.log(.01), 20*math.log(.5),
            10*math.log(.01)+10*math.log(.99)]
    top = max(logs)
    pl = [-10*(x-top)/math.log(10) for x in logs]
    assert pl[1] == 0 and pl[0] > 140
    return dict(ref_reads=10, alt_reads=10, error_probability=eps,
                log_likelihoods=logs, relative_PL=pl, best_GT='0/1',
                note='Two-allele independent-read error model; not caller GQ')


def sv_evidence():
    return dict(ref_support=14, alt_support=16, observed_support_fraction=16/30,
                ideal_CN_from_depth=2*15/30, log2_depth_ratio=math.log2(15/30),
                note='Ideal pure diploid depth example, not a CNV/SV caller')


def matching():
    # A tiny maximum cardinality matching, not Truvari's scoring/selection engine.
    bases = [1000, 1100]
    comps = [1050]
    used = {}
    def visit(b, seen):
        for c in range(len(comps)):
            if abs(bases[b]-comps[c]) <= 60 and c not in seen:
                seen.add(c)
                if c not in used or visit(used[c], seen):
                    used[c] = b
                    return True
        return False
    tp = sum(visit(b, set()) for b in range(len(bases)))
    assert tp == 1
    return dict(base_positions=bases, comp_positions=comps, tolerance=60,
                TP=tp, FP=len(comps)-tp, FN=len(bases)-tp, precision=1, recall=.5,
                naive_independent_base_hits=2)


def merge():
    pos = [1000, 1060, 1120]
    return dict(positions=pos, threshold=70, adjacent_distances=[60, 60],
                extreme_distance=pos[-1]-pos[0],
                single_link_cluster_count=1,
                note='Transitive proximity can merge endpoints farther than the threshold',
                absent_record_gt=None, observed_force_gt='0/0', force_GQ=0)


def benchmark():
    base, comp = {'a', 'b', 'c'}, {'a', 'b', 'd'}
    tp = len(base & comp)
    return dict(TP=tp, FP=len(comp-base), FN=len(base-comp),
                precision=tp/len(comp), recall=tp/len(base),
                intersection_only_agreement=1, no_truth_recall=None,
                note='Predefined artificial event IDs, not VCF/haplotype matching')


def platform():
    coverage, width = 30, 1200  # 1kb event plus two 100bp flanks
    expected = {str(r):coverage*max(0, 1-width/r) for r in [150, 16000]}
    assert expected['150'] == 0 and math.isclose(expected['16000'], 27.75)
    return dict(coverage=coverage, event_plus_flanks_bp=width,
                expected_spanning_reads=expected,
                note='Uniform fixed-length read-start model; no sequencing/mapping error')


def assembly():
    ref, query = 'AAAACCCC', 'AAAAGGCCCC'
    return dict(reference=ref, query=query, forward_length_difference=len(query)-len(ref),
                reverse_length_difference=len(ref)-len(query),
                anchor_POS_1based=4, REF='A', ALT='AGG',
                note='Known simple insertion; no genome alignment or synteny inference')


def leakage():
    # Each subject has two identical artificial records. Row holdout retains
    # another record of every test subject; subject holdout cannot use it.
    train_seen = {'A':0, 'B':10, 'C':20, 'D':30}
    row_mae = sum(abs(v-train_seen[k]) for k, v in train_seen.items())/4
    subject_train = {'A':0, 'B':10, 'C':20}
    unseen_prediction = sum(subject_train.values())/len(subject_train)
    assert row_mae == 0 and unseen_prediction == 10
    return dict(row_holdout_MAE=row_mae, new_subject_D_prediction=unseen_prediction,
                new_subject_D_truth=30, new_subject_D_absolute_error=20,
                note='Repeated-subject lookup toy, not a general CV performance claim')


def regularization():
    xtx, xty = 2, 4
    ridge = {str(l):xty/(xtx+l) for l in [0, 2, 20]}
    lasso = max(xty-1, 0)/xtx  # loss .5*SSE + lambda*abs(a), lambda=1
    assert ridge['2'] == 1 and lasso == 1.5
    return dict(ridge_loss='SSE + lambda*a^2', ridge_coefficients=ridge,
                lasso_loss='0.5*SSE + lambda*abs(a)', lasso_lambda=1, lasso_coefficient=lasso)


CASES = {'fst':fst, 'background':background, 'gblup':gblup, 'reml':reml,
         'ld':ld, 'roh':roh, 'ancestry':ancestry, 'pedigree':pedigree,
         'lmm':lmm, 'history':history, 'structure':structure, 'selection':selection,
         'bayes':bayes, 'genotype':genotype, 'sv-evidence':sv_evidence,
         'matching':matching, 'merge':merge, 'benchmark':benchmark,
         'platform':platform, 'assembly':assembly, 'leakage':leakage,
         'regularization':regularization}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--list', action='store_true', help='Print available case names and exit')
    parser.add_argument('--case', choices=['all']+list(CASES), default='all')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.list:
        print('\n'.join(CASES))
        return
    selected = CASES if args.case == 'all' else {args.case:CASES[args.case]}
    result = {'scope':'Artificial teaching arithmetic, not production statistical software',
              'cases':{name:func() for name, func in selected.items()}}
    text = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n'
    if args.output:
        if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
            parser.error('Use a new or empty output directory; existing results are preserved')
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output/'summary.json').write_text(text, encoding='utf-8')
    print(text, end='')


if __name__ == '__main__':
    main()
