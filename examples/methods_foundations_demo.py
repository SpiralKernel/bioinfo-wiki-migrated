#!/usr/bin/env python3
"""Worked examples accompanying the Chinese methods chapters; artificial inputs only.

Run from project root with the project Python environment. NumPy/SciPy are needed
for the corresponding cases; --list and standard-library cases do not import them.
Code blocks are kept in sync with the publication via recorded SHA256 checks.
"""
import argparse
from contextlib import redirect_stdout
from io import StringIO
import hashlib
import json
from pathlib import Path
import sys

PROJECT = Path(__file__).resolve().parents[1]
if str(PROJECT) not in sys.path:
    sys.path.insert(0, str(PROJECT))

def case_blast():
    from examples.alignment_demo import seed_extend, smith_waterman
    hit = seed_extend('ABCD', 'XXABCDYY', k=2)[0]
    print('simple_hit', hit)
    print('seed_miss', seed_extend('ABCD', 'AXCY', k=2))
    print('exact_local', smith_waterman('ABCD', 'AXCY')['score'])
    assert hit['score'] == 8
    assert seed_extend('ABCD', 'AXCY', k=2) == []
    assert smith_waterman('ABCD', 'AXCY')['score'] > 0


def case_pca():
    import numpy as np
    X = np.array([[1.,2.],[2.,4.],[3.,6.]])
    mean = X.mean(axis=0)
    C = X - mean
    U, s, Vt = np.linalg.svd(C, full_matrices=False)
    scores = C @ Vt.T
    ratio = s**2 / (s**2).sum()
    print('singular_values', np.round(s,6))
    print('explained_ratio', np.round(ratio,6))
    print('PC1_absolute_scores', np.round(abs(scores[:,0]),6))
    reconstructed = scores[:,:1] @ Vt[:1,:] + mean
    print('reconstruction_error', round(float(abs(reconstructed-X).max()),12))
    assert np.allclose(ratio,[1.,0.])
    assert np.allclose(reconstructed,X)


def case_smith_waterman():
    from examples.alignment_demo import smith_waterman
    r = smith_waterman('ABCD', 'XXABCDYY')
    print(r)
    assert r['score'] == 8
    assert r['start'] == [0, 2] and r['end'] == [4, 6]
    assert r['query'] == r['subject'] == 'ABCD'
    gapped = smith_waterman('ABCDE', 'ABXCDE')
    print('gap_example', gapped['score'], gapped['query'], gapped['subject'])
    assert gapped['score'] == 8
    assert smith_waterman('AAAA', 'BBBB')['score'] == 0


def case_t_tests():
    import numpy as np
    from scipy import stats
    A = np.array([1.,2.,3.,4.])
    B = np.array([2.,4.,5.,7.])
    effect = B.mean()-A.mean()
    se = np.sqrt(A.var(ddof=1)/len(A)+B.var(ddof=1)/len(B))
    welch = stats.ttest_ind(B,A,equal_var=False,alternative='two-sided')
    paired = stats.ttest_rel(B,A,alternative='two-sided')
    print('effect_SE_t', round(effect,6), round(se,6), round(float(welch.statistic),6))
    print('Welch_p_df', round(float(welch.pvalue),6), round(float(welch.df),6))
    print('paired_t_p', round(float(paired.statistic),6), round(float(paired.pvalue),6))
    assert np.isclose(welch.statistic,effect/se)
    assert np.isclose(paired.statistic,np.sqrt(24))


def case_mle_map():
    from math import log, log1p
    n, k, a, b = 10, 8, 2, 2
    mle = k/n
    map_estimate = (k+a-1)/(n+a+b-2)
    posterior_mean = (k+a)/(n+a+b)
    def ll(p):
        return k*log(p)+(n-k)*log1p(-p)
    print('MLE_MAP_mean', round(mle,6), round(map_estimate,6), round(posterior_mean,6))
    print('likelihood_prefers_MLE', ll(mle) > ll(map_estimate))
    assert abs(map_estimate-.75) < 1e-12
    assert ll(mle) > ll(.5)


def case_hypothesis():
    from math import comb
    n, observed = 10, 8
    prob = [comb(n,k)/2**n for k in range(n+1)]
    p_greater = sum(prob[observed:])
    p_two = sum(v for v in prob if v <= prob[observed]+1e-15)
    print('one_sided', p_greater)
    print('probability_ordered_two_sided', p_two)
    assert p_greater == 56/1024 and p_two == 112/1024


def case_classification():
    pairs = [(0.9,1),(0.8,0),(0.7,1),(0.1,0)]
    t = .75
    TP = sum(y==1 and s>=t for s,y in pairs)
    FP = sum(y==0 and s>=t for s,y in pairs)
    FN = sum(y==1 and s<t for s,y in pairs)
    TN = sum(y==0 and s<t for s,y in pairs)
    pos = [s for s,y in pairs if y==1]
    neg = [s for s,y in pairs if y==0]
    auc = sum((a>b)+.5*(a==b) for a in pos for b in neg)/(len(pos)*len(neg))
    tp, ap = 0, 0.
    for rank,(score,y) in enumerate(pairs,1):
        if y:
            tp += 1
            ap += (tp/rank)/len(pos)
    print('TP_FP_FN_TN',TP,FP,FN,TN,'AUC',auc,'AP',round(ap,6))
    assert (TP,FP,FN,TN)==(1,1,1,1) and auc==.75
    assert abs(ap-5/6)<1e-12


def case_contingency():
    import numpy as np
    from scipy.stats import chi2_contingency, fisher_exact
    T = np.array([[8,2],[1,9]])
    chi = chi2_contingency(T,correction=False)
    fisher = fisher_exact(T,alternative='two-sided')
    print('expected', chi.expected_freq.tolist())
    print('chi2_df_p', round(float(chi.statistic),6), chi.dof, round(float(chi.pvalue),6))
    print('OR_Fisher_p', float(fisher.statistic), round(float(fisher.pvalue),6))
    print('risk_ratio', (8/10)/(1/10))
    assert np.allclose(chi.expected_freq,[[4.5,5.5],[4.5,5.5]])
    assert fisher.statistic == 36


def case_global_alignment():
    from examples.basic_algorithms_demo import global_alignment
    score, left, right = global_alignment('ACGT', 'AGT')
    print(score, left, right)
    rescored = sum(-2 if '-' in (a,b) else (2 if a == b else -1)
                   for a,b in zip(left,right))
    assert left.replace('-','') == 'ACGT'
    assert right.replace('-','') == 'AGT'
    assert score == rescored == 4
    print('recomputed_score', rescored)


def case_bfs():
    from collections import deque
    graph = {'A':['B','C'], 'B':['D'], 'C':['D'], 'D':['A'], 'Z':[]}
    distance, parent, queue = {'A':0}, {'A':None}, deque(['A'])
    while queue:
        node = queue.popleft()
        for neighbor in graph[node]:
            if neighbor not in distance:
                distance[neighbor] = distance[node] + 1
                parent[neighbor] = node
                queue.append(neighbor)
    print('distances', sorted(distance.items()))
    print('unreachable', sorted(set(graph) - set(distance)))
    assert distance['D'] == 2 and 'Z' not in distance


def case_glm():
    from math import exp, log
    beta0, beta1 = -2., log(2)
    def probability(x):
        return 1/(1+exp(-(beta0+beta1*x)))
    p0, p1 = probability(0), probability(1)
    OR = (p1/(1-p1))/(p0/(1-p0))
    print('probabilities',round(p0,6),round(p1,6))
    print('OR_RR_RD',round(OR,6),round(p1/p0,6),round(p1-p0,6))
    rate = .1
    print('Poisson_means_at_exposure_10_20', [round(exp(log(e)+log(rate)),6) for e in (10,20)])
    assert abs(OR-2)<1e-12 and p1/p0 < 2


def case_search_hash():
    from bisect import bisect_left, bisect_right
    from collections import Counter
    positions = [1, 3, 3, 8]
    left, right = bisect_left(positions, 3), bisect_right(positions, 3)
    print('bounds', left, right)
    print('closed_3_8', positions[left:bisect_right(positions, 8)])
    s, k = 'ATATAT', 2
    counts = Counter(s[i:i+k] for i in range(len(s)-k+1))
    print('overlapping_kmers', sorted(counts.items()))
    assert (left, right) == (1, 3)
    assert counts == {'AT': 3, 'TA': 2}


def case_anova():
    from statistics import mean
    from scipy.stats import f
    G = [[1,2,3],[3,4,5],[5,6,7]]
    all_y = [x for group in G for x in group]
    grand = mean(all_y)
    ssb = sum(len(g)*(mean(g)-grand)**2 for g in G)
    ssw = sum(sum((x-mean(g))**2 for x in g) for g in G)
    F = (ssb/(len(G)-1))/(ssw/(len(all_y)-len(G)))
    print('SSB_SSW_F_P',ssb,ssw,F,round(float(f.sf(F,2,6)),6))
    control, treated = [10.,12.], [13.,17.]
    interaction = (treated[1]-control[1])-(treated[0]-control[0])
    print('difference_of_effects',interaction)
    assert (ssb,ssw,F)==(24,6,12) and interaction==2


def case_gradient():
    def optimize(step, iterations=12):
        b, losses = 0., []
        for _ in range(iterations):
            losses.append((b-3)**2/2)
            b -= step*(b-3)
        return b, losses
    b, good = optimize(.5)
    _, bad = optimize(2.1)
    print('coefficient', round(b,6))
    print('good_loss', round(good[0],6), round(good[-1],6))
    print('bad_increases', bad[-1] > bad[0])
    assert abs(b-3) < .001
    assert all(v<=u for u,v in zip(good,good[1:]))
    assert bad[-1] > bad[0]


def case_diagnostics():
    from statistics import mean, variance
    low, high = [-1.,0.,1.], [-4.,0.,4.]
    print('residual_means',mean(low),mean(high))
    print('sample_variances',variance(low),variance(high))
    print('variance_ratio',variance(high)/variance(low))
    from scipy.stats import shapiro, levene
    s = shapiro(low+high)
    l = levene(low,high,center='median')
    print('Shapiro_P_small_n',round(float(s.pvalue),6))
    print('Levene_P_small_n',round(float(l.pvalue),6))
    assert variance(low)==1 and variance(high)==16


def case_confounding():
    import numpy as np
    x = np.array([1.,2.,3.,8.,9.,10.])
    y = np.array([11.,12.,13.,1.,2.,3.])
    group = np.array([0.,0.,0.,1.,1.,1.])
    simple = np.linalg.lstsq(np.column_stack([np.ones(6),x]),y,rcond=None)[0]
    adjusted = np.linalg.lstsq(np.column_stack([np.ones(6),x,group]),y,rcond=None)[0]
    print('pooled_slope',round(float(simple[1]),6))
    print('adjusted_intercept_x_group',np.round(adjusted,6))
    assert simple[1] < 0 and np.isclose(adjusted[1],1)


def case_rank_tests():
    from scipy.stats import mannwhitneyu, binomtest, wilcoxon
    A, B = [1,2,3], [4,5,6]
    UA = sum(a>b for a in A for b in B)+.5*sum(a==b for a in A for b in B)
    r = mannwhitneyu(A,B,alternative='two-sided',method='exact')
    print('U_A_exact_p', UA, r.pvalue)
    d = [1,2,-1,3,4]
    sign = binomtest(sum(v>0 for v in d),len(d),.5,alternative='two-sided')
    print('sign_test', sign.pvalue)
    w = wilcoxon([1,2,3,4],alternative='two-sided',method='exact')
    print('signed_rank', w.statistic, w.pvalue)
    assert UA == r.statistic == 0
    assert abs(r.pvalue-.1)<1e-12


def case_complexity():
    n = 1_000_000
    print('scan_checks', n)
    print('binary_upper_steps', n.bit_length())
    rows, columns, bytes_per_value = 10_000, 1_000_000, 8
    size = rows * columns * bytes_per_value
    print('dense_GB', size / 10**9)
    print('dense_GiB', round(size / 2**30, 2))
    small = [1, 3, 3, 8, 10]
    print('closed_interval', [x for x in small if 3 <= x <= 8])
    assert [x for x in small if 3 <= x <= 8] == [3, 3, 8]


def case_resampling():
    from itertools import combinations
    values = [1.,2.,3.,4.,5.,6.]
    observed = 3.
    null = []
    for Aidx in combinations(range(6),3):
        selected = set(Aidx)
        A = [x for i,x in enumerate(values) if i in selected]
        B = [x for i,x in enumerate(values) if i not in selected]
        null.append(sum(B)/3-sum(A)/3)
    p = sum(abs(t)>=abs(observed)-1e-12 for t in null)/len(null)
    print('allocations_exact_p', len(null), p)
    assert len(null)==20 and p==.1
    from random import Random
    from statistics import mean
    x, rng = [1.,2.,3.,4.], Random(20261026)
    boot = [mean(rng.choices(x,k=len(x))) for _ in range(2000)]
    print('bootstrap_draws', len(boot), 'mean_near', round(mean(boot),3))
    print('range', min(boot), max(boot))
    assert len(boot)==2000 and all(1<=v<=4 for v in boot)


def case_clustering():
    points, centers = [0.,1.,9.,10.], [0.,10.]
    labels = [min(range(2), key=lambda k:(x-centers[k])**2) for x in points]
    updated = [sum(x for x,c in zip(points,labels) if c==k)/labels.count(k)
               for k in range(2)]
    sse = sum((x-updated[c])**2 for x,c in zip(points,labels))
    print('labels', labels, 'centers', updated, 'SSE', sse)
    assert labels == [0,0,1,1] and sse == 1
    from examples.basic_algorithms_demo import kmeans
    import numpy as np
    C, labels2, history = kmeans(np.array(points)[:,None],2)
    print('final_SSE', history[-1], 'nonincreasing', bool((np.diff(history)<=1e-10).all()))
    assert np.isclose(history[-1],1)


def case_distance():
    from math import dist, sqrt
    A, B, C = (0,0), (1,1000), (3,0)
    print('raw', round(dist(A,B),3), dist(A,C))
    scale = (1,1000)
    def transform(x):
        return tuple(v/s for v,s in zip(x,scale))
    print('rescaled', round(dist(transform(A),transform(B)),3),
          dist(transform(A),transform(C)))
    x, y = (1,2), (2,4)
    cosine = sum(a*b for a,b in zip(x,y)) / sqrt(sum(a*a for a in x)*sum(b*b for b in y))
    print('proportional_vectors', dist(x,y), round(cosine,6))
    assert dist(A,B) > dist(A,C)
    assert dist(transform(A),transform(B)) < dist(transform(A),transform(C))


def case_classifiers():
    train = [(0.,'A'),(1.,'A'),(9.,'B'),(10.,'B')]
    neighbors = sorted(train,key=lambda row:abs(row[0]-8.5))[:3]
    prediction = max(('A','B'),key=lambda c:sum(label==c for _,label in neighbors))
    a_score, b_score = .5*.8*.7, .5*.2*.3
    posterior_a = a_score/(a_score+b_score)
    gini_parent = 1-(.5**2+.5**2)
    print('kNN', neighbors, prediction)
    print('NB_P_A', round(posterior_a,6))
    print('tree_Gini_gain', gini_parent-0)
    assert prediction == 'B' and abs(posterior_a-28/31)<1e-12
    assert gini_parent == .5


def case_em():
    from math import exp, log, pi
    x, means, weights = [0.,1.,4.,5.], [0.,5.], [.5,.5]
    def density(v,m):
        return exp(-.5*(v-m)**2)/(2*pi)**.5
    def likelihood(ms,ws):
        return sum(log(sum(w*density(v,m) for w,m in zip(ws,ms))) for v in x)
    responsibilities = []
    for v in x:
        a = [w*density(v,m) for w,m in zip(weights,means)]
        responsibilities.append([z/sum(a) for z in a])
    Nk = [sum(row[k] for row in responsibilities) for k in range(2)]
    new_weights = [n/len(x) for n in Nk]
    new_means = [sum(row[k]*v for row,v in zip(responsibilities,x))/Nk[k] for k in range(2)]
    print('left_responsibilities',[round(row[0],6) for row in responsibilities])
    print('new_means',[round(m,6) for m in new_means],'weights',new_weights)
    print('loglik_before_after',round(likelihood(means,weights),6),round(likelihood(new_means,new_weights),6))
    assert all(abs(sum(row)-1)<1e-12 for row in responsibilities)
    assert likelihood(new_means,new_weights) >= likelihood(means,weights)-1e-12


def case_hmm():
    from itertools import product
    pi = [.5,.5]
    A = [[.8,.2],[.2,.8]]
    emission = [{'A':.9,'B':.1},{'A':.2,'B':.8}]
    obs = 'AAB'
    delta = [pi[j]*emission[j][obs[0]] for j in range(2)]
    paths = [[j] for j in range(2)]
    for symbol in obs[1:]:
        new_delta, new_paths = [], []
        for j in range(2):
            i = max(range(2),key=lambda h:delta[h]*A[h][j])
            new_delta.append(delta[i]*A[i][j]*emission[j][symbol])
            new_paths.append(paths[i]+[j])
        delta, paths = new_delta, new_paths
    best = max(range(2),key=lambda j:delta[j])
    def joint(path):
        p = pi[path[0]]*emission[path[0]][obs[0]]
        for t in range(1,len(obs)):
            p *= A[path[t-1]][path[t]]*emission[path[t]][obs[t]]
        return p
    all_paths = list(product(range(2),repeat=len(obs)))
    print('Viterbi',paths[best],round(delta[best],6))
    print('total_observation_probability',round(sum(joint(p) for p in all_paths),6))
    assert paths[best] == [0,0,1]
    assert abs(delta[best]-max(map(joint,all_paths))) < 1e-12


def case_fm_index():
    T = 'banana$'
    SA = sorted(range(len(T)),key=lambda i:T[i:])
    BWT = ''.join(T[(i-1)%len(T)] for i in SA)
    C = {c:sum(x<c for x in T) for c in set(T)}
    def occ(c,k):
        return BWT[:k].count(c)
    l, r = 0, len(T)
    for c in reversed('ana'):
        l, r = C[c]+occ(c,l), C[c]+occ(c,r)
        print('prepend',c,'interval',l,r)
    hits = sorted(SA[l:r])
    print('SA',SA,'BWT',BWT,'positions',hits)
    assert BWT=='annb$aa' and hits==[1,3]
    assert hits==[i for i in range(len(T)) if T.startswith('ana',i)]


def case_debruijn():
    from collections import Counter
    sequence, k = 'ATGACG', 3
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    edges = [(word[:-1],word[1:]) for word in kmers]
    counts = Counter(edges)
    assembled = edges[0][0]+''.join(v[-1] for u,v in edges)
    print('kmers',kmers)
    print('edges',edges)
    print('assembled',assembled,'length',len(assembled))
    assert assembled==sequence
    branch_words = ['ATG','ATC']
    branches = {w[1:] for w in branch_words if w[:-1]=='AT'}
    print('AT_outneighbors',sorted(branches))
    assert len(branches)==2


CASES = {
    'blast': ('BLAST原理与种子延伸', case_blast),
    'pca': ('PCA原理与实现', case_pca),
    'smith-waterman': ('Smith-Waterman局部比对', case_smith_waterman),
    't-tests': ('t检验配对设计与Welch检验', case_t_tests),
    'mle-map': ('似然最大似然与MAP估计', case_mle_map),
    'hypothesis': ('假设检验与双侧单侧选择', case_hypothesis),
    'classification': ('分类指标ROC与PR曲线', case_classification),
    'contingency': ('列联表卡方与Fisher精确检验', case_contingency),
    'global-alignment': ('动态规划与序列比对基础', case_global_alignment),
    'bfs': ('图与最短路径的生信入门', case_bfs),
    'glm': ('广义线性模型与逻辑回归', case_glm),
    'search-hash': ('排序查找与哈希在生信中的应用', case_search_hash),
    'anova': ('方差分析交互作用与对比', case_anova),
    'gradient': ('梯度下降与数值优化', case_gradient),
    'diagnostics': ('正态性方差齐性与模型诊断', case_diagnostics),
    'confounding': ('相关回归与混杂', case_confounding),
    'rank-tests': ('秩检验符号检验与非参数方法', case_rank_tests),
    'complexity': ('算法输入输出复杂度与数据结构', case_complexity),
    'resampling': ('置换检验与Bootstrap实践', case_resampling),
    'clustering': ('聚类Kmeans与层次聚类', case_clustering),
    'distance': ('距离标准化与相似性', case_distance),
    'classifiers': ('近邻朴素贝叶斯与决策树入门', case_classifiers),
    'em': ('EM算法与隐变量估计', case_em),
    'hmm': ('隐马尔可夫模型与Viterbi解码', case_hmm),
    'fm-index': ('BWT与FM索引的精确匹配', case_fm_index),
    'debruijn': ('deBruijn图与短读长组装原理', case_debruijn),
}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--list', action='store_true', help='List case keys without running')
    parser.add_argument('--case', choices=sorted(CASES), help='Run one case; default runs all')
    parser.add_argument('--output', type=Path, default=PROJECT/'examples/output/methods-foundations')
    args = parser.parse_args()
    if args.list:
        for key, (page, _) in CASES.items():
            print(key + '\t' + page)
        return
    out = args.output.resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        parser.error('Output must be a new or empty directory; use --output to choose another')
    selected = [args.case] if args.case else list(CASES)
    rows = []
    for key in selected:
        page, function = CASES[key]
        buffer = StringIO()
        with redirect_stdout(buffer):
            function()
        rows.append(dict(case=key, page=page, stdout=buffer.getvalue(), checks='passed'))
    summary = dict(scope='Artificial arithmetic and toy algorithms; not production or real-data performance',
                   python=sys.version.split()[0], cases=len(rows), results=rows,
                   script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    out.mkdir(parents=True, exist_ok=True)
    (out/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
