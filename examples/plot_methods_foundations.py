#!/usr/bin/env python3
"""Two original teaching figures: a global-alignment matrix and rank-one PCA.

Artificial inputs are identical to the corresponding methods chapters.
Requires NumPy and Matplotlib; output must be a new or empty directory.
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=PROJECT/'examples/output/methods-foundations-plots')
    args = parser.parse_args()
    out = args.output.resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        parser.error('Use a new or empty output directory; previous results are preserved')
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size':11, 'axes.spines.top':False, 'axes.spines.right':False})
    # Use a transparent recurrence, independently checked against the chapter table.
    a,b='ACGT','AGT'
    F=np.zeros((len(a)+1,len(b)+1),int)
    F[:,0]=-2*np.arange(len(a)+1);F[0,:]=-2*np.arange(len(b)+1)
    for i in range(1,len(a)+1):
        for j in range(1,len(b)+1):
            F[i,j]=max(F[i-1,j-1]+(2 if a[i-1]==b[j-1] else -1),F[i-1,j]-2,F[i,j-1]-2)
    assert F.tolist()==[[0,-2,-4,-6],[-2,2,0,-2],[-4,0,1,-1],[-6,-2,2,0],[-8,-4,0,4]]
    fig,ax=plt.subplots(figsize=(6.6,5.4),layout='constrained')
    ax.imshow(F,cmap='Blues',vmin=-9,vmax=7)
    for (i,j),v in np.ndenumerate(F):ax.text(j,i,str(v),ha='center',va='center',color='#1b2838',fontsize=13)
    # One valid optimal path; arrows follow traceback, not forward filling.
    path=[(4,3),(3,2),(2,1),(1,1),(0,0)]
    for (i,j),(ii,jj) in zip(path,path[1:]):
        ax.annotate('',xy=(jj,ii),xytext=(j,i),arrowprops={'arrowstyle':'->','color':'#d55e00','lw':2.2,'shrinkA':11,'shrinkB':11})
    ax.set_xticks(range(4),['empty',*b]);ax.set_yticks(range(5),['empty',*a])
    ax.set(xlabel='Prefix of sequence B',ylabel='Prefix of sequence A',
           title='Global alignment: match +2, mismatch -1, gap -2')
    ax.text(.5,-.21,'Traceback: ACGT / A-GT     Score = 3 x 2 - 2 = 4',transform=ax.transAxes,ha='center',fontsize=11)
    fig.savefig(out/'global-alignment.png',dpi=180,facecolor='white');plt.close(fig)
    X=np.array([[1.,2.],[2.,4.],[3.,6.]])
    mean=X.mean(axis=0);C=X-mean
    U,s,Vt=np.linalg.svd(C,full_matrices=False)
    direction=Vt[0].copy()
    if direction[0]<0:direction*=-1
    scores=C@direction
    assert np.allclose(scores,[-np.sqrt(5),0,np.sqrt(5)])
    fig,axes=plt.subplots(1,2,figsize=(10,4.7),layout='constrained',gridspec_kw={'width_ratios':[1.1,1]})
    ax=axes[0];ax.axhline(0,color='.85');ax.axvline(0,color='.85')
    ax.plot([-2,2],[-4,4],color='#0072b2',alpha=.4,label='PC1 direction')
    ax.scatter(C[:,0],C[:,1],s=65,c=['#0072b2','#009e73','#d55e00'],zorder=3)
    for i,(x,y) in enumerate(C):ax.annotate(f'S{i+1}',(x,y),xytext=(9,0),textcoords='offset points')
    ax.annotate('',xy=tuple(direction*3.5),xytext=(0,0),arrowprops={'arrowstyle':'->','color':'#1b2838','lw':2})
    ax.set(xlim=(-2,2),ylim=(-4.1,4.1),xlabel='Centered feature 1',ylabel='Centered feature 2',title='Center: subtract (2, 4)')
    ax.set_aspect('equal',adjustable='box')
    axes[1].axhline(0,color='.7',lw=1)
    axes[1].scatter(scores,np.zeros(3),s=75,c=['#0072b2','#009e73','#d55e00'],zorder=3)
    for i,v in enumerate(scores):axes[1].annotate(f'S{i+1}\n{v:.3f}',(v,0),xytext=(0,12),textcoords='offset points',ha='center')
    axes[1].set(xlim=(-3,3),ylim=(-1,1),xlabel='PC1 score',yticks=[],title='Project: 100% of variance on PC1')
    axes[1].spines[['left','right','top']].set_visible(False)
    axes[1].text(.5,.15,'A sign flip gives the same axis.\nReconstruction error is approximately zero.',transform=axes[1].transAxes,ha='center',fontsize=10)
    fig.savefig(out/'pca-projection.png',dpi=180,facecolor='white');plt.close(fig)
    print(out)


if __name__=='__main__':
    main()
