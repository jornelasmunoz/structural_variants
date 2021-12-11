import os, sys
sys.path.insert(0, 'structural_variants/lib/')
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
from scipy import sparse
#import sympy
import random


def generate_diploid_data(params):
    q = np.random.permutation(params['n'])
    startVal = int(params['k']*params['pctNovel']); #print(startVal)
    endVal = int(startVal +params['k']) ; #print(endVal)
    similarity = int(params['pct_similarity']* params['k']) # pct_similarity * number of SVs

    signals = ['f_p2', 'f_c']
    for i,letter in enumerate(['p', 'h', 'n']):
        signals.append('f_%s'%letter)
        signals.append('z_%s'%letter)
        signals.append('y_%s'%letter)

    # Initialize signals
    d= {}
    for signal in signals: d[signal] = np.zeros((params['n'],1), dtype=np.int32)

    # parent signals: 
    #        f_p  - k elements will be 1s and 2s randomly selected
    #        f_p2 - floor of %similarity*k elements will be the same as f_p and the rest will be random 1s and 2s
    for i in q[:params['k']]: d['f_p'][i] = np.random.randint(1,3) 
    d['f_p2'][q[0:similarity]] = d['f_p'][q[0:similarity]]
    d['f_p2'][random.choices(q, k=params['k'] -similarity)]= np.random.randint(1,3) 

    # child signal
    #     inherited
    for i in np.arange(d['f_p'].shape[0]):
        if   (d['f_p'][i]==2 and d['f_p2'][i]==2): d['f_c'][i]= 2
        elif (d['f_p'][i]==1 and d['f_p2'][i]==1): d['f_c'][i]= np.random.randint(0,3)
        elif (d['f_p'][i]==2 and d['f_p2'][i]==0) or (d['f_p'][i]==0 and d['f_p2'][i]==2): d['f_c'][i]= 1
        elif (d['f_p'][i]==2 and d['f_p2'][i]==1) or (d['f_p'][i]==1 and d['f_p2'][i]==2): d['f_c'][i]= np.random.randint(1,3)
        elif (d['f_p'][i]==1 and d['f_p2'][i]==0) or (d['f_p'][i]==0 and d['f_p2'][i]==1): d['f_c'][i]= np.random.randint(0,2)
        
    #     novel        
    d['f_n'][random.choices(q, k=params['k'] -similarity)]= np.random.randint(1,3) 
    d['f_c'] = d['f_h'] +d ['f_n']
    
    # convert signals to indicators
    for j,letter in enumerate(['p','h','n']):
        for i in np.arange(d['f_c'].shape[0]):
            if   d['f_%s'%letter][i]==2: d['z_%s'%letter][i]=1
            elif d['f_%s'%letter][i]==1: d['y_%s'%letter][i]=1
            
    for i, letter in enumerate(['p','c']):
        d['A_%s'%letter]   = (params['lambda_%s'%letter] - params['erreps'])*sparse.eye(params['n'])
        d['mu_%s'%letter]  = np.matmul(d['A_%s'%letter].toarray(), d['f_%s'%letter]) + params['erreps']
        d['var_%s'%letter] = d['mu_%s'%letter] +(1/params['r'])*(d['mu_%s'%letter]**2)
        d['s_%s'%letter]   = np.random.negative_binomial(d['mu_%s'%letter]/(d['var_%s'%letter]-d['mu_%s'%letter]),d['mu_%s'%letter]/d['var_%s'%letter])
        
    return d
params = {
    'r': 1,
    'n': 10**2,
    'k': 10,
    'lambda_c': 4,
    'lambda_p': 8,
    'pctNovel': 0.15,
    'erreps'  : 1e-2,
    #'suffix'  : ['p','c'],
    'pct_similarity': 0.6}
data =generate_diploid_data(params)