#import os, sys
#sys.path.insert(0, 'structural_variants/lib/')
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.io import savemat
#import sympy
import random


def generate_diploid_data(params, prnt=True, seed=None):
    '''
    Generate simulated data for a one parent, one child Structural Variant (SV) analysis
    Args: A dictionary containing the following parameters as keys
        n: size of data vectors (signals)
        k: total number of structural variants
        pctNovel: percent of novel structural variants in [0,1] (biological reality- very small %)
        lambda_p, lambda_c: sequence coverage of child and parent, respectively
        erreps: error (>0) incurred by sequencing and mapping process
        r: dispersion parameter for Negative Binomial distribution; we use 1 to maximize variance under Neg Bin
    
    The output is provided for different individuals. 
    Here, p2 = parent 2, p = parent 1, c = child, h = inherited, and n = novel. 
    We will calculate necessary variables for parent 1 (p) and the child's inherited (h) and novel (n) 
    SV's, and only use parent 2 (p2) for the logical implementation of inheritance of the child's signals. 
    
    Output: A dictionary containing the following data elements as keys
        Input parameters
        f_{p2,c,p}: nx1 true signal of an individual
        z_{p,h,c}:  nx1 indicator vector of homogeneous structural variants
        y_{p,h,c}:  nx1 indicator vector of heterogeneous structural variants
        mu_{p,c}:   Mean sequence coverage; mu_{} = A_{} * f_{};
        var_{p,c}:  Variance of sequence coverage; var_{} = mu_{} + (1/r) mu_{}^2, where r=1
        s_{p,c}:    nx1 random vector drawn from Negative binomial distribution
        A_z,{c,p}:   (lambda_c - erreps) I_n  
        A_y,{c,p}: (2*lambda_c - erreps) I_n
                    Note: The matrices A_{} are sparse diagonal nxn matrix and I_n is nxn identity matrix
        
    '''
    np.random.seed(seed=seed)
    if prnt: print('Using seed {} \nGenerating data...'.format(seed))
    # First, we randomly permute a sequence (0,1,2,3,...n-1) to be used for indices
    q = np.random.permutation(np.arange(params['n']))
    similarity = int(params['pct_similarity']* params['k']) # pct_similarity * number of SVs

    # Create signal variables and initialize signal vectors with all zeros
    signals = ['f_p2', 'f_c']
    for i,letter in enumerate(['p', 'h', 'n']):
        signals.append('f_%s'%letter)
        signals.append('z_%s'%letter)
        signals.append('y_%s'%letter)

    d= {}
    for signal in signals: d[signal] = np.zeros((params['n'],1), dtype=np.int32)

    # Parent signals: 
    #        f_p  - k elements will be 1s and 2s randomly selected
    #        f_p2 - floor of %similarity*k elements will be the same as f_p and the rest will be random 1s and 2s
    #               (i.e. the parents will only share a given percentage of SV's)

    # Insert k number of 1's and 2's in the first parent
    for i in q[:params['k']]: d['f_p'][i] = np.random.randint(1,3) 

    # parent 2 shares similarity number of SV's with parent 1, 
    # q[0:similarity] is the positions for which both parents will share SV's
    # the remaining k - similarity positions will be chosen randomly as to not overlap with existing SVs
    rand_choices = np.random.choice(np.setdiff1d(q, q[0:similarity]), size=params['k'] -similarity, replace= False); #rand_choices.sort(); q[0:similarity].sort()
    d['f_p2'][q[0:similarity]] = d['f_p'][q[0:similarity]]
    d['f_p2'][rand_choices]= np.random.randint(1,3) 

    # verify that both parents have the same amount of nonzero entries
    if np.count_nonzero(d['f_p2']) != np.count_nonzero(d['f_p']): 
        print('PARENTS DO NOT HAVE EQUAL AMOUNT OF NONZERO ENTRIES !!')
        print('f_p: {} \nf_p2: {}'.format(np.count_nonzero(d['f_p']),np.count_nonzero(d['f_p2'])))


    # Child signal
    #     First, we defined the inherited SVs through a logical implementation 
    #     of inheritance using parent 1 (f_p) and parent 2 (f_p2)
    for i in np.arange(d['f_p'].shape[0]):
        if   (d['f_p'][i]==2 and d['f_p2'][i]==2): d['f_h'][i]= 2
        elif (d['f_p'][i]==1 and d['f_p2'][i]==1): d['f_h'][i]= np.random.randint(0,3)
        elif (d['f_p'][i]==2 and d['f_p2'][i]==0) or (d['f_p'][i]==0 and d['f_p2'][i]==2): d['f_h'][i]= 1
        elif (d['f_p'][i]==2 and d['f_p2'][i]==1) or (d['f_p'][i]==1 and d['f_p2'][i]==2): d['f_h'][i]= np.random.randint(1,3)
        elif (d['f_p'][i]==1 and d['f_p2'][i]==0) or (d['f_p'][i]==0 and d['f_p2'][i]==1): d['f_h'][i]= np.random.randint(0,2)

    #    Next, we define the novel SVs 
    #    define inherited indices and novel indices, make sure they do not overlap
    inherited_pos = d['f_h'].nonzero()[0]; inherited_pos.sort()
    novel_pos = np.random.choice(np.setdiff1d(q, inherited_pos), size=int(params['k'] *params['pctNovel']), replace= False)
    d['f_n'][novel_pos]= np.random.randint(1,3) 

    #    Lastly, we define the complete child signal
    d['f_c'] = d['f_h'] + d['f_n']

    # convert signals to indicators
    for letter in ['p','h','n']:
        for i in np.arange(d['f_c'].shape[0]):
            if   d['f_%s'%letter][i]==2: d['z_%s'%letter][i]=1
            elif d['f_%s'%letter][i]==1: d['y_%s'%letter][i]=1
    
    # define coverage matrices, mean, variance, and observation signal for parent and child
    for i, letter in enumerate(['p','c']):
        d['A_z%s'%letter]   = (2*params['lambda_%s'%letter] - params['erreps'])*sparse.eye(params['n'])
        d['A_y%s'%letter]   = (params['lambda_%s'%letter] - params['erreps'])*sparse.eye(params['n'])
        d['mu_%s'%letter]  = np.matmul((d['A_z%s'%letter]+d['A_y%s'%letter]).toarray(), d['f_%s'%letter]) + params['erreps']
        d['var_%s'%letter] = d['mu_%s'%letter] +(1/params['r'])*(d['mu_%s'%letter]**2)
        d['s_%s'%letter]   = np.random.negative_binomial(d['mu_%s'%letter]/(d['var_%s'%letter]-d['mu_%s'%letter]),d['mu_%s'%letter]/d['var_%s'%letter])

    data = {**d, **params}
    if prnt:
        print('Done!')
        print()
        print('Using parameters:')
        for key, val in params.items():
            print('\t', key, ': ', val)  
    return data

def generate_haploid_data(params):
    '''
    Generate simulated data for a one parent, one child Structural Variant analysis
    Args: A dictionary containing the following parameters as keys
        n: size of data vectors (signals)
        k: total number of structural variants
        pctNovel: percent of novel structural variants in [0,1] (biological reality- very small %)
        lambda_p, lambda_c: sequence coverage of child and parent, respectively
        erreps: error (>0) incurred by sequencing and mapping process
        r: dispersion parameter for Negative Binomial distribution
    
    Output: A dictionary containing the following data elements as keys
        A_c: (lambda_c - erreps) I_n, sparse diagonal nxn matrix. I_n is nxn identity matrix
        A_p: (lambda_p - erreps) I_n, sparse diagonal nxn matrix. I_n is nxn identity matrix
        mu_p, var_p: Mean and variance sequence coverage for parent; mu_p = A_p * f_p 
        mu_c, var_c: Mean and variance sequence coverage for child;  mu_c = A_c * f_c 
        s_p: nx1 random vector drawn from Negative binomial distribution (for parent)
        s_c: nx1 random vector drawn from Negative binomial distribution (for child)
        TODO: add mu and var
        for i in {P (parent), H (inherited), N (novel)}:
        z_i: nx1 indicator vector of homogeneous structural variants
        y_i: nx1 indicator vector of heterogeneous structural variants
        
    '''
    q = np.random.permutation(params['n'])
    #print(q)
    startVal = int(params['k']*params['pctNovel']); #print(startVal)
    endVal = int(startVal +params['k']) ; #print(endVal)

    f_p, f_c, f_h, f_n = np.zeros((params['n'],1), dtype=np.int8),np.zeros((params['n'],1), dtype=np.int8), np.zeros((params['n'],1), dtype=np.int8), np.zeros((params['n'],1), dtype=np.int8)
    f_p[q[: params['k']]], f_c[q[startVal:endVal]] = 1,1
    f_h[q[startVal:params['k']+1]], f_n[q[params['k']+1:endVal]] = 1,1
    
    
    d = {}
    d['f_p'] = f_p; d['f_h'] = f_h; d['f_n'] = f_n; d['f_c'] = f_h + f_n; 
    
    for i, letter in enumerate(['p','c']):
        d['A_%s'%letter]   = (params["lambda_%s"%letter] - params['erreps'])*sparse.eye(params['n'])
        d['mu_%s'%letter]  = np.matmul(d['A_%s'%letter].toarray(), d['f_%s'%letter]) + params['erreps']
        d['var_%s'%letter] = d['mu_%s'%letter] +(1/params['r'])*(d['mu_%s'%letter]**2)
        d['s_%s'%letter]   = np.random.negative_binomial(d['mu_%s'%letter]/(d['var_%s'%letter]-d['mu_%s'%letter]),d['mu_%s'%letter]/d['var_%s'%letter])
    
    
    return d