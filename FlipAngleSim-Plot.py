# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 11:45:26 2024

@author: p15094oj
"""

import numpy as np
import scipy
import matplotlib.pyplot as plt
import seaborn as sns

def make_TR_range(TR_min,TR_max,N):
    #TR = ms
    TR_min = TR_min/1000 #s
    TR_max = TR_max/1000 #s
    return(np.linspace(TR_min,TR_max,N))

def make_R10_range(T10_min,T10_max,N):
    #T1 = ms
    T10_min = T10_min/1000 #s
    T10_max = T10_max/1000 #s
    R10_min = 1/T10_min #s-1
    R10_max = 1/T10_max #s-1
    return(np.linspace(R10_min,R10_max,N))

def S(FA,R1,TR):
    #R1 = s-1
    #TR = s-1
    #FA = deg
    FA = FA * np.pi/180
    return(np.sin(FA)*(1-np.exp(-TR*R1)))/(1-(np.cos(FA)*np.exp(-TR*R1)))

def dS(FA,R10,TR,r1,C):
    #FA = deg
    #R1 = s-1
    #TR = s-1
    #r1 = s-1mM-1
    #C = mM
    dR1 = r1*C + R10
    return(S(FA,dR1,TR) - S(FA,R10,TR))

# TR = 0.003 #s
# R10 = 1 # s-1
# C = 0.01 #mM
# r1 = 3.4 #s-1mM-1
# FA_opt = scipy.optimize.fmin(lambda x: -dS(x,R10,TR,r1,C),10)

# dSig = []
# fa_list = [2,4,6,8,10,12,14,16,18,20]
# for fa in fa_list:
#     dSig.append(dS(fa,R10,TR,r1,C))

#%% Vary T10 and TR
# Constants
N = 10 # No samples
C_fixed = 0.01 #mM
r1 = 3.4 #s-1mM-1

# Range of T10 and TR
T10_min = 800 # ms - some NAWM
T10_max = 1500 # ms - some GM
R10_range = make_R10_range(T10_min,T10_max,N)
TR_min = 1 #ms
TR_max = 20 #ms
TR_range = make_TR_range(TR_min,TR_max,N)

FA_opt = np.zeros((N,N))#
delta_sig = np.zeros((N,N))
TR_grid = np.zeros((N,N))
R10_grid = np.zeros((N,N))
# for all TR and T10
for i in range(N):
    for j in range(N): 
        R10 = R10_range[i]
        TR = TR_range[j]
        # find FA that will maximise dS
        FA_opt_temp = scipy.optimize.fmin(lambda x: -dS(x,R10,TR,r1,C_fixed),0,full_output=True)
        FA_opt[i,j] = FA_opt_temp[0]
        delta_sig[i,j] = -FA_opt_temp[1]
        TR_grid[i,j] = TR
        R10_grid[i,j] = R10

xlabs = np.char.mod('%.0f',TR_range*1000).tolist() # TR in ms
ylabs = np.char.mod('%.0f',(1/R10_range)*1000).tolist() # T10 in ms
fig = plt.figure()
ax = sns.heatmap(FA_opt,square=True,cmap='jet',annot=True,annot_kws={'font':'Arial','fontsize':12},cbar_kws={'label':'\nOptimal Flip Angle [degrees]'},xticklabels=xlabs,yticklabels=ylabs)
ax.set_xlabel('TR [ms]')
ax.set_ylabel('T$_{10}$ [ms]')
ax.tick_params(bottom=True, left=True)
plt.title('Optimal Flip Angle for Maximum Sensitivity to Change in T1')
plt.savefig('C:/Users/p15094oj/Documents/Python_scripts/FlipAngle_Opt/OptFA_T10vTR.png')
plt.show()

ax = sns.heatmap(delta_sig,square=True,cmap='jet',annot=False,annot_kws={'font':'Arial','fontsize':12},cbar_kws={'label':'\nSignal change at optimal flip angle'},xticklabels=xlabs,yticklabels=ylabs)
ax.set_xlabel('TR [ms]')
ax.set_ylabel('T$_{10}$ [ms]')
ax.tick_params(bottom=True, left=True)
plt.title('Optimal Flip Angle for Maximum Sensitivity to Change in T1')
plt.savefig('C:/Users/p15094oj/Documents/Python_scripts/FlipAngle_Opt/OptFA_delta_sig.png')
plt.show()

#%% Vary C and TR
# Constants
N = 10 # No samples
R10_fixed = 1 #ms
r1 = 3.4 #s-1mM-1

# Range of C and TR
C_range = np.linspace(0.001,0.01,N)
TR_min = 1 #ms
TR_max = 20 #ms
TR_range =  make_TR_range(TR_min,TR_max,N)

FA_opt = np.zeros((N,N))
TR_grid = np.zeros((N,N))
C_grid = np.zeros((N,N))
# for all TR and T10
for i in range(N):
    for j in range(N): 
        C = C_range[i]
        TR = TR_range[j]
        # find FA that will maximise dS
        FA_opt_temp = scipy.optimize.fmin(lambda x: -dS(x,R10_fixed,TR,r1,C),0,full_output=True)
        FA_opt[i,j] = FA_opt_temp[0]
        TR_grid[i,j] = TR
        C_grid[i,j] = C

xlabs = np.char.mod('%.0f',TR_range*1000).tolist() # TR in ms
ylabs = np.char.mod('%.3f',C_range).tolist() # T10 in ms
fig = plt.figure()
ax = sns.heatmap(FA_opt,square=True,cmap='jet',annot=True,annot_kws={'font':'Arial','fontsize':12},cbar_kws={'label':'\nOptimal Flip Angle [degrees]'},xticklabels=xlabs,yticklabels=ylabs)
ax.set_xlabel('TR [ms]')
ax.set_ylabel('C [mM]')
ax.tick_params(bottom=True, left=True)
plt.title('Optimal Flip Angle for Maximum Sensitivity to Change in T1')
plt.savefig('C:/Users/p15094oj/Documents/Python_scripts/FlipAngle_Opt/OptFA_CvTR.png')
plt.show()

#%% Max Sig
# Constants
N = 10 # No samples
# Range of T10 and TR
T10_min = 800 # ms - some NAWM
T10_max = 1500 # ms - some GM
R10_range = make_R10_range(T10_min,T10_max,N)
TR_min = 1 #ms
TR_max = 20 #ms
TR_range = make_TR_range(TR_min,TR_max,N)

FA_opt = np.zeros((N,N))
delta_sig = np.zeros((N,N))
TR_grid = np.zeros((N,N))
R10_grid = np.zeros((N,N))
# for all TR and T10
for i in range(N):
    for j in range(N): 
        R10 = R10_range[i]
        TR = TR_range[j]
        # find FA that will maximise dS
        FA_opt_temp = scipy.optimize.fmin(lambda x: -S(x,R10,TR),0,full_output=True)
        FA_opt[i,j] = FA_opt_temp[0]
        delta_sig[i,j] = -FA_opt_temp[1]
        TR_grid[i,j] = TR
        R10_grid[i,j] = R10

xlabs = np.char.mod('%.0f',TR_range*1000).tolist() # TR in ms
ylabs = np.char.mod('%.0f',(1/R10_range)*1000).tolist() # T10 in ms
fig = plt.figure()
ax = sns.heatmap(FA_opt,square=True,cmap='jet',annot=True,annot_kws={'font':'Arial','fontsize':12},cbar_kws={'label':'\nOptimal Flip Angle [degrees]'},xticklabels=xlabs,yticklabels=ylabs)
ax.set_xlabel('TR [ms]')
ax.set_ylabel('T$_{10}$ [ms]')
ax.tick_params(bottom=True, left=True)
plt.title('Optimal Flip Angle for Maximum Signal')
plt.savefig('C:/Users/p15094oj/Documents/Python_scripts/FlipAngle_Opt/OptFA_maxsig_T10vTR.png')
plt.show()

ax = sns.heatmap(delta_sig,square=True,cmap='jet',annot=False,annot_kws={'font':'Arial','fontsize':12},cbar_kws={'label':'\nMaximum signal at optimal flip angle'},xticklabels=xlabs,yticklabels=ylabs)
ax.set_xlabel('TR [ms]')
ax.set_ylabel('T$_{10}$ [ms]')
ax.tick_params(bottom=True, left=True)
plt.title('Optimal Flip Angle for Maximum Signal')
plt.savefig('C:/Users/p15094oj/Documents/Python_scripts/FlipAngle_Opt/OptFA_maxsig.png')
plt.show()