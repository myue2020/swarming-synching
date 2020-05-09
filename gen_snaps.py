import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt

import sys
import os

idir=sys.argv[1]
odir='snap_pngs'

if not os.path.exists(odir):
    os.makedirs(odir)

def f(theta):
    return 0.45*(1+np.cos(theta))

def p2rgb(t):
    return (f(t), f(t - 2*np.pi/3), f(t + 2*np.pi/3))

snaps = [0,1800]

for snap in snaps:
    x = []
    y = []
    col = []

    data = pd.read_csv(idir+'/fr.'+str(snap), header=None, skiprows=[0], delim_whitespace=True, index_col=False)

    x = data[0].values
    y = data[1].values
    t = data[2].values

    col = [p2rgb(p) for p in t]

    fig = plt.figure()

    plt.scatter(x, y, c=col)

    plt.title('t=' + str(snap*600.0/1800) + ' Snapshot')
    plt.savefig(odir+'/snap_' + str(snap) + '.png')
