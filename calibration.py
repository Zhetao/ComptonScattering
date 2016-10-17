# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from math import *
import os
import numpy as np
from scipy.optimize import curve_fit

def energyConv(x, a, b):
    return a*x+b

script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
f = open(os.path.join(script_dir,'CalibrationData.txt'),'r');
caList = []

for line in f:
    temp = line.rstrip().split();
    temp = [float(elem) for elem in temp];
    caList.append(temp);
#print caList

yList = [term[0] for term in caList];
xList = [np.mean(term[1:]) for term in caList];
stdList = [np.std(term[1:]) for term in caList];
pointY = [term[0]  for term in caList for i in range(5)]
pointX = [term[i]  for term in caList for i in range(1,6)]

popt, pcov = curve_fit(energyConv, pointX, pointY);
[a, b] = popt;
perr = np.sqrt(np.diag(pcov))
#Plot data
plt.clf();
fig = plt.figure()
fig.suptitle('Calibration: Energy-Channel Number Relationship', fontsize=14, fontweight='bold')

x = np.arange(1,16000, 0.2);
plt.plot(x, a*x+b);
plt.xlabel('Channel Number')
plt.ylabel('Energy(keV)')
plt.errorbar(xList,yList, xerr=stdList, capsize =6 ,ecolor='r',fmt = 'none');
plt.show();
plt.savefig(os.path.join(script_dir,'calibration.png'))

print popt, perr
f.close();