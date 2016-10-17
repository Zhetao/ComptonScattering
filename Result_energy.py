# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from math import *
import os
import numpy as np
from scipy.optimize import curve_fit
from uncertainties.umath import *
from uncertainties import ufloat

Omega = pi*(7.95*0.01/2)**2;
n = 9.22*10**20; #2.7*10^3*(pi*0.005^2*0.015)*13/27*Na
#[ 0.08564978  3.87254083] [  1.72473518e-04   1.63431991e+00]
A= 0.08564978;
B = 3.87254083;
dA = 1.72473518e-04;
dB = 1.63431991e+00;

uA = ufloat(0.08564978, 1.72473518e-04)
uB = ufloat(3.87254083, 1.63431991e+00)


def convFunc(x):
    return A*x+B;

def uconvFunc(x):
    return uA*x+uB;

script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
f = open(os.path.join(script_dir,'Fitting.txt'),'r');
dataList = [];

for line in f:
    temp = line.rstrip().split();
    temp = [float(elem) for elem in temp];
    dataList.append(temp);

pErrList = [];
nErrList = [];
for term in dataList:
    if term[0] == 0:
        zeroChannel = term[1];
        zeroDiff = term[2];
    if term[0] > 0:
        pErrList.append(sqrt(term[3]**2*A**2+term[1]**2*dA**2+dB**2));
    if term[0] < 0:
        nErrList.append(sqrt(term[3]**2*A**2+term[1]**2*dA**2+dB**2));        

pAngleList = [abs(term[0]) for term in dataList if term[0] > 0];
pChanList = [term[1] for term in dataList if term[0] > 0];
pmuErrList = [term[3] for term in dataList if term[0] > 0];
nAngleList = [abs(term[0]) for term in dataList if term[0] < 0];
nChanList = [term[1] for term in dataList if term[0] < 0];
nmuErrList = [term[3] for term in dataList if term[0] < 0];
diffList = [term[2] for term in dataList];

upEnergyList = [];
unEnergyList = [];
for i in range(len(pChanList)):
    upChannel = ufloat(pChanList[i], pmuErrList[i]);  
    upEnergyList.append(uconvFunc(upChannel))
    
for i in range(len(pChanList)):
    unChannel = ufloat(nChanList[i], nmuErrList[i]);  
    unEnergyList.append(uconvFunc(unChannel))
print 'new method';
print upEnergyList, unEnergyList;   

pEnergyList = [convFunc(channel) for channel in pChanList];
nEnergyList = [convFunc(channel) for channel in nChanList];
print 'origin'
print pErrList,nErrList

plt.clf();
fig = plt.figure()
fig.suptitle('Scattered Photon Energy-Angle relationship', fontsize=14, fontweight='bold')

x1 = np.arange(0, 90, 0.2);
plt.plot(x1, 1/(1.0/661+1.954*10**(-3)*(1-np.cos(x1/180*np.pi))), label = 'theoretical value');
posi = plt.errorbar(pAngleList,pEnergyList, yerr=pErrList,ecolor='r',fmt = 'none', label = 'positive angle');
nega = plt.errorbar(nAngleList,nEnergyList, yerr=nErrList,ecolor='g',fmt = 'none', label = 'negative angle');
plt.legend(handles=[posi, nega])
plt.show();
plt.savefig(os.path.join(script_dir,'energy.png'))
