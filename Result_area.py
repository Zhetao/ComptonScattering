# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from math import *
import os
import numpy as np
from scipy.optimize import curve_fit
import codecs
from uncertainties.umath import *
from uncertainties import ufloat

def findClosest(numList, goal):
    tempList = [term-goal for term in numList];
    idx = (np.abs(tempList)).argmin()
    #print goal, effArr[idx];
    return idx
r = ufloat(7.95*0.01/2,0.05*0.01)
R = ufloat(39.4*0.01, 0.5*0.01)
Omega = pi*r**2/R**2;
L = ufloat(0.015, 0.001)
n = 2.7*10**3*13/27*1000*6.02*10**23*L; #2.7*10^3*13/27*1000*Na*0.015


def convFunc(x):
    return (0.08564978*x + 3.87254083)*1000

script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
f = open(os.path.join(script_dir,'Fitting.txt'),'r');
f2 = codecs.open(os.path.join(script_dir,'DetectionEff.txt'), encoding='utf-16');
f3 = open(os.path.join(script_dir,'DetectorEff.txt'),'w');
#f2 = open(,'r');
dataList = [];
energyArr = [];
effArr = [];

for line in f2:
    temp = line.split();
    energyArr.append(float(temp[0]));
    effArr.append(float(temp[1]));
    f3.write(temp[0]+' '+temp[1]+'\n')

for line in f:
    temp = line.rstrip().split();
    temp = [float(elem) for elem in temp];
    dataList.append(temp);

pdiffList =[];
ndiffList =[];
for term in dataList:
    if term[0] == 0:
        zeroChannel = term[1];
        zeroDiff = term[2];
    if term[0] > 0:
        pdiffList.append(ufloat(term[2],term[4]));
    if term[0] < 0: 
        ndiffList.append(ufloat(term[2],term[4]));   
            
zeroIdx = findClosest(energyArr,convFunc(zeroChannel)/1000000.0);
zeroDiff = zeroDiff/effArr[zeroIdx];
dzeroDiff = sqrt(zeroDiff);

pAngleList = [abs(term[0]) for term in dataList if term[0] > 0];
pChanList = [term[1] for term in dataList if term[0] > 0];
nAngleList = [abs(term[0]) for term in dataList if term[0] < 0];
nChanList = [term[1] for term in dataList if term[0] < 0];
#angleList = [abs(term[0]) for term in dataList if term[0] != 0];
#pdiffList = [term[2] for term in dataList if term[0] > 0];
#ndiffList = [term[2] for term in dataList if term[0] < 0];

#pWavelenList = [1240.0/convFunc(channel) for channel in pChanList];
#nWavelenList = [1240.0/convFunc(channel) for channel in nChanList];
pareaList = [];
nareaList = [];
pErrList = [];
nErrList = [];
for i in range(len(pdiffList)):
    tempIdx = findClosest(energyArr,convFunc(pChanList[i])/1000000.0);
    eps = effArr[tempIdx];
    ueps = ufloat(eps, 0.15 * eps);
    crossSection = pdiffList[i]/zeroDiff/Omega/n/ueps
    pareaList.append(crossSection.nominal_value);
    pErrList.append(crossSection.std_dev);
for i in range(len(ndiffList)):
    tempIdx = findClosest(energyArr,convFunc(nChanList[i])/1000000.0);
    eps = effArr[tempIdx];
    ueps = ufloat(eps, 0.15 * eps);
    crossSection = ndiffList[i]/zeroDiff/Omega/n/ueps
    nareaList.append(crossSection.nominal_value);
    nErrList.append(crossSection.std_dev);
print pareaList

plt.clf();
fig = plt.figure()
fig.suptitle('Differential cross section-scattering angle relationship', fontsize=14, fontweight='bold')
x = np.arange(0, 100, 0.2);
g = 1.293
cos = np.cos(x/180*np.pi);
y = (2.82*10**-15)**2*(1+cos**2)/2.0*(1/((1+g*(1+cos)))**2)*(1+(g**2*((1-cos)**2)/((1+cos**2)*(1+g*(1-cos)))))

x1 = 1/(1+1.294*(1-cos))
y1 = 3.967*10**-30*x1**2*(x1+1.0/x1-1+cos**2);
plt.plot(x, y1);
plt.xlabel('Scattering angle')
plt.ylabel('Differential Cross Section')
posi = plt.errorbar(pAngleList,pareaList, yerr=pErrList, capsize =6,ecolor='r',fmt = 'none', label = 'positive angle');
nega = plt.errorbar(nAngleList,nareaList, yerr=nErrList, capsize =6,ecolor='g',fmt = 'none', label = 'negative angle');
plt.legend(handles=[posi, nega])
plt.show();
plt.xlim([-5,100])
plt.ylim([-max(pareaList),max(pareaList)*2.3])
plt.savefig(os.path.join(script_dir,'crossSection.png'))