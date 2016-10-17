# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from math import *
import os
import numpy as np
from scipy.optimize import curve_fit
from uncertainties.umath import *
from uncertainties import ufloat
import codecs
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

uEnergyList = [];    
AngleList = [term[0] for term in dataList if term[0] != 0];
ChanList = [term[1] for term in dataList if term[0] != 0];
muErrList = [term[3] for term in dataList if term[0] != 0]; 
 
for i in range(len(ChanList)):
    uChannel = ufloat(ChanList[i], muErrList[i]);  
    uEnergyList.append(uconvFunc(uChannel))
print uEnergyList

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

f3 = codecs.open(os.path.join(script_dir,'DetectionEff.txt'), encoding='utf-16');
#f2 = open(,'r');

energyArr = [];
effArr = [];

for line in f3:
    temp = line.split();
    energyArr.append(float(temp[0]));
    effArr.append(float(temp[1]));

diffList = [];
for term in dataList:
    if term[0] == 0:
        zeroChannel = term[1];
        zeroDiff = term[2];
    if term[0] != 0:
        diffList.append(ufloat(term[2],term[4]));   
zeroIdx = findClosest(energyArr,convFunc(zeroChannel)/1000000.0);
zeroDiff = zeroDiff/effArr[zeroIdx];
dzeroDiff = sqrt(zeroDiff);

AngleList = [term[0] for term in dataList if term[0] != 0];
ChanList = [term[1] for term in dataList if term[0] != 0];
areaList = [];
ErrList = [];
for i in range(len(diffList)):
    tempIdx = findClosest(energyArr,convFunc(ChanList[i])/1000000.0);
    eps = effArr[tempIdx];
    ueps = ufloat(eps, 0.15 * eps);
    crossSection = diffList[i]/zeroDiff/Omega/n/ueps
    areaList.append(crossSection.nominal_value*10**30);
    ErrList.append(crossSection.std_dev*10**30);
print areaList, ErrList
#res = open(os.path.join(script_dir,'CombinedResult.txt'),'w');
egyAna = open(os.path.join(script_dir,'EnergyAnalysis.txt'),'w');
csAna = open(os.path.join(script_dir,'CSAnalysis.txt'),'w');

for i in range(len(AngleList)): 
#    res.write(str(AngleList[i]) + ' ' +'$'+ '%.1f' %(uEnergyList[i].nominal_value)+  '\pm' + '%.1f' %(uEnergyList[i].std_dev)+'$'+
#    ' '+'$'+  '%.2f' %(areaList[i])+  '\pm' +  '%.2f' %(ErrList[i])+'$'+ '\n');#*mergeSize) + '\n');
    thValue = 1/(1.0/662+1.954*10**(-3)*(1-np.cos(AngleList[i]/180*np.pi)));
    sigmaAway = abs(uEnergyList[i].nominal_value-thValue)/uEnergyList[i].std_dev;
    print sigmaAway;
    egyAna.write(str(AngleList[i]) + ' ' +'$'+ '%.1f' %(uEnergyList[i].nominal_value)+  '\pm' + '%.1f' %(uEnergyList[i].std_dev)+'$'+
    ' '+ str(int(thValue))+ ' ' +'%.1f' %sigmaAway+ '\n');

for i in range(len(AngleList)): 
#    res.write(str(AngleList[i]) + ' ' +'$'+ '%.1f' %(uEnergyList[i].nominal_value)+  '\pm' + '%.1f' %(uEnergyList[i].std_dev)+'$'+
#    ' '+'$'+  '%.2f' %(areaList[i])+  '\pm' +  '%.2f' %(ErrList[i])+'$'+ '\n');#*mergeSize) + '\n');
    x = AngleList[i];

    cos = np.cos(x/180*np.pi);
    x1 = 1/(1+1.294*(1-cos))
    thValue = 3.967*x1**2*(x1+1.0/x1-1+cos**2);
    sigmaAway = abs(areaList[i]-thValue)/ErrList[i];
    print sigmaAway;
    csAna.write(str(AngleList[i]) + ' ' +'$'+ '%.1f' %(areaList[i])+  '\pm' + '%.1f' %(ErrList[i])+'$'+
    ' '+ '%.1f' %thValue+ ' ' +'%.1f' %sigmaAway+ '\n');