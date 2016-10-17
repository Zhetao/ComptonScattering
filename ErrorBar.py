# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from math import *
import os
import numpy as np
from scipy.optimize import curve_fit

def noisePoly(x, a, b, c):
    return a*x**2+b*x+c

def sigGaussian(x, a0, miu, sig):
    return a0*np.exp(-(x-miu)**2/(2*sig**2))


##Read Data and Write Data 
script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
timeList = {0:600, 10:300,20:600,30:800,40:1300,50:1500,60:1800,70:2400,80:7200,90:10800}
sigTime = {0:[140,145,175,195],-10:[90,130,170,230],-20:[70,120,170,230],-30:[60,105,160,220],-40:[60,112,150,210],-50:[45,90,130,220],-60:[45,75,110,200],-70:[25,60,100,150],-80:[25,58,92,135],-90:[25,52,80,140],
                10:[90,130,170,230],20:[70,130,175,230],30:[60,120,165,220],40:[60,100,155,210],50:[45,90,130,220],60:[45,75,115,200]}
f = open(os.path.join(script_dir,'Fitting.txt'),'w');
f2 = open(os.path.join(script_dir,'CountingRate.txt'),'w');

for num in [0, 10,-10,20,-20,30,-30, 40,-40,50,-50,60,-60,-70,-80,-90]:
    f_without = open(os.path.join(script_dir,'without/' + str(num) + '.TKA'), 'r');
    f_with = open(os.path.join(script_dir,'with/' + str(num) + '.TKA'), 'r');
    L_wo = [];
    L_w = [];
    for line in f_without:
        L_wo.append(int(line.strip('\n')));
    for line in f_with:
        L_w.append(int(line.strip('\n')));

##Bin the data of binned size and compute error
    mergeSize = 50;
    #mergeWithout = [];
    #mergeWith = [];
    mergeDiff = [];
    mergeErr = [];
    for j in range(int(len(L_wo)/float(mergeSize))):
        withoutValue = sum(L_wo[j*mergeSize:(j+1)*mergeSize-1]);
        #mergeWithout.append(withoutValue/float(timeList[abs(num)]));    
        
        withValue = sum(L_w[j*mergeSize:(j+1)*mergeSize-1]);
        if num == 0:
            withoutValue = 0;
        #print withValue
        #mergeWith.append(withValue/float(timeList[abs(num)]));
        #print withValue, withoutValue,(withValue-withoutValue)/float(timeList[abs(num)])
        mergeDiff.append((withValue-withoutValue)/float(timeList[abs(num)]));
        mergeErr.append(sqrt(abs(withValue+withoutValue))/float(timeList[abs(num)]));
    if num == 0:
        print max(mergeDiff);    
    
    [truncStart, sigStart, sigEnd, truncEnd] = sigTime[num];
    index = range(len(mergeDiff));
    noiseVal = mergeDiff[truncStart:sigStart]+mergeDiff[sigEnd:truncEnd];
    noiseIdx = index[truncStart:sigStart]+index[sigEnd:truncEnd];
    noiseErr = mergeErr[truncStart:sigStart]+mergeErr[sigEnd:truncEnd];
    sigVal = mergeDiff[sigStart:sigEnd];
    sigIdx = index[sigStart:sigEnd];    
    sigErr = mergeErr[sigStart:sigEnd];
    guess = [1.0, (sigStart+sigEnd)/2,(sigEnd-sigStart)*1/4];
    #polyGuess = [-1/60, 10, -364];
     
##Fit polynomial noise
    #if num == 0:
    #    popt, pcov = curve_fit(noisePoly, noiseIdx, noiseVal, sigma = noiseErr, p0 = polyGuess);
    #else:
    popt, pcov = curve_fit(noisePoly, noiseIdx, noiseVal, sigma = noiseErr);
    [a,b,c] = popt;

##Fit signal
    popt, pcov = curve_fit(sigGaussian, sigIdx, sigVal, sigma = sigErr, p0 = guess);
    [a0, miu, sig] = popt;
    perr = np.sqrt(np.diag(pcov))
    print perr
    [a0err,miuerr, sigerr] = perr;
    

##Integration over difference
    counter = 0;
    countererr = 0;
    for i in range(sigStart,sigEnd+1):
        y1 = a*i**2+b*i+c;
        y2 = a0*np.exp(-(i-miu)**2/(2*sig**2));
        if (y2-y1) > 0:
            counter += y2-y1;
            countererr += y1+y2

##Write in file
    f.write(str(num) + ' ' + str(miu*mergeSize) + ' ' + str(counter)+  ' ' + str(miuerr)+ ' ' + str(sqrt(countererr))+ '\n');#*mergeSize) + '\n');
    counter = counter * timeList[abs(num)];
    countererr = countererr * timeList[abs(num)];
    f2.write(str(num) + ' ' +'$'+ str(int(counter))+  '\pm' + str(int(sqrt(countererr)))+'$'+  ' '+str(timeList[abs(num)]) +' '+'$'+  str(int(float(counter)/timeList[abs(num)]))+  '\pm' + "%.2f" %(sqrt(countererr)/timeList[abs(num)])+'$'+ '\n');#*mergeSize) + '\n');

##Plot figure with error bar 
    plt.clf();
    x1 = np.arange(truncStart, truncEnd, 0.2);
    plt.plot(x1, a*x1**2+b*x1+c)
    x2 = np.arange(sigStart, sigEnd, 0.2);
    plt.plot(x2, a0*np.exp(-(x2-miu)**2/(2*sig**2)))
    #plt.plot(mergeDiff,'.');
    plt.errorbar(index,mergeDiff, yerr=mergeErr,ecolor='r',fmt = 'none');
    plt.xlim([0.0,300.0])
    if num == 0:
        plt.ylim([0,150]);
    else:
        plt.ylim([-0.2,1.5]);

    plt.show();
    #plt.figure();
    plt.savefig(os.path.join(script_dir,str(num)+'.png'))

f.close();
f2.close();