"""
I'm going to make a training set with a lot of points, then split it somehow, and then train individual emulators and then combine their results, comparing this with a single emulator.

To sean: apologies ahead of time for not commenting this. I made it in, like, an hour.
"""

import numpy as np
import os, sys
import emulator
import matplotlib.pyplot as plt
np.random.seed(85719)

#Basic training data
N = 120 #Multiple of 2,3,4,5,6
x = np.linspace(0,10,num=N)
yerr = np.fabs(0.1 + 0.5*np.random.rand(N))
y = np.cos(x) + 0.5*np.random.randn(N) + 1

#Create a SINGLE emulator
true_emu = emulator.Emulator(name="truth",xdata=x, ydata=y, yerr=yerr)
true_emu.train()

#Predict
xstar = np.linspace(min(x)-1,max(x)+1,N*2)
ystar, ystarvar = true_emu.predict(xstar)
ystarerr = np.sqrt(ystarvar)

def plot_emu_lines(x,y,yerr,c,errs=True):
    plt.plot(x,y,ls='-',c=c)
    if errs:
        plt.plot(x,y+yerr,ls="--",c=c)
        plt.plot(x,y-yerr,ls="--",c=c)

plt.rc('text',usetex=True, fontsize=20)
plt.errorbar(x,y,np.fabs(yerr),ls='',marker='o',color='k',ms=8,label="f")
plot_emu_lines(xstar,ystar,ystarerr,c='r')
#Labels
plt.xlabel(r"$x$",fontsize=24)
plt.ylabel(r"$y$",fontsize=24)
plt.subplots_adjust(bottom=0.15,left=0.15)
#plt.show()
plt.clf()

#Loop over
cs = ["r","g","b","y","purple","cyan","pink"]
for i in xrange(2,4): #i is the number of emulators
    plot_emu_lines(xstar,ystar,ystarerr,c='k')
    stepsize = N/i #integer division
    emu_list = []
    #Predict with all of them and figure out
    ystar_list = []
    ystarvar_list = []
    ystarerr_list = []
    xj_list = []
    yj_list = []
    yerrj_list = []
    for j in xrange(0,i): #Loop over emulators
        xj = x[j*stepsize:(j+1)*stepsize]
        yj = y[j*stepsize:(j+1)*stepsize]
        xj_list.append(xj)
        yj_list.append(yj)
        yerrj = yerr[j*stepsize:(j+1)*stepsize]
        yerrj_list.append(yerrj)
        jemu = emulator.Emulator(name="i%d_j%d"%(i,j),xdata=xj,ydata=yj,yerr=yerrj)
        jemu.train()
        emu_list.append(jemu)
        ystarj, ystarvarj = emu_list[j].predict(xstar)
        ystarerrj = np.sqrt(ystarvarj)
        ystar_list.append(ystarj)
        ystarvar_list.append(ystarvarj)
        ystarerr_list.append(ystarerrj)

    for j in xrange(0,i): #Loop over emulators and plot
        plt.errorbar(xj_list[j],yj_list[j],yerrj_list[j],marker='o',ls='',c=cs[j],alpha=0.25)
        plot_emu_lines(xstar,ystar_list[j],ystarerr_list[j],c=cs[j],errs=False)
    plt.show()
    plot_emu_lines(xstar,ystar,ystarerr,c='k')
    for j in xrange(0,i): #Loop over emulators and plot
        plt.errorbar(xj_list[j],yj_list[j],yerrj_list[j],marker='o',ls='',c=cs[j],alpha=0.25)
        plot_emu_lines(xstar,ystar_list[j],ystarerr_list[j],c=cs[j])
    plt.show()

    #Do the combining
    ystar_combined = np.zeros_like(xstar)
    ystarvar_combined = np.zeros_like(xstar)
    for j in xrange(0,i):
        ystarvar_combined += ystarerr_list[j]**-2
    ystarvar_combined = ystarvar_combined**-1
    for j in xrange(0,i):
        ystar_combined += ystarvar_combined * ystar_list[j]/ystarerr_list[j]**2
    plot_emu_lines(xstar,ystar,ystarerr,c='k')
    plot_emu_lines(xstar,ystar_combined,np.sqrt(ystarvar_combined),c='r')
    plt.show()
