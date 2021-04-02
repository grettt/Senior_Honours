#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 16:44:38 2021

@author: gregorywilcox
"""
import numpy
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import math


data = pd.read_csv('vizier_J_A+A_638_A104_tablee1_20210131.csv', sep=",")
nump = data.to_numpy()

data2 = pd.read_csv('catalogue3.txt', sep=" ")
nump2 = data2.to_numpy()

dataJ = pd.read_csv('Sgr_model.dat', sep=",")
numpJ = dataJ.to_numpy()
print(numpJ)
   
from matplotlib.pyplot import figure

#This method calculates the Angular momentum of the data and imports the most useful data from the set before
#Tailor made for the BHB and K-Giant data supplied
def Ang_mom(array, positive):
    i = 0
    LEN=len(array)
    Ang_array = np.full((LEN, 12), 0.0)
    while (i < LEN):
        if abs(array[i, 8]) > -1:#Leftover from reproducing original Ang Mom Method
            if positive == True:# Same
                if abs(array[i, 20]) > 20:
                    Ang_array[i, 0] = array[i, 0]
                    Ang_array[i, 1] = array[i, 1]
                    Ang_array[i, 2] = array[i, 2]
                    a = array[i, 1]*array[i, 5] - array[i, 2]*array[i, 4]
                    Ang_array[i, 3] = a
                    b = array[i, 2]*array[i, 3] - array[i, 0]*array[i, 5]
                    Ang_array[i, 4] = b
                    c = array[i, 0]*array[i, 4] - array[i, 1]*array[i, 3]
                    Ang_array[i, 5] = c
                    Ang_array[i, 6] = array[i, 19]
                    Ang_array[i, 7] = array[i, 20]
                    Ang_array[i, 8] = array[i, 8] 
                    Ang_array[i, 9] = array[i, 12] 
                    Ang_array[i, 10] = array[i, 6]
                    Ang_array[i, 11] = array[i, 7] 
                   
            else:
                if abs(array[i, 20]) <= 20:
                    Ang_array[i, 0] = array[i, 0]
                    Ang_array[i, 1] = array[i, 1]
                    Ang_array[i, 2] = array[i, 2]
                    a = array[i, 1]*array[i, 5] - array[i, 2]*array[i, 4]
                    Ang_array[i, 3] = a
                    b = array[i, 2]*array[i, 3] - array[i, 0]*array[i, 5]
                    Ang_array[i, 4] = b
                    c = array[i, 0]*array[i, 4] - array[i, 1]*array[i, 3]
                    Ang_array[i, 5] = c
                    Ang_array[i, 6] = array[i, 19]
                    Ang_array[i, 7] = array[i, 20]
                    Ang_array[i, 8] = array[i, 8] 
                    Ang_array[i, 9] = array[i, 12] 
                    Ang_array[i, 10] = array[i, 6]
                    Ang_array[i, 11] = array[i, 7] 
            
        i = i+1
        
    i = 0
    j = 0    
    while (i < LEN-j):#deletes if filled values are zero (Beta >20)
        if abs(Ang_array[i, 3]) < 0.001:
            if abs(Ang_array[i, 4]) < 0.001:
                if abs(Ang_array[i, 5]) < 0.001:
                    array = np.delete(Ang_array, i, 0)
                    j = j+1
        else:
            i = i+1
            # Returns array in format [x, y, z, Lx, Ly, Lz, lambda_sgr, beta_sgr, distance, distance error, lambda, beta]
 
    return Ang_array
    
def select_sgr(array):
    
    LEN = len(array)
    i = 0
    j = 0    
    while (i < LEN-j):
        diff_Lx = 605 - array[i, 3]#slope cutoffs for selecting sagittarius cluster
        diff_Ly = -4515 - array[i, 4]
        diff_Lz = -1267 - array[i, 5]
        slope_Lxy_1 = diff_Ly + 0.21* diff_Lx #upper horizontal
        slope_Lxy_2 = -0.08* diff_Ly + diff_Lx  #right vert
        slope_Lxy_3 = 0.4*diff_Ly + diff_Lx #left vert
        slope_Lxy_4 = diff_Ly - 0.4* diff_Lx #lower horizontal 
        slope_Lxz_1 = diff_Lz - 0.3* diff_Lx #lower horiz
        slope_Lxz_2 = diff_Lz - 0.25* diff_Lx #upper horiz
        slope_Lzy_1 = diff_Ly - 2.8* diff_Lz
        slope_Lzy_2 = diff_Ly + 1.5* diff_Lz
        slope_Lzy_3 = diff_Ly - 0.1*diff_Lz
        slope_Lzy_4 = diff_Ly - 1.1*diff_Lz
        
        #distance cutoffs for selecting sgr cluster, prob not very efficient
        if slope_Lzy_1 < -4000:
            array = np.delete(array, i, 0)
            j = j+1
        else:
            if slope_Lzy_2 < -2000:
                array = np.delete(array, i, 0)
                j = j+1
            else:
                if slope_Lzy_3 > 3000:
                    array = np.delete(array, i, 0)
                    j = j+1
                else:
                    if slope_Lxy_1 < -1800:
                        array = np.delete(array, i, 0)
                        j = j+1
                    else:
                        if slope_Lxy_2 < -1600:
                            array = np.delete(array, i, 0)
                            j = j+1
                        else:
                            if slope_Lxy_3 > 2400:
                                array = np.delete(array, i, 0)
                                j = j+1
                            else:
                                if slope_Lxy_4 > 3000:
                                    array = np.delete(array, i, 0)
                                    j = j+1
                                else:
                                    if slope_Lzy_4 > 3100:
                                        array = np.delete(array, i, 0)
                                        j = j+1
                                    else:
                                        if slope_Lxz_1 > 2000:
                                            array = np.delete(array, i, 0)
                                            j = j+1
                                        else:
                                            if slope_Lxz_2 < -1300:
                                                array = np.delete(array, i, 0) #deletes if outside very specific sgr zone
                                                j = j+1
                                            else:
                                                i = i+1
        
            
    LEN = len(array)
    i = 0
    j = 0    
    while (i < LEN):
        if array[i, 6] > 180:
            array[i, 6] = array[i, 6]-360#rearranges data so it's in range -180 -> 180
        
        i = i+1
            # Returns array in same format as earlier
    
    
    return array

def OLD_select_sgr(array):
    
    LEN = len(array)#the old angular momentum method, range of 3000 kpckm/s
    i = 0
    j = 0    
    while (i < LEN-j):
        diff_Lx = 605 - array[i, 3]
        diff_Ly = -4515 - array[i, 4]
        diff_Lz = -1267 - array[i, 5]
        
        diff_ang_mom= math.sqrt(diff_Lx**2 + diff_Ly**2 + diff_Lz**2)
        
        if diff_ang_mom > 3000:
            array = np.delete(array, i, 0)
            j = j+1
        else:
            i = i+1
        
            
    LEN = len(array)
    i = 0
    j = 0    
    while (i < LEN):
        if array[i, 6] > 180:
            array[i, 6] = array[i, 6]-360#rearranges data so it's in range -180 -> 180
        
        i = i+1
    
    
    return array

def comp_array(array1, array2):#attempted comparison array
    
#this was supposed to take two arrays and compare every star with every other one and produce a plot with
#only the stars that were not on both

#it did not work
    LEN = len(array1)
    LEN2 = len(array2)
    i = 0
    j = 0  
    k = 0
    while (i < LEN-j):
        while (k < LEN2-j):
            diff_x = array1[i, 0] - array2[k, 0]
            diff_y = array1[i, 1] - array2[k, 1]
            diff_z = array1[i, 2] - array2[k, 2]
            diff_Lx = array1[i, 3] - array2[k, 3]
            diff_Ly = array1[i, 4] - array2[k, 4]
            diff_Lz = array1[i, 5] - array2[k, 5]
        
            if abs(diff_x) < 1:
                if abs(diff_y) < 1:
                    if abs(diff_z) < 1:
                        if abs(diff_Lx) < 1:
                            if abs(diff_Ly) < 1:
                                if abs(diff_Lz) < 1:
                                    array1 = np.delete(array1, i, 0)
                                    j = j+1
                                    k = 0
                                else:
                                    k = k+1
                            else:
                                k = k+1
                        else:
                            k = k+1
                    else:
                        k = k+1
                else:
                    k = k+1
            else:
                k = k+1
        i = i+1
        k = 0
                                    
        
            
    LEN = len(array1)
    i = 0
    j = 0    
    while (i < LEN):
        if array1[i, 6] > 180:
            array1[i, 6] = array1[i, 6]-360
        
        i = i+1
           
    
    
    return array1


def move_180(array, col):#simple program, moves coordinates of whatever column to between -180 and 180 if not already
    LEN = len(array)
    i = 0  
    while (i < LEN):
        if array[i, col] < -180:
            array[i, col] = array[i, col]+360
            i = i-1
        if array[i, col] > 180:
            array[i, col] = array[i, col]-360
            i = i-1
        i = i+1
            
    #prob makes some other stuff earlier redundant but oh well
    
    return array

def galac_to_sgr(array, lam, bet, d):#tried several different ways of converting coordinates that didn't seem to work
#closest I got was this
    LEN = len(array)
    i = 0  
    while (i < LEN):
        array[i, lam] = -array[i, lam] #+273.8
        array[i, bet] = -array[i, bet] #+ 13.5
        array[i, d] = array[i, d]# +25
        i = i+1
            
    return array
 
def thin_graph(array):#simple method to delete every other star from a selection
    LEN = len(array)
    i = 0  
    j = 0
    while (i < LEN-j):
        if i % 2 == 1:
            array = np.delete(array, i, 0)
            j = j+1
            i = i+1
        else:
            i = i+1
            
    
    
    return array   

def remove_2wrap(array):# method to remove second wrap from the Jorge data
    LEN = len(array)
    i = 0  
    j = 0
    while (i < LEN-j):
        if array[i, 0] == 2:
            array = np.delete(array, i, 0)
            j = j+1
            i = i+1
        else:
            i = i+1
            
    
    return array   
    
data3 = pd.read_csv('BHB_edr3_sgr.txt', sep=" ")
#print (data3)
nump3 = data3.to_numpy()


data4 = pd.read_csv('KGiant_edr3_sgr.txt', sep=" ")

nump4 = data4.to_numpy()
#read in data





B = Ang_mom(nump3, False)

K = Ang_mom(nump4, False)
#convert to ang mom


nump = move_180(nump, 16)

numpJ = galac_to_sgr(numpJ, 11, 12, 8)
numpJ = move_180(numpJ, 9)
numpJ2 = move_180(numpJ, 11)

numpJ = remove_2wrap(numpJ)
numpJ2 = remove_2wrap(numpJ2)
# a couple of attempts at working Jorge data, still not sure where I went wrong


figure(num=None, figsize=(4, 4), dpi=1280, facecolor='w', edgecolor='k')
plt.plot(B[:,3], B[:,4], 'bo', label = 'BHB stars', markersize=0.3)
plt.plot(K[:,3], K[:,4], 'ro', label = 'K-Giant stars',markersize=0.3)
plt.axis([-10000, 10000, -10000, 10000])
plt.xlabel('Angular momentum along x axis [Lx] (kpckm/s)')
plt.ylabel('Angular momentum along y axis [Ly] (kpckm/s)') 
plt.title('Angular Momentum Lx vs. Ly graph ') 
plt.legend() 
plt.savefig('Lx-Ly-uncut.png')



figure(num=None, figsize=(4, 4), dpi=1280, facecolor='w', edgecolor='k')
plt.plot(B[:,3], B[:,5], 'bo', label = 'BHB stars', markersize=0.3)
plt.plot(K[:,3], K[:,5], 'ro', label = 'K-Giant stars',markersize=0.3)
plt.axis([-10000, 10000, -10000, 10000])
plt.xlabel('Angular momentum along x axis [Lx] (kpckm/s)')
plt.ylabel('Angular momentum along z axis [Lz] (kpckm/s)') 
plt.title('Angular Momentum Lx vs. Lz graph ') 
plt.legend() 
plt.savefig('Lx-Lz-uncut.png')


figure(num=None, figsize=(4, 4), dpi=1280, facecolor='w', edgecolor='k')
plt.plot(B[:,5], B[:,4], 'bo', label = 'BHB stars', markersize=0.3)
plt.plot(K[:,5], K[:,4], 'ro', label = 'K-Giant stars',markersize=0.3)
plt.axis([-10000, 10000, -10000, 10000])
plt.xlabel('Angular momentum along z axis [Lz] (kpckm/s)')
plt.ylabel('Angular momentum along y axis [Ly] (kpckm/s)') 
plt.title('Angular Momentum Lz vs. Ly graph ') 
plt.legend() 
plt.savefig('Lz-Ly-uncut.png')


B2 = select_sgr(B)
K2 = select_sgr(K)
#select sgr stars


B2 = move_180(B2, 10)
K2 = move_180(K2, 10)



nump = move_180(nump, 10)

B3 = OLD_select_sgr(B)
K3 = OLD_select_sgr(K)
#old graph for comparison


COMP_B = comp_array(B2, B3)
COMP_K = comp_array(K2, K3)
COMP_B2 = comp_array(B3, COMP_B)
COMP_K2 = comp_array(K3, COMP_K)
COMP_B2 = comp_array(COMP_B2, B2)
COMP_K2 = comp_array(COMP_K2, K2)
#attempts to compare star selection, big failure

print('new sgr')
print(len(K2)+len(B2))
print('old sgr')
print(len(K3)+len(B3))
#gives the number of stars in each


print('')

print('comp sgr')
print(len(COMP_B2)+len(COMP_K2))
#give compared stars-ignore


#producing lots of graphs, lots of comments left over from copying the plot code
figure(num=None, figsize=(4, 4), dpi=1280, facecolor='w', edgecolor='k')
plt.plot(B2[:,3], B2[:,4], 'bv', label = 'BHB stars', markersize=0.3)
plt.plot(K2[:,3], K2[:,4], 'rs', label = 'K-Giant stars',markersize=0.3)
plt.axis([-10000, 10000, -10000, 10000])
plt.xlabel('Angular momentum along x axis [Lx] (kpckm/s)')
plt.ylabel('Angular momentum along y axis [Ly] (kpckm/s)') 
plt.title('Angular Momentum Lx vs. Ly graph ') 
plt.legend() 
plt.savefig('Lx-Ly.png')


figure(num=None, figsize=(4, 4), dpi=1280, facecolor='w', edgecolor='k')
plt.plot(B2[:,3], B2[:,5], 'bv', label = 'BHB stars', markersize=0.3)
plt.plot(K2[:,3], K2[:,5], 'rs', label = 'K-Giant stars',markersize=0.3)
plt.axis([-10000, 10000, -10000, 10000])
plt.xlabel('Angular momentum along x axis [Lx] (kpckm/s)')
plt.ylabel('Angular momentum along z axis [Lz] (kpckm/s)') 
plt.title('Angular Momentum Lx vs. Lz graph ') 
plt.legend() 
plt.savefig('Lx-Lz.png')


figure(num=None, figsize=(4, 4), dpi=1280, facecolor='w', edgecolor='k')
plt.plot(B2[:,5], B2[:,4], 'bv', label = 'BHB stars', markersize=0.3)
plt.plot(K2[:,5], K2[:,4], 'rs', label = 'K-Giant stars',markersize=0.3)
plt.axis([-10000, 10000, -10000, 10000])
plt.xlabel('Angular momentum along z axis [Lz] (kpckm/s)')
plt.ylabel('Angular momentum along y axis [Ly] (kpckm/s)') 
plt.title('Angular Momentum Lz vs. Ly graph ') 
plt.legend() 
plt.savefig('Lz-Ly.png')

figure(num=None, figsize=(4, 4), dpi=1280, facecolor='w', edgecolor='k')
plt.plot(B3[:,3], B3[:,4], 'bv', label = 'BHB stars', markersize=0.3)
plt.plot(K3[:,3], K3[:,4], 'rs', label = 'K-Giant stars',markersize=0.3)
plt.axis([-10000, 10000, -10000, 10000])
plt.xlabel('Angular momentum along x axis [Lx] (kpckm/s)')
plt.ylabel('Angular momentum along y axis [Ly] (kpckm/s)') 
plt.title('Angular Momentum Lx vs. Ly graph ') 
plt.legend() 
plt.savefig('OLDLx-Ly.png')


figure(num=None, figsize=(4, 4), dpi=1280, facecolor='w', edgecolor='k')
plt.plot(B3[:,3], B3[:,5], 'bv', label = 'BHB stars', markersize=0.3)
plt.plot(K3[:,3], K3[:,5], 'rs', label = 'K-Giant stars',markersize=0.3)
plt.axis([-10000, 10000, -10000, 10000])
plt.xlabel('Angular momentum along x axis [Lx] (kpckm/s)')
plt.ylabel('Angular momentum along z axis [Lz] (kpckm/s)') 
plt.title('Angular Momentum Lx vs. Lz graph ') 
plt.legend() 
plt.savefig('OLDLx-Lz.png')


figure(num=None, figsize=(4, 4), dpi=1280, facecolor='w', edgecolor='k')
plt.plot(B3[:,5], B3[:,4], 'bv', label = 'BHB stars', markersize=0.3)
plt.plot(K3[:,5], K3[:,4], 'rs', label = 'K-Giant stars',markersize=0.3)
plt.axis([-10000, 10000, -10000, 10000])
plt.xlabel('Angular momentum along z axis [Lz] (kpckm/s)')
plt.ylabel('Angular momentum along y axis [Ly] (kpckm/s)') 
plt.title('Angular Momentum Lz vs. Ly graph ') 
plt.legend() 
plt.savefig('OLDLz-Ly.png')

figure(num=None, figsize=(8, 6), dpi=1280, facecolor='w', edgecolor='k')
plt.plot( nump2[:,21], nump2[:,14], 'g.', label = '5D Vasiliev data', markersize=0.4)
plt.plot( nump[:,16], nump[:,2], 'g.', label = '5D Ramos data', markersize=0.4)
plt.plot(K2[:,6], K2[:,8], 'rs', label = '6D K-Giant Stars', markersize=0.4)
plt.plot(B2[:,6], B2[:,8], 'bv', label = '6D BHB Stars', markersize=0.4)
#plt.axis([-10000, 10000, -10000, 10000])
#nump2[:,21], nump2[:,14], 'bo',
#plt.show()
plt.xlabel('Galactocentric star coordinates in the Sagittarius frame [Lambda_Sgr] (degrees)')
plt.ylabel('Distance [D] (kpc)') 
plt.title('Comparison of 5D Sagittarius model and 6D stars selected using angular momentum method') 
plt.legend() 
plt.savefig('PRODUCED_GRAPH.png')


figure(num=None, figsize=(8, 6), dpi=1280, facecolor='w', edgecolor='k')
plt.plot( nump2[:,21], nump2[:,22], 'g.', label = '5D Vasiliev data', markersize=0.4)
#plt.plot( nump[:,16], nump[:,17], 'g*', label = '5D Ramos data - reduced for clarity', markersize=0.4)
plt.plot(K2[:,6], K2[:,7], 'rv', label = '6D K-Giant Stars', markersize=0.4)
plt.plot(B2[:,6], B2[:,7], 'bs', label = '6D BHB Stars', markersize=0.4)
#plt.axis([-10000, 10000, -10000, 10000])
#nump2[:,21], nump2[:,14], 'bo',
#plt.show()
plt.xlabel('Galactocentric star coordinates in the Sagittarius frame [Lambda_Sgr] (degrees)')
plt.ylabel('Distance [D] (kpc)') 
plt.title('Comparison of 5D Sagittarius model and 6D stars selected using angular momentum method') 
plt.legend() 
plt.savefig('PRODUCED_GRAPH_reduced_lambet.png')

figure(num=None, figsize=(8, 6), dpi=1280, facecolor='w', edgecolor='k')
plt.plot( nump[:,16], nump[:,2], 'b.', label = '5D Ramos data Model', markersize=0.4)
plt.plot( nump2[:,21], nump2[:,14], 'b.', label = '5D Vasiliev Model', markersize=0.4)
plt.plot( numpJ[:,11], numpJ[:,8], 'g.', label = '6D Peñarrubia Model', markersize=0.4)
plt.plot(K2[:,6], K2[:,8], 'rv', label = '6D K-Giant Stars', markersize=0.4)
plt.plot(B2[:,6], B2[:,8], 'rs', label = '6D BHB Stars', markersize=0.4)
plt.axis([-200, 200, -5, 100])
#nump2[:,21], nump2[:,14], 'bo',
#plt.show()
plt.xlabel('Galactocentric star coordinates in the Sagittarius frame [Lambda_Sgr] (degrees)')
plt.ylabel('Distance [D] (kpc)') 
plt.title('Comparison of 5D Sagittarius model and 6D stars selected using angular momentum method') 
plt.legend() 
plt.savefig('PRODUCED_GRAPH_Jorge4.png')

figure(num=None, figsize=(8, 6), dpi=1280, facecolor='w', edgecolor='k')
plt.plot( nump[:,16], nump[:,2], 'b.', label = '5D Ramos data', markersize=0.4)
plt.plot( nump2[:,21], nump2[:,14], 'b.', label = '5D Vasiliev data', markersize=0.4)
#plt.plot( numpJ[:,9], numpJ[:,8], 'go', label = '6D Peñarrubia Model', markersize=0.4)
plt.plot(K2[:,6], K2[:,8], 'rv', label = '6D K-Giant Stars', markersize=0.4)
plt.plot(B2[:,6], B2[:,8], 'rs', label = '6D BHB Stars', markersize=0.4)
plt.axis([-200, 200, -5, 100])
#nump2[:,21], nump2[:,14], 'bo',
#plt.show()
plt.xlabel('Galactocentric star coordinates in the Sagittarius frame [Lambda_Sgr] (degrees)')
plt.ylabel('Distance [D] (kpc)') 
plt.title('Comparison of 5D Sagittarius model and 6D stars selected using angular momentum method') 
plt.legend() 
plt.savefig('PRODUCED_GRAPH_Jorge3.png')


#thinning the graphs for visibility
nump = thin_graph(nump)
nump = thin_graph(nump)
nump = thin_graph(nump)
nump = thin_graph(nump)

nump2 = thin_graph(nump2)
nump2 = thin_graph(nump2)
nump2 = thin_graph(nump2)

figure(num=None, figsize=(8, 6), dpi=1280, facecolor='w', edgecolor='k')
plt.plot( nump2[:,21], nump2[:,14], 'g.', label = '5D Vasiliev Model', markersize=0.4)
plt.plot(K2[:,6], K2[:,8], 'rv', label = '6D K-Giant Stars', markersize=0.4)
plt.plot(B2[:,6], B2[:,8], 'bs', label = '6D BHB Stars', markersize=0.4)
#plt.axis([-10000, 10000, -10000, 10000])
#nump2[:,21], nump2[:,14], 'bo',
#plt.show()
plt.xlabel('Galactocentric star coordinates in the Sagittarius frame [Lambda_Sgr] (degrees)')
plt.ylabel('Distance [D] (kpc)') 
plt.title('Comparison of 5D Sagittarius model and 6D stars selected using angular momentum method') 
plt.legend() 
plt.savefig('PRODUCED_GRAPH_reduced.png')

figure(num=None, figsize=(8, 6), dpi=1280, facecolor='w', edgecolor='k')
plt.plot( numpJ[:,9], numpJ[:,8], 'g.', label = '6D Peñarrubia Model', markersize=0.4)
plt.plot(K2[:,10], K2[:,8], 'rv', label = '6D K-Giant Stars', markersize=0.4)
plt.plot(B2[:,10], B2[:,8], 'bs', label = '6D BHB Stars', markersize=0.4)
#plt.axis([-10000, 10000, -10000, 10000])
#nump2[:,21], nump2[:,14], 'bo',
#plt.show()
plt.axis([-200, 200, -5, 100])
plt.xlabel('Galactocentric star coordinates in longitude (degrees)')
plt.ylabel('Distance [D] (kpc)') 
plt.title('Comparison of 5D Sagittarius model and 6D stars selected using angular momentum method') 
plt.legend() 
plt.savefig('PRODUCED_GRAPH_Jorge.png')

figure(num=None, figsize=(8, 6), dpi=1280, facecolor='w', edgecolor='k')
plt.plot( numpJ[:,9], numpJ[:,8], 'g.', label = '6D Peñarrubia Model', markersize=0.4)
plt.plot( nump[:,10], nump[:,2], 'b.', label = '5D Ramos Model - reduced ', markersize=0.4)
plt.plot(K2[:,10], K2[:,8], 'rv', label = '6D K-Giant Stars', markersize=0.4)
plt.plot(B2[:,10], B2[:,8], 'rs', label = '6D BHB Stars', markersize=0.4)
#plt.axis([-10000, 10000, -10000, 10000])
#nump2[:,21], nump2[:,14], 'bo',
#plt.show()
plt.axis([-200, 200, -5, 100])
plt.xlabel('Galactocentric star coordinates in longitude (degrees)')
plt.ylabel('Distance [D] (kpc)') 
plt.title('Comparison of 5D Sagittarius model and 6D stars selected using angular momentum method') 
plt.legend() 
plt.savefig('PRODUCED_GRAPH_Jorge_EVR.png')


figure(num=None, figsize=(8, 6), dpi=1280, facecolor='w', edgecolor='k')
plt.plot( nump[:,16], nump[:,2], 'g*', label = '5D Ramos Model - reduced for clarity', markersize=0.4)
plt.plot(K3[:,6], K3[:,8], 'rv', label = '6D K-Giant Stars', markersize=0.4)
plt.errorbar(K3[:,6], K3[:,8], yerr=K3[:,9], fmt='.k', ecolor='red', elinewidth=0.1, capsize=0, markersize=0);
plt.plot(B3[:,6], B3[:,8], 'bs', label = '6D BHB Stars', markersize=0.4)
plt.errorbar(B3[:,6], B3[:,8], yerr=B3[:,9], fmt='.k', ecolor='blue', elinewidth=0.1, capsize=0, markersize=0);
#plt.axis([-10000, 10000, -10000, 10000])
#nump2[:,21], nump2[:,14], 'bo',
#plt.show()
plt.xlabel('Galactocentric star coordinates in the Sagittarius frame [Lambda_Sgr] (degrees)')
plt.ylabel('Distance [D] (kpc)') 
plt.title('Comparison of 5D Sagittarius model and 6D stars with error bars') 
plt.legend() 
plt.savefig('PRODUCED_GRAPH_reduced3.png')


figure(num=None, figsize=(8, 6), dpi=1280, facecolor='w', edgecolor='k')
plt.plot( nump[:,16], nump[:,2], 'g*', label = '5D Ramos Model - reduced for clarity', markersize=0.4)
plt.plot(K2[:,6], K2[:,8], 'rv', label = '6D K-Giant Stars', markersize=0.4)
plt.errorbar(K2[:,6], K2[:,8], yerr=K2[:,9], fmt='.k', ecolor='red', elinewidth=0.1, capsize=0, markersize=0);
plt.plot(B2[:,6], B2[:,8], 'bs', label = '6D BHB Stars', markersize=0.4)
plt.errorbar(B2[:,6], B2[:,8], yerr=B2[:,9], fmt='.k', ecolor='blue', elinewidth=0.1, capsize=0, markersize=0);
#plt.axis([-10000, 10000, -10000, 10000])
#nump2[:,21], nump2[:,14], 'bo',
#plt.show()
plt.xlabel('Galactocentric star coordinates in the Sagittarius frame [Lambda_Sgr] (degrees)')
plt.ylabel('Distance [D] (kpc)') 
plt.title('Comparison of 5D Sagittarius model and 6D stars with error bars') 
plt.legend() 
plt.savefig('PRODUCED_GRAPH_reduced2.png')

figure(num=None, figsize=(8, 6), dpi=1280, facecolor='w', edgecolor='k')
plt.plot(K2[:,6], K2[:,8], 'rv', B2[:,6], B2[:,8], 'bs', markersize=0.4)
#plt.axis([-10000, 10000, -10000, 10000])
#plt.show()
plt.xlabel('Galactocentric star coordinates in the Sagittarius frame [Lambda_Sgr] (degrees)')
plt.ylabel('Distance [D] (kpc)') 
plt.title('6D stars selected using the new angular momentum method') 
plt.savefig('PRODUCED_GRAPH_no_blue.png')

figure(num=None, figsize=(8, 6), dpi=1280, facecolor='w', edgecolor='k')
plt.plot(K3[:,6], K3[:,8], 'rv', B3[:,6], B3[:,8], 'bs', markersize=0.4)
#plt.axis([-10000, 10000, -10000, 10000])
#plt.show()
plt.xlabel('Galactocentric star coordinates in the Sagittarius frame [Lambda_Sgr] (degrees)')
plt.ylabel('Distance [D] (kpc)') 
plt.title('6D stars selected using the original angular momentum method') 
plt.savefig('PRODUCED_GRAPH_no_blue_2.png')

figure(num=None, figsize=(8, 6), dpi=1280, facecolor='w', edgecolor='k')
plt.plot( nump2[:,21], nump2[:,22], 'go', label = '5D Vasiliev Model', markersize=0.4)
plt.plot( nump[:,16], nump[:,17], 'g*', label = '5D Ramos Model - reduced for clarity', markersize=0.4)
plt.plot(K2[:,6], K2[:,7], 'rv', label = '6D K-Giant Stars', markersize=0.4)
plt.plot(B2[:,6], B2[:,7], 'bs', label = '6D BHB Stars', markersize=0.4)
#plt.axis([-10000, 10000, -10000, 10000])
#nump2[:,21], nump2[:,14], 'bo',
#plt.show()
plt.xlabel('Galactocentric star coordinates in the Sagittarius frame [Lambda_Sgr] (degrees)')
plt.ylabel('Galactocentric star coordinates in the Sagittarius frame [Beta_Sgr] (degrees)') 
plt.title('Comparison of 5D Sagittarius model and 6D stars in the Lambda-Beta frame') 
plt.legend() 
plt.savefig('PRODUCED_GRAPH_reduced_lambet.png')

figure(num=None, figsize=(8, 6), dpi=1280, facecolor='w', edgecolor='k')
plt.plot(COMP_K[:,3], COMP_K[:,4], 'rv', COMP_B[:,3], COMP_B[:,4], 'bs', markersize=0.4)
plt.axis([-10000, 10000, -10000, 10000])
#plt.show()
plt.xlabel('Angular momentum along x axis [Lx] (kpckm/s)')
plt.ylabel('Angular momentum along y axis [Ly] (kpckm/s)') 
plt.title('comparison of two methods') 
plt.savefig('PRODUCED_GRAPH_compar.png')

figure(num=None, figsize=(8, 6), dpi=1280, facecolor='w', edgecolor='k')
plt.plot(COMP_K2[:,3], COMP_K2[:,4], 'rv', COMP_B2[:,3], COMP_B2[:,4], 'bs', markersize=0.4)
plt.axis([-10000, 10000, -10000, 10000])
#plt.show()
plt.xlabel('Angular momentum along x axis [Lx] (kpckm/s)')
plt.ylabel('Angular momentum along y axis [Ly] (kpckm/s)') 
plt.title('comparison of two methods') 
plt.savefig('PRODUCED_GRAPH_compar2.png')
