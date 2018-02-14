#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 11:45:32 2017

@author: gudeharikishan
"""


import csv 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from peakdetect import peakdetect
from numpy import matrix
import datetime
import matplotlib.pyplot as plt2
import scipy
from scipy import stats

data = []
ifile  = open('M000005389.csv', "r")
read = csv.reader(ifile)
for row in read :
    data.append(row)

print(len(data))
arraydata = data

time=[]

for x in range(len(data)):
    time.append(data[x][0])
    
def parse_date(text):    
    for fmt in ('%m/%j/%y  %H:%M:%S %p', '%m/%j/%y  %H:%M:%S'):
        try:
            return datetime.datetime.strptime(text, fmt)
        except ValueError:
            pass
    raise ValueError('no valid date format found')

tm = []           
for x in range(len(time)):
    z = parse_date(time[x])
    tm.append(z)
tm = np.array(tm)   
 

mn = []
for x in range(len(tm)):
    v = str(tm[x])
    f = v.split(" ")
    g = f[1].split(":")
    mn.append(g[1])        
    
# identify and delete fhr and uc dates
    
sc = []
for x in range(len(mn)):
    r = np.arange(240)
    sc.append(r)

sec =[]
vbn = 0
for x in range(len(mn)):
    for y in range(240):
        sec.append(vbn+y)
    vbn =vbn+240
        
print(len(sec))
print(len(sc))
    
print(tm[1])
dt = str(tm[10])
dt = dt.split(" ")
print(dt[1])
st = str(dt[1])
st = st.split(":")
print(st[1])    


ds = arraydata

x = 0
while(x < (len(arraydata))):
    if arraydata[x][1] == 'HR2' and arraydata[x+1][1] == 'UA':
        x = x+2
        continue
    else:
        arraydata.pop(x)

x = 1        
while(x < (len(arraydata))):
    if arraydata[x][1] == 'UA' and arraydata[x-1][1] == 'HR2':
        x = x+2
        continue
    else:
        arraydata.pop(x)        
    
fhr=[]
uc=[]


for x in range(len(arraydata)):
    if arraydata[x][1] == 'HR2':
        fhr.append(arraydata[x])
    else:
        uc.append(arraydata[x])

for x in range(len(fhr)):
    fhr[x].pop(0)
    fhr[x].pop(0)
    fhr[x].pop(0)

for x in range(len(uc)):
    uc[x].pop(0)
    uc[x].pop(0) 
    uc[x].pop(0)
      
        
def flatten(l):
    flatList = []
    for elem in l:
        # if an element of a list is a list
        # iterate over this list and add elements to flatList 
        if type(elem) == list:
            for e in elem:
                flatList.append(e)
        else:
            flatList.append(elem)
    return flatList
    
fhr = flatten(fhr)
uc = flatten(uc)
list(fhr)
list(uc)
fhr = list(map(int, fhr))
uc = list(map(int, uc))
      
fhr = np.array(fhr)
uc = np.array(uc)
        
c= 1
ind = []
for x in range(len(fhr)):
    if fhr[x] < 100:
        ind.append(x)
        c= c+1
        continue      

fhr = np.delete(fhr,ind)
uc = np.delete(uc,ind)

# split windows
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]
        

fhr1 = fhr.tolist()
fhr1 = list(chunks(fhr, 4800))
fhr1 = np.array(fhr1)
fhrz = fhr1[4]
a = fhrz
uc1 = uc.tolist()
uc1 = list(chunks(uc, 4800))
uc1 = np.array(uc1)
ucz = uc1[4]
b =ucz

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
#plt.plot(yhat = savgol_filter(fhr1[0], 11, 3),color = 'red')

N = 60
fhrz = pd.rolling_mean(fhrz, N)[N-1:]

mdn = np.mean(fhrz)
mn = np.median(fhrz)
j=0

count = 0

while(abs(mdn-j) > 0.15):    
    peaks = peakdetect(fhrz, lookahead=300)
    peaks = np.array(peaks)
    maxs = np.array(peaks[0])
    mins = np.array(peaks[1])
    maxs = matrix(maxs).transpose()[0].getA()[0]
    mins = matrix(mins).transpose()[0].getA()[0]
    #print(maxs) #Maximums
    #print(mins) #Minimums
    
    ends = []
    starts = []
    
    pk = []
    m15 = mdn+10
    for k in range(len(maxs)):
        if fhrz[int(maxs[k])] > m15:
            pk.append(maxs[k])
    
    # Removing values above 170        
    pk = list(filter(lambda x : (int(fhrz[int(x)])) < 200 , pk))
    #print(pk)
    
    ec = 0
    sc = 0
    #Finding starting and ending points
    for x in range (len(pk)):
        t =int(len(fhrz)-1)
        s = int(pk[x])
        n = sc
        p = ec
        for y in range (s,t):
            if fhrz[y]-fhrz[y+1]<0:
                ends.append(y)
                ec = ec+1
                break
            if fhrz[y] == mdn:
                ends.append(y)
                ec = ec+1
                break
        for k in range (0,s-1):
            if fhrz[s-k]-fhrz[s-k-1]<0:
                starts.append(s-k)
                sc = sc +1
                break
            if fhrz[s-k] == mdn:
                starts.append(s-k)
                sc = sc+1
                break
        if sc == n:
            starts.append('NA')
        if ec == p:
            ends.append('NA')
        
        
            
            
            
    # add statement: 
        
    ev = 0
    sv = 0        
    mends = []
    mstarts = []       
    lw = []
    n15 = mdn-9     
    for k in range(len(mins)):
        if fhrz[int(mins[k])] < n15:
            lw.append(mins[k])
            
    for x in range (len(lw)):
        t =int(len(fhrz)-1)
        s = int(lw[x])
        o = sv
        r = ev
        for y in range (s,t):
            if fhrz[y]-fhrz[y+1]>0:
                mends.append(y)
                ev = ev+1
                break
            if fhrz[y] == mdn:
                mends.append(y)
                ev = ev+1
                break
        for k in range (0,s-1):
            if fhrz[s-k]-fhrz[s-k-1]>0:
                mstarts.append(s-k)
                sv = sv+1
                break
            if fhrz[s-k] == mdn:
                mstarts.append(s-k)
                sv = sv+1
                break
        if sv == o:
            mstarts.append('NA')
        if ev == r:
            mends.append('NA')
     
    np.array(starts)
    np.array(mstarts)
    np.array(ends)
    np.array(mends)
    np.array(pk)
    np.array(lw)
    
    count = 0
    check = []
    #check duration
    for x in range(len(starts)):
        if starts[x] ==  'NA':
            np.delete(starts, x)
            np.delete(ends, x)
            np.delete(pk, x)
            break
            
        if ends[x] =='NA':
            np.delete(starts, x)
            np.delete(ends, x)
            np.delete(pk, x)
            break
            
        check.append(starts[x] - ends[x])
        if ends[x] - starts[x] > 40: 
            count = count+1
        else:
            np.delete(starts, x)
            np.delete(ends, x)
            np.delete(pk, x)
            
    for x in range(len(mstarts)):
        if mstarts[x] ==  'NA':
            np.delete(mstarts, x)
            np.delete(mends, x)
            np.delete(lw, x)
            break
            
        if mends[x] =='NA':
            np.delete(mstarts, x)
            np.delete(mends, x)
            np.delete(lw, x)
            break
            
        check.append(mstarts[x] - mends[x])
        if mends[x] - mstarts[x] > 40: 
            count = count+1
        else:
            np.delete(mstarts, x)
            np.delete(mends, x)
            np.delete(lw, x)        
            
    
        
    #check value for baseline:
    index = []
    for x in range(len(starts)):
        if starts[x] ==  'NA':
            break
        if ends[x] =='NA':
            break
        else:
            for y in range(int(starts[x]),int(ends[x])):
                index.append(y)
                
    for x in range(len(mstarts)):
        if mstarts[x] ==  'NA':
            break
        if mends[x] =='NA':
            break
        else:
            for y in range(mstarts[x],mends[x]):
                index.append(y)        
        
    new_x = np.delete(fhrz, index)
    j = mdn
    mdn = np.mean(new_x)
    count = count+1
    
    


N = 50
ucz = pd.rolling_mean(ucz, N)[N-1:]
z1 = peakdetect(ucz, lookahead=300)
z = np.array(z1[0])
cont = matrix(z).transpose()[0].getA()[0]
mdn =np.median(ucz)
#plt.plot(x, y, '.-')
#plt.plot(x, ys)
#plt.show()

cends = []
cstarts = []
ec = 0
sc = 0


cn = []
for k in range(len(cont)):
    if ucz[int(cont[k])] > 20:
        cn.append(cont[k])
    
cont = np.delete(cont,cn)

#Finding starting and ending points
for x in range (len(cont)):
    t =int(len(ucz)-1)
    s = int(cont[x])
    n = sc
    p = ec
    for y in range (s,t):
        if ucz[y]-ucz[y+1]<0:
            cends.append(y)
            ec = ec+1
            break
    for k in range (0,s-1):
        if ucz[s-k]-ucz[s-k-1]<0:
            cstarts.append(s-k)
            sc = sc +1
            break
    if sc == n:
        cstarts.append('NA')
    if ec == p:
        cends.append('NA')
de = []        
for x in range(len(cont)):
    if cstarts[x] ==  'NA':
        de.append(x)
        break
    if cends[x] =='NA':
        de.append(x)
        break
cont = np.delete(cont,de)


mnz = []
for x in range(len(fhrz)):
    mnz.append(j)
    
cn = np.array(cn)
cn = cn.astype(np.int64)
lw = np.array(lw)
lw = lw.astype(np.int64)
pk = np.array(pk)
pk = pk.astype(np.int64)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)   

major_xticks = [240,480,720,960,1200,1440,1680,1920,2160,2400,2640,2880,3120,3360,3600,3840,4080,4320,4560,4800]                                              
minor_xticks = np.arange(0,4800,40)      
z = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
   
# Either one of the below snippets, processed or raw.                                    
plt.plot(ucz)
plt.plot(fhrz)
bd = plt.plot(cn,ucz[cn],'bo')
rd = plt.plot(lw,fhrz[lw],'ro')
kd = plt.plot(pk,fhrz[pk],'ko')
plt.xticks(major_xticks, z)
plt.legend([kd,rd,bd],["Accelerations", "Decelerations", "Contractions"])


#plt.plot(a)
#plt.plot(b)
#plt.xticks(major_xticks, z)

plt.grid(True)

ax.set_title('FHR and UC for 20 minutes')
ax.set_xlabel('minutes')
ax.set_xticks(major_xticks)                                                       
ax.set_xticks(minor_xticks, minor=True)                                           
ax.grid(which='minor', alpha=0.2)                                                
ax.grid(which='major', alpha=0.5)                                      

        
plt.show()                                    

    


    
