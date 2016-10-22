#!/usr/bin/python

import numpy, sys

inf=open(sys.argv[1],'r')

lines=inf.readlines()

n=len(lines)

r=numpy.empty((n))
v=numpy.empty((n))
f=numpy.empty((n))

for i in range(n):
    data=lines[i].split()
    r[i]=float(data[1])
    v[i]=float(data[2])
    f[i]=float(data[3])

for i in range(1,n-1):

    dr=r[i+1]-r[i-1]
    dv=v[i+1]-v[i-1]
    fnum=-dv/dr

    print r[i],f[i],fnum


