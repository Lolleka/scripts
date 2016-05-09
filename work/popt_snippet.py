#!/usr/bin/python
f = open ('./popt.array','r')
f.readline()
s_params = f.readline()
a_params = s_params.split('\t')
POPT_PARAM=a_params[2]
START=int(a_params[4])
END=float(a_params[5])
NSTEPS=int(a_params[6])
STEP=float(a_params[8])
