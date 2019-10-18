# egfm processing completation in python
# defination of class and method in egfm.py
from egfm import *
import matplotlib.pyplot as plt
mm = 8192 # maximum number of data points for spectra calculating
nn = 30 # maximum number of subfaults alone strike or dip direction
pi = np.pi
# ground motion data
d = np.zeros((20000, ))  # tempory array for waveform data
di = np.zeros((20000, )) # integrated array of d
dii = np.zeros((20000, )) # integrated array of di
numPoints = 15000
dm = np.zeros((numPoints, )) # acceleration data of mainshock
da = np.zeros((numPoints, )) # acceleration data of aftershock 
ds = np.zeros((numPoints, )) # acceleration data of synthetic shock
dmm = np.zeros((numPoints, )) # velocity data of mainshock 
daa = np.zeros((numPoints, )) # velocity data of aftershock
dss = np.zeros((numPoints, )) # velocity data of synthetic shock
dmmm = np.zeros((numPoints, )) # displacement data of mainshock 
daaa = np.zeros((numPoints, )) # displacement data of aftershock
dsss = np.zeros((numPoints, )) # displacement data of synthetic shock
# spacitial variables of the vector from subfault to station
rl = np.zeros((nn, nn)) # length of the vector, km
ph = np.zeros((nn, nn)) 
# ph, the angle from strike to the project of vector to x-y plane, phy
th = np.zeros((nn, nn)) # the angle from z-axis to vector, theta
gs = np.zeros((nn, nn)) # initial rupture time at subfaults
wg = np.zeros((nn, nn)) 
# factor for distance and s-wave radiation pattern correction of subfaluts
jdl = np.zeros((nn, nn))
# time lag for the waveform from each subfault 
# when they are superposed in time domian
#
spc = np.zeros((20000, ))
spcs = np.zeros((20000, ))
#
amp = np.zeros((10000, ))
fre = np.zeros((10000, ))

# read parameters 
fread = input('Enter filename of input parameters ')
# fobs = input('Enter filename for saving observed data ')
# como = input('Enter a comment for observed data ')
# fsyn = input('Enter a filename for saving synthetic data ')
# coms = input('Enter a comment for synthetic data ')
# number of padding zeros
nzero = int(0)
fileInput = open(fread, 'r')
tmp = fileInput.read()
# 61, number of parameters fread
tmp1 = tmp.split(maxsplit=61) # list of parameters, string
# kiban = int(tmp1[0])
# kkiban = int(tmp1[1])
nevnt = int(tmp1[2])
# dtl = float(tmp1[3])
ifilt = int(tmp1[4]) # index of filter
fl = float(tmp1[5]) # low frequency of filter
fh = float(tmp1[6]) # high frequency of filter
fs = float(tmp1[7]) # stop frequency
ap = float(tmp1[8]) # attenuation at pass frequency, db
astop = float(tmp1[9]) # attenuation at stop frequency, db
# ut = float(tmp1[10]) # parameters for ploting
# q = np.asfarray(tmp1[11:16])
# spatial variables of mainshock
strike = float(tmp1[17]) # strike of mainshock, degree
dip = float(tmp1[18]) # dip of mainshock
rake = float(tmp1[19]) # rake of slip for mainshock
z = float(tmp1[20]) # depth of mainshock, km
# transfer from degrees to radian
pu = np.pi / 180
strike *= pu
dip *= pu
rake *= pu
# number of componets need to synthesize
ncmpt = int(tmp1[21])
# spatial variables of aftershock
stra = float(tmp1[22])
dipa = float(tmp1[23])
raka = float(tmp1[24])
za = float(tmp1[25])
# parameters of subfault's geometry
dx = float(tmp1[26]) # length of subfaults along strike, km
dw = float(tmp1[27]) # width of subfaults along dip, km
tra = float(tmp1[28]) # rise time of after shock, second
# used in case of multiple aftershock
dx0, dw0, tr0 = np.asfarray(tmp1[29:32])
# nx and nw are number of subfaults along strike and dip
nx, nw = int(tmp1[32]), int(tmp1[33])
# number of element events summed up at each subfaults
nt = int(tmp1[34])
# in order to avoid the artificial periodicity that
# appears when NT element events are summed up with
# the interval tra
ntt = int(tmp1[35])
# nsx and nsw is the location of the start point of rupture
nsx, nsw = int(tmp1[36]), int(tmp1[37])
# NS is used in case of a multiple shock.
# CFACTOR is a parameter used to correct the difference in stress drop
# between the element event and the mainshock. 
# CFACTOR=(stress drop large event)/(stress drop small event)
ns = int(tmp1[38])
cfactor = float(tmp1[39])
# variables about propagation
vs = float(tmp1[40]) # S-Wave velocity, km/s
vr = float(tmp1[41]) # rupture velocity, km/s
ird = int(tmp1[42]) # option for style of rupture propagation
ipfm = int(tmp1[43]) # index for porpagation (radiation) pattern mechanism correction
# transfer from degree to radian
stra *= pu
dipa *= pu
raka *= pu
# filename of main and after shock
fmain = tmp1[44]
faft = tmp1[45]
# 
rm = float(tmp1[46]) # epicentral distance of mainshock, km
pm = float(tmp1[47]) # azimuth of the mainshock, degree
ra = float(tmp1[48]) # epicentral distance of aftershock, km
pa = float(tmp1[49]) # azimuth of the aftershock, degree
# transfer from degree to radian
pm *= pu
pa *= pu
icmp = int(float(tmp1[50])) # indication of the component of data
imdl = int(tmp1[51]) # the shape of slip time function
# start and end time index of main and after shocks
ksm = int(tmp1[52])
kem = int(tmp1[53])
ksa = int(tmp1[54])
kea = int(tmp1[55])
# station and magnitude 
station = tmp1[56]
magnitude = tmp1[57]
# spectral variables
# istart = int(tmp1[58]) # start data point for the spectra
# iend = int(tmp1[59])   # end data point for the spectra
# iwind = int(tmp1[60])  # window length for calculating the spectra
# close the input file
fileInput.close()
# define parameters and variables


# read main shock, get waveform data and time interval
# utilize functions in egfm.py
d, dt = readTH(fmain)
ave, d = zero(d, 0, d.size) # ave, average of data in d
di, maxVel, maxAcc = inac(dt, d.size, d)
ave, di = zero(di, 0, di.size)
dii, maxDis, maxVel = inac(dt, di.size, di)
ave, dii = zero(dii, 0, dii.size)
# bandpass acceleration, velocity and displacement of main shock
d = bandPass(d, dt, fl, fh, fs, ap, astop, nzero)
di = bandPass(di, dt, fl, fh, fs, ap, astop, nzero)
dii = bandPass(dii, dt, fl, fh, fs, ap, astop, nzero)
# stock main shock data in dm
dm = d
dmm = di
dmmm = dii
# ns, number of waveform data points for synthetic 
ns = numPoints
# reduction mean value again, after bandpass and integrate
ave, dm = zero(dm, 0, numPoints)
print('\n\n Mean of the acceleration of mainshock is {:.2f}'.format(ave), ' gal')
ave, dmm = zero(dmm, 0, numPoints)
print('Mean of the velocity of mainshock is {:.2f}'.format(ave), ' cm/sec')
ave, dmmm = zero(dmmm, 0, numPoints)
print('Mean of the velocity of mainshock is {:.2f}'.format(ave), ' cm')

# read after shock, get waveform data and time interval
# utilize functions in egfm.py
d, dt = readTH(faft)
ave, d = zero(d, 1, d.size) # ave, average of data in d
di, maxVel, maxAcc = inac(dt, d.size, d)
ave, di = zero(di, 1, di.size)
dii, maxDis, maxVel = inac(dt, di.size, di)
ave, dii = zero(dii, 1, dii.size)
# bandpass acceleration, velocity and displacement of main shock
d = bandPass(d, dt, fl, fh, fs, ap, astop, nzero)
di = bandPass(di, dt, fl, fh, fs, ap, astop, nzero)
dii = bandPass(dii, dt, fl, fh, fs, ap, astop, nzero)
# stock after shock data in da
da = d * cfactor
daa = di * cfactor
daaa = dii * cfactor
# number of points in after shock wavefor data file
na = da.size
# reduction mean value again, after bandpass and integrate
ave, da = zero(da, 0, numPoints)
print('\n\n Mean of the acceleration of aftershock is {:.2f}'.format(ave), ' gal')
ave, daa = zero(daa, 0, numPoints)
print('Mean of the velocity of aftershock is {:.2f}'.format(ave), ' cm/sec')
ave, daaa = zero(daaa, 0, numPoints)
print('Mean of the velocity of aftershock is {:.2f}'.format(ave), ' cm')


# calculation of rupture initiation time: gs(i, j) (sec)
# (x, y) coordinate of station along and normal to strike direction
# origin is the mainshock's epicenter
x, y = crdChg(rm, pm, strike)
# (r, theta, phi) of vectors from subfaults to station
rl, th, ph = rlCal(x, y, z, dx, dw, nx, nw, nsx, nsw, dip, dx0, dw0)
gs = gsCal(dx, dw, nx, nw, nsx, nsw, ird, vr)

# calculation of wg(i, j) and jdl(i, j)*dt (sec)
# wg: factor for distance and S-Wave radiation pattern correction
# jdl * dt: time lag for the waveform from each subfault used when
# the are superposed in time domine
tha = np.pi - np.arctan2(ra, za) 
# take-off angle of vector from hypocenter of aftershock to station
rr = np.sqrt(rm ** 2 + z ** 2) # hypocentral distance of mainshock
tmp = np.sqrt(ra ** 2 + za ** 2) 
# hypocentral distance of aftershock
ra = tmp
# modify factor for propagation, include distance and radiation pattern
# radiation pattern of aftershock
rdna = calRdn(ipfm, icmp, stra, dipa, raka, pa, tha)
# radiation pattern of subfaults
# suppose evch subfault has the same strike, dip and rake as major fault
rdn = calRdn(ipfm, icmp, \
        strike, dip, rake, \
        ph, th)
# rard: modify factor for radiation
if (ipfm == 1):
    rard = rdn / rdna  # rard, radiation correct factor
if (ipfm == 2):
    rard = rdn / np.abs(rdn)
if (ipfm == 0):
    rard = np.ones(ph.shape)

# really start to calculate jdl now
# wg: modify factor for distance and radiation
wg = ra * np.ones(rl.shape) / rl * rard
jdl = (gs + (rl - rr) / vs) / dt + 0.5
jdl0 = jdl.min()

# ntt: de-period factor, periodity come from same rise time
ntt = int(nt * tra / dt)
trad = tra * float(nt) / float(ntt)

# weight for superpose
weight = calWeight(wg, ntt, nt, imdl)
for i in range(nx):
    for j in range(nw):
        wgt = wg[i, j] # current wg?
        idl0 = 10000
        for k in range(ntt):
            idl = int(jdl[i, j] + (k - 1) * trad / dt + 0.5 - jdl0)
            ds = synt(ds, ns, da, na, idl, wgt, weight[i, j, k])
            dss = synt(dss, ns, daa, na, idl, wgt, weight[i, j, k])
            dsss = synt(dsss, ns, daaa, na, idl, wgt, weight[i, j, k])
        # end  k-loop
    # end j-loop
# end i-loop

ave, ds = zero(ds, 0, ns)
print('\n\n Average of synthetic acclearation is {:.2f}'.format(ave), ' gal')
ave, dss = zero(dss, 0, ns)
print('Average of synthetic velocity is {:.2f}'.format(ave), ' cm/sec')
ave, dsss = zero(dsss, 0, ns)
print('Average of synthetic displacement is {:.2f}'.format(ave), ' cm')

# calculation of correlation between observed & synthesized waveform
# and the residual and the time when correlation is maximum
nCal = 800
ca, ra, ma = cor2(dm, ds, ksm, ksa, nCal)
print('\n\n Correlation between observed and synthesized acceleration')
print('Waveforms is {:.2f}'.format(ca))
print('Residual of the acceleration waveforms is {:.2f}'.format(ra))
print('And maximum correlation time is {:^5d}'.format(ma))
cv, rv, mv = cor2(dmm, dss, ksm, ksa, nCal)
print('\n\n Correlation between observed and synthesized velocity')
print('Waveforms is {:.2f}'.format(cv))
print('Residual of the velocity waveforms is {:.2f}'.format(rv))
print('And maximum correlation time is {:^5d}'.format(mv))
cd, rd, md = cor2(dmmm, dsss, ksm, ksa, nCal)
print('\n\n Correlation between observed and synthesized displacement')
print('Waveforms is {:.2f}'.format(cd))
print('Residual of the displacement waveforms is {:.2f}'.format(rd))
print('And maximum correlation time is {:^5d}'.format(md))



