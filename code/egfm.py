# Empirical Green's Function Modify
# complete in python
import numpy as np
from scipy import signal

def readTH(fn):
    '''
    read time history of acceleration file
    get data and timeinterval
    '''
    numData = 0
    with open(fn, 'r') as f:
        while (f.readline() != ''):
            numData += 1
        f.seek(0)
        tmp = f.read()
        ms = numData * 2
        data = np.zeros((numData, ))
        tmp1 = tmp.split(maxsplit=ms)
        tmp2 = np.asfarray(tmp1)
        dt = tmp2[0]
        for i in range(int(tmp2.size / 2)):
            data[i] = tmp2[2 * i + 1]
    return data, dt

def zero(data, ns, ne):
  '''
  sub-functions like command rmean in SAC
  ns: index of start point
  ne: index of end point
  ave: average of data
  '''
  ave = np.mean(data[ns:ne])
  data -= ave
  return ave, data

def inac(dt, ndata, acc):
    '''
    INtegration of ACceleration
    dt: time interval
    ndata: number of data points
    acc: acceleration (gal)
    vel: velocity (cm/sec)
    maxVel: maximum absolute velocity
    maxAcc: maximum absolute acceleration
    '''
    vel = np.zeros(acc.shape)
    maxVel = 0
    maxAcc = 0
    vel[0] = acc[0] * dt
    for i in range(1, ndata):
        vel[i] = vel[i - 1] + (acc[i - 1] + acc[i]) / 2 * dt
        maxVel = np.maximum(maxVel, np.abs(vel[i]))
        maxAcc = np.maximum(maxAcc, np.abs(acc[i]))
    return vel, maxVel, maxAcc




def crdChg(r, theta, strike):
  '''
  coordinate change function (r, theta) --> (x, y) 
  r  : epicentral distance                                          
  theta:  azimth from source to receiver (radian) clockwise            
  strike   (radian)                                            
  x,y :  new corrdinate of station, km
  x: along strike
  y: normal to strike
  epicenter as origin point in X-Y plane
  '''
  x = r * np.cos(theta - strike)
  y = r * np.sin(theta - strike)
  return x, y

#     note: z-axis down to the earth
#     theta: angle between vector and z-axis (radian)
#     phi: angle between vector and strike (radian)
def rlCal(x, y, z, dx, dw, nl, nw, nsx, nsw, dip, dx0, dw0):
  '''
  function for calculation of (r, theta, phi) for the vector         
  from each subfault to station                                     
  output : r(i,j), theta(i,j), phi(i,j)
  '''
  # initial value of matrix r, theta, phi
  r = np.zeros((nl, nw))
  theta = np.zeros((nl, nw)) 
  phi = np.zeros((nl, nw))
  # nl, nw, number of subfaults along strike and dip
  for l in range(nl):
    for m in range(nw):
      # nsx, nsw, index of rupture start point along strike and dip
      # dx, dw, length of subfault along strike and dip
      # dx0, dw0, multiple asperities
      dxl = (l + 1 - nsx) * dx - dx0
      dwm = (m + 1 - nsw) * dw - dw0
      # dxl, dwm, subfault center point coordinate in fault surface
      # (dxl, dwm), epicenter (0, 0)
      # general coordinate systme: North-East-Down
      # change to Strike-Normal-Down system is
      # subfault center: (dxl, dwm * cos(dip), z + dwm * sin(dip))
      # station: (x, y, 0)
      # Then, vector from subfault to station is:
      # (x - dxl, y - dwm * cos(dip), -(z + dwm * sin(dip)))
      xx = x - dxl
      yy = y - dwm * np.cos(dip)
      zz = 0.0 - (z + dwm * np.sin(dip))
      rr = np.sqrt(xx ** 2 + yy ** 2)
      r[l][m] = np.sqrt(xx ** 2 + yy ** 2 + zz ** 2)
      phi[l][m] = np.arctan2(yy, xx)
      theta[l][m] = np.arctan2(rr, zz)
  return r, theta, phi

def gsCal(dx, dw, nx, nw, nsx, nsw, ird, vr):
  '''
  function for calculation of rupture initiation time
  Gisin Start Calculation !!!
  vr: rupture velocity on fault surface
  ird: index of rupture radiation
  1: along strike direction
  2: circle
  '''
  gs = np.zeros((nx, nw))
  for i in range(nx):
    for j in range(nw):
      if (ird == 1):
          gs[i][j] = np.abs(i + 1 - nsx) * dx / vr
      elif (ird == 2):
          x = np.abs(i + 1 - nsx) * dx
          y = np.abs(j + 1 - nsw) * dw
          gs[i][j] = np.sqrt(x ** 2 + y ** 2) / vr
      else :
          raise ValueError("Wrong index of rupture radiation")
  return gs

def calWeight(mfp, ntt, nt, imdl):
    '''
    calculate weight for superpose in time domine
    mfp: modify factor matrix for subfaults
    ntt: number of times need to stack for de-periodity
    nt: number of elements in each subfaults
    imdl: index of modify 
    weight: weight matrix for supperpose
    '''
    tmpTrans = np.ones((ntt, ))
    if (imdl == 1):
        for k in range(ntt):
            tmpTrans[k] = float(nt) / float(ntt)
        # end for
    elif (imdl == 2):
        for k in range(1, ntt):
            tmpTrans[k] = float(nt - 1) / float(ntt)
    elif (imdl == 3):
        for k in range(1, ntt):
            tmpTrans[k] = float(nt - 1) / float(ntt) * \
                    np.exp(1 + (1 - k) / float(ntt)) / (np.exp(1) - 1)
    elif (imdl == 4):
        for k in range(1, ntt):
            tmpTrans[k] = float(nt - 1) / (2 * (ntt * (k - 1)) ** 0.5)
    # end if
    tmpTrans = tmpTrans.reshape((1, ntt))
    weight = np.expand_dims(mfp, 2) * tmpTrans
    return weight




def synt(dsyn, nsyn, da, na, ndelay, wg, wei):
    '''
    function for waveform synthesis, time domine superposition
    output: dsyn
    ndelay: number of data points to delay
    da: data of after shock, element earthquake
    na: number of after shock data points
    wg: weight of subfalut
    wei: weight of element waveform in subfault 
    '''
    npractical = na
    if ((ndelay + na) >= nsyn): 
        npractical = nsyn - ndelay
    dsyn[ndelay:(ndelay + npractical - 1)] += wg * wei * da[0:npractical - 1]
    return dsyn

def calRdn(ipfm, icmp, strike, dip, rake, azimuth, theta):
    '''
    calculate modify factor for propagation
    include distance and radiation pattern
    ipfm: Index of Propagation Factor Modify
    fault surface spatial variables: strike, dip
    slip angle: rake
    azimuth: angle from north to vecotr between epicenter and station
    theta: take-off angle of vector between hypocenter and station
    rdna: radiation of aftershock
    '''
    if (ipfm >= 1):
        # radiation pattern of SH and SV
        rdsh, rdsv = rdatn(strike, dip, rake, azimuth, theta)
        if (icmp == 1000):
            rdn = rdsv * np.sin(theta)
        elif (icmp == -1000):
            rdn = -rdsv * np.sin(theta)
        else :
            rdn = rdsv * np.cos(theta) * np.cos(icmp - azimuth) + \
                rdsh * np.sin(icmp - azimuth)
    else :
        rdn = np.ones(theta.shape)
    return rdn

def calJdl(strike, dip, rake, ra, rl, theta, phi, gs):
    '''
    calculate maxtrix of time lag for subfaults: jdl
    major fault's spatial variables: strike, dip, rake
    subfault's spatial variables: rl, theta, phi
    gs: start time of rupture for subfaults
    ra: hypocenter distance of aftershock
    '''
    


    

def rdatn(strike, dip, rake, az, theta):
    '''
    calculation of the radiation term induced by the direct body wave 
    strike  (radian)                                            
    dip : dip angle                                                   
    rake : rake angle                                                  
    az  : azimth from sorce to receiver (clockwise from north)        
    theta  : take-off angle (from down)                                  
    rdSh: radiation for sh-wave (output)                              
    rdSv: radiation for sv-wave (output)                              
    '''
    sr = np.sin(rake); cr = np.cos(rake)
    sd = np.sin(dip); cd = np.cos(dip)
    st = np.sin(theta); ct = np.cos(theta)
    ss = np.sin(az - strike); cs = np.cos(az - strike)

    rdSv = sr * (cd ** 2 - sd ** 2) * (ct ** 2 - st ** 2) * ss - \
    cr * cd * (ct ** 2 - st ** 2) * cs + \
    cr * sd *st * ct * 2 * ss * cs - \
    sr * sd * cd * 2 * st * ct * (1 + ss ** 2)
    rdSh = cr * cd *ct * ss + cr * sd * (cs ** 2 - ss ** 2) + \
    sr * (cd ** 2 - sd ** 2) * ct * cs - \
    sr * sd * cd * st * 2 * ss *cs
    return rdSh, rdSv

def cor2(dMain, dSyn, ksm, ksa, nCal):
    '''
    function for calculation of residual and correlation
    by move average method
    dMain: data of mainshock
    dSyn: data of synthetic target event
    ksm: sequiral number of start points of mainshock
    ksa: sequiral number of start points of aftershock
    nCal: number of points used to calculate residual and correlation
    maxTime: when the correlation is maximum
    '''
    nsum = 10000
    maxTime = 0
    dm = dMain
    ds = dSyn
    # pad zeros after data
    if (dMain.size < nsum):
        dm = np.concatenate((dMain, np.zeros((nsum - dMain.size))))
    if (dSyn.size < nsum):
        ds = np.concatenate((dSyn, np.zeros((nsum - dSyn.size))))
    nt = 200 # 1/50th of length of dMain
    c = 0.0 # initial correlationon
    r = 0.0 # initial residual
    # calculate correlation and maxTime
    for ntime in range(-nt, nt + 1):
        x, y, z = 0.0, 0.0, 0.0
        for i in range(1, nCal + 1):
            x += dm[ksm - 1 + i] * ds[ksa - 1 + i + ntime]
            y += dm[ksm - 1 + i] ** 2
            z += ds[ksa - 1 + i + ntime] ** 2
            # end for loop
        co = x / (np.sqrt(y) * np.sqrt(z))
        if (co > c):
            c = co
            maxTime = ntime
            # end if
    # calculate residual
    s, p, q = 0.0, 0.0, 0.0
    for i in range(1, nCal + 1):
        p += dm[ksm - 1 + i] ** 2
        q += ds[ksa - 1 + i] ** 2
        s += (dm[ksm - 1 + i] - ds[ksa - 1 + i + maxTime]) ** 2
        # end for loop
    r = s / (np.sqrt(p) * np.sqrt(q))
    return c, r, maxTime

#------------------filter section--------------------------------
def recFil(x, y, n, h, nml):
    '''
    recursive filering, frequency domain: 
    f(z) = (1 + a * z + aa * z ** 2) / (1 + b * z + bb * z ** 2)
    time domain: y(n) = x(n) + a * x(n - 1) + aa * x(n - 2)
    - b * y(n - 1) - bb * y(n - 2)
    c  x: input time series                                                 
    y: output time series                                                 
    n: length of x and y                                                  
    h: filter coefficients ; h(1)=a, h(2)=aa, h(3)=b, h(4)=bb             
    nml: >0; for normal direction filtering                              
         <0; for reverse direction filtering                             
    '''
    if (n <= 0):
        raise ValueError("Number of data points in raw is Wrong!")
    if (nml >= 0):
        j, jd = 0, 1  # index of data point and step
    else :
        j, jd = n - 1, -1

    a, aa, b, bb = h[0], h[1], h[2], h[3]
    u1, u2, v1, v2 = 0.0, 0.0, 0.0, 0.0
    for i in range(n):
        u3 = u2
        u2 = u1
        u1 = x[j]
        v3 = v2
        v2 = v1
        v1 = u1 + a * u2 + aa * u3 - b * v2 - bb * v3
        y[j] = v1
        j += jd
        # end for loop
    return x, y


def bandPass(x, dt, fl, fh, fs, ap, astop, nzero):
    '''
    band pass fileter, utilize scipy.signal
    x: raw data
    dt: time interval of raw data, sec
    fl: low frequency, Hz
    fh: high frequency, Hz
    fs: stop frequency, Hz
    ap: attenuation amplitude of pass frequency
    astop: attenuation amplitude of stop frequency
    nzero: number of zeros pad before raw data
    '''
    # pad zero before x
    y = np.concatenate((np.zeros((nzero, )), x))
    fn = int(1 / dt / 2)
    if (fl <= 0):
        sos = signal.iirdesign(wp=fh/fn, ws=fs/fn, gpass=ap, gstop=astop, \
                ftype='cheby2', output='sos')
    else :
        wpl = fl / fn; wph = fh / fn;
        wsl = fl * 0.9 / fn; wsh = fs / fn;
        sos = signal.iirdesign(wp=[wpl, wph], ws=[wsl, wsh], gpass=ap, gstop=astop, \
            ftype='cheby2', output='sos')
        # end if
    z = signal.sosfilt(sos, y)
    w = z[nzero:]
    return w



#------------------end filer section-----------------------------

def zoo2th(fzoo, fth):
    '''
    transfer file from ZOOFORMAT to timehistory
    '''
    with open(fzoo, 'r') as f:
        for i in range(4):
            f.readline()
        sampFreq = int(f.readline())
        numPoints = int(f.readline())
        f.seek(0)
        for i in range(10):
            f.readline()
        factor = float(f.readline())
        f.seek(0)
        for i in range(33):
            f.readline()
        tmp = f.read()
        tmp1 = tmp.split(maxsplit=numPoints)
        tmp2 = np.asfarray(tmp1)
        acc = tmp2 * factor
    dt = 1 / sampFreq
    t = np.linspace(dt, dt * numPoints, numPoints)
    with open(fth, 'w') as g:
        for i in range(numPoints):
            g.write(str(t[i]) + ' '  + str(acc[i]) + '\n')

