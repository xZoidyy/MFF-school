import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.optimize import curve_fit
import scipy.integrate as integrate

data = pd.read_csv('deuteron.txt', sep = '\s+', header=None)
data = pd.DataFrame(data)

Hist = data[1]

velikost = len(Hist)

#########################
#####    Binning  ########
#########################
bins = np.arange(0, (velikost+1), 1)
binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)]) ## prostredek binu pro fit

#plt.hist(bins[:-1], bins=bins, weights=Hist, label='hist') alternativa plotovani


fig = plt.figure()
ax1 = fig.add_subplot(111)
#plt.xlabel('E $ [keV] $')
plt.xlabel('Channel')
plt.ylabel("Counts")

##### DEF FUNKCI ########

def Gaus(x, m, s):
  N = s * np.sqrt(np.pi/2) * (erf((m-r1)/np.sqrt(2)/s) - erf((m-r2)/np.sqrt(2)/s))
  return 1./N * np.exp(-(x-m)**2/(2*s**2))

def Gauss(x, m, s, N):
  return N * np.exp(-1.0 * (x - m)**2 / (2 * s**2))

def Lin(x, K, B):
  return K*x + B  #(1 + K*(x-(r1+r2)/2) ) / (r2 - r1)

def Bkg(x, m, s, K, A, B):
  return A * np.exp(-1.0 * (x - m)**2 / (2 * s**2)) + K*x + B #B * ((1 + K*(x-(r1+r2)/2) ) / (r2 - r1))
##############################################

mu = []
sig = []

#### ROI + FITTING ####
rstart = [955, 1770, 2796, 3500] #3529
rstop = [974, 1795, 2830, 3530] #3569

for j in range(0, len(rstart)): #len(rstart)
  r1, r2 = 0, 0
  r1 = rstart[j]
  r2 = rstop[j]
  val, cov = 0, 0

  ROIbins = []
  ROIhist = []

  for i in range(len(binscenters)-1):
    if binscenters[i] >= r1 and binscenters[i] <= r2:
      ROIbins.append(binscenters[i])
      ROIhist.append(Hist[i])

  val, cov = curve_fit(Bkg, xdata=ROIbins, ydata=ROIhist, p0=[(r1+r2)/2, 1, -0.1, 10000, 700])
  SE = np.sqrt(np.diag(cov))
  #print(val)
  mu.append(val[0])
  sig.append(SE[0])

  xp = np.linspace(r1, r2, 10000)

  ##### pro plot jen jednoho peaku se spektrem #######
  #plt.bar(ROIbins, ROIhist, width=bins[1] - bins[0], color='navy', label=r'Naměřená data')

  #plt.plot(xp, Bkg(xp, *val), color='darkorange', linewidth=1, label=r'Fitovaná závislost')

  ###### text s mean value Gaussiánu u peaků ######
  #plt.text(val[0]-6, val[3], '$\mu$ = ' + str(round(val[0], 2)), horizontalalignment='center', verticalalignment='bottom', color='darkorange', label='Mean of Gaussian')
  #plt.text(val[0]-6, 1860, '$\sigma_{\mu}$ = ' + str(round(SE[0], 2)), horizontalalignment='center', verticalalignment='bottom', color='darkorange', label='Mean of Gaussian')

###################################################

print(F'The value of Mean of Gaussian is {mu}.')
print(F'Standard error of the Mean value of Gaussian is {sig}.')

#plt.bar(binscenters, Hist, width=bins[1] - bins[0], color='navy', label=r'Naměřená data')

plt.yscale('log')
leg = ax1.legend()
plt.show()

############################################
###     KALIBRACE     ####
###########################################

fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.ylabel('E $ [keV] $')
plt.xlabel('Střední hodnota gaussiánu ($\mu$)')

E = [609.34, 1120.29, 1764.50, 2204.10]

def fun(x, a, b):
  return a*x + b

def fun2(E, a, b): # protože curve_fit umí jen chyby u y-osy tak to invertuji
  return (E-b)/a

parameters, covariance = curve_fit(fun2, E, mu, sigma=sig)
#parameters, covariance = curve_fit(fun, mu, E)

fit_a = parameters[0]
fit_b = parameters[1]

#print(fit_c) # hodnota fitu

SE = np.sqrt(np.diag(covariance))
SE_a = SE[0] #odchylka fitu
SE_b = SE[1] #odchylka fitu

print(F'The value of a is {fit_a} with standard error of {SE_a:}.')
print(F'The value of b is {fit_b} with standard error of {SE_b:}.')

xx = np.linspace(950, 3536, 10000)
xd = np.linspace(600, 2210, 10000)

#####plt.plot(mu, E, 'g.', linewidth=1, label=r'')
#plt.plot(xx, fun(xx, fit_a, fit_b), 'g--', linewidth=1, label=r'Fitovaná závislost: $E_{\gamma} = a \cdot Channel + b$')

####plt.plot(xd, fun2(xd, fit_a, fit_b), 'g--', linewidth=1, label=r'')
#plt.errorbar(mu, E, yerr=sig, fmt='o', markersize=4, capsize=4, color='black', label='Střední hodnoty gaussiánů ($\mu$) přiřazené \ntabelovaným hodnotám charakteristického \nzáření $^{226}$Ra')

leg = ax1.legend()
plt.show()


################################################
##### FINAL HISTOGRAM (CALIBRATED) ############
###############################################

a = 0.624
b = 7.673

fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.xlabel('E $ [keV] $')
plt.ylabel('Četnost')

Ebins = []

for i in range(0, len(bins)):
  Ebins.append(fit_a*bins[i]+fit_b)

Ebinscenters = np.array([0.5 * (Ebins[i] + Ebins[i+1]) for i in range(len(Ebins)-1)])

#vyska1 = np.linspace(0, 46218, 100000)
#point1 = np.full(100000, 468.82)
#vyska2 = np.linspace(0, 46218, 100000)
#point2 = np.full(100000, 485.85)

#vyska3 = np.linspace(0, 730, 100000)
#point3 = np.full(100000, 4323.3)
#vyska4 = np.linspace(0, 730, 100000)
#point4 = np.full(100000, 4539.4)

#plt.bar(Ebins, Hist, width=Ebins[1] - Ebins[0], color='navy', label=r'Naměřená data')

##plt.plot(point1, vyska1, 'r-')
##plt.text(455.70, 29500, '$E_{min}^{Li}$ = 468,82 keV', horizontalalignment='center', verticalalignment='bottom', color='red', label='Mean of Gaussian')
##plt.plot(point2, vyska2, 'r-')
#plt.text(499.05, 29500, '$E_{max}^{Li}$ = 485,85 keV', horizontalalignment='center', verticalalignment='bottom', color='red', label='Mean of Gaussian')

#plt.plot(point3, vyska3, 'r-')
#plt.text(4323.3, 900, '$E_{min}^{C}$ = 4323,3 keV', horizontalalignment='center', verticalalignment='bottom', color='red', label='Mean of Gaussian')
#plt.plot(point4, vyska4, 'r-')
#plt.text(4539.4, 900, '$E_{max}^{C}$ = 4539,4 keV', horizontalalignment='center', verticalalignment='bottom', color='red', label='Mean of Gaussian')


mu = []
sig = []
K = []
B = []
N = []

#### ROI + FITTING ####
rstart = [2208] #3529
rstop = [2236] #3569

for j in range(0, len(rstart)): #len(rstart)
  r1, r2 = 0, 0
  r1 = rstart[j]
  r2 = rstop[j]
  val, cov = 0, 0

  ROIbins = []
  ROIhist = []

  for i in range(len(binscenters)-1):
    if Ebins[i] >= r1 and Ebins[i] <= r2:
      ROIbins.append(Ebins[i])
      ROIhist.append(Hist[i])

  val, cov = curve_fit(Bkg, xdata=ROIbins, ydata=ROIhist, p0=[(r1+r2)/2, 1, -0.1, 10000, 700])
  SE = np.sqrt(np.diag(cov))
  #print(val)
  mu.append(val[0])
  sig.append(SE[0])
  K.append(val[2])
  B.append(val[4])

  xp = np.linspace(r1, r2, 10000)

  ##### pro plot jen jednoho peaku se spektrem #######
  plt.bar(ROIbins, ROIhist, width=bins[1] - bins[0], color='navy', label=r'Naměřená data')

  plt.plot(xp, Bkg(xp, *val), color='darkorange', linewidth=2, label=r'Fitovaná závislost')
  #plt.plot(xp, Lin(xp, K[0], B[0]), 'g-', linewidth=2, label=r'Fitovaná závislost')
  #plt.plot(xp, Gauss(xp, K[0], B[0]), 'g-', linewidth=2, label=r'Fitovaná závislost')


  ###### text s mean value Gaussiánu u peaků ######
  plt.text(val[0]-6, 4200, '$\mu$ = ' + str(round(val[0], 2)) + ' $\pm$ ' + str(round(SE[0], 2)), horizontalalignment='center', verticalalignment='bottom', color='darkorange', label='Mean of Gaussian')
  plt.text(val[0]-6, 3900, '$\sigma$ = ' + str(round(val[1], 2)) + ' $\pm$ ' + str(round(SE[1], 2)), horizontalalignment='center', verticalalignment='bottom', color='darkorange', label='Mean of Gaussian')

###################################################

print(F'The value of Mean of Gaussian is {mu}.')
print(F'Standard error of the Mean value of Gaussian is {sig}.')

print(F'The value of sigma of Gaussian is {val[1]}.')
print(F'Standard error of the sigma of Gaussian is {SE[1]}.')

AreaBkg = integrate.quad(lambda x: Bkg(x, *val), r1, r2)
AreaLin = integrate.quad(lambda x: Lin(x, K[0], B[0]), r1, r2)

Area = AreaBkg[0] - AreaLin[0]
print(AreaBkg)
print(AreaLin)
erorArea = np.sqrt(AreaBkg[1]**2 + AreaLin[1]**2)
print(F'Area of peak is {Area} with error {erorArea}')

plt.text(val[0] - 6, 3600, 'Area = ' + str(round(Area, 3)), horizontalalignment='center', verticalalignment='bottom', color='darkorange', label='Mean of Gaussian')

#plt.yscale('log')
leg = ax1.legend()
plt.show()

