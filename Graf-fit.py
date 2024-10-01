import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

data = pd.read_csv('A9-365.txt', sep='\s+', header=None)
data = pd.DataFrame(data)

x1 = data[1]
y1 = data[0]
#vzorec
xdata = []
ydata = []

for i in range(0, len(x1)-1):
    if x1[i] > 1.5 and x1[i] < 2.2:
        xdata.append(x1[i])
        ydata.append(y1[i])
    i = i + 1

#print(ydata)
#print(xdata)

fig = plt.figure()

ax1 = fig.add_subplot(111)

# cislo = stupen polynomu pri fitovani
z1, cov1 = np.polyfit(x1, y1, 5, cov=True)

#print(z1) # koeficienty prislusne danemu stupni

#print(np.sqrt(np.diag(cov1))) #vypisuje chybu fitu

p1 = np.poly1d(z1)


xp = np.linspace(0, 3, 100)

plt.xlabel(r'$ V \: (V) $ ')
plt.ylabel(r'$ I $ (nA) ')

def fun(x, a):
    y = a*(x+3)/(x+3)
    return y

def nul(x):
    y = 0*x
    return y

def fit(x, a, b):
    y = a*np.exp(b*x)
    return y

def slope(x, a, b):
    y = a*b*np.exp(b*x)
    return y



parameters, covariance = curve_fit(fit, xdata, ydata)
fit_a = parameters[0]
fit_b = parameters[1]


SE = np.sqrt(np.diag(covariance))
#print(np.diag(covariance))
SE_a = SE[0] #odchylka fitu
SE_b = SE[1] #odchylka fitu

xx = np.linspace(1.5, 2.2, 10000)

def klesajici(x):
    return -0.001*x

DataFit = fit(xx, fit_a, fit_b) - 0.01
index = []
for i in range(0, len(DataFit)):
    if DataFit[i] < 0:
        DataFit[i] = DataFit[i] - DataFit[i]
    i = i + 1

for i in xx:
    if ((fit(i, fit_a, fit_b)+fun(i, min(y1))) < 0.01 and ((fit(i, fit_a, fit_b)+fun(i, min(y1))) > -0.01)):
        bz = i
        break

x = bz
y = nul(x)

for i in xx:
    if ((fit(i, fit_a, fit_b)+fun(i, min(y1))) < 0.01 + min(y1)  and ((fit(i, fit_a, fit_b)+fun(i, min(y1))) > -0.01 + min(y1))):
        fn = i
        break

xbod = fn
ybod = fun(xbod, min(y1))
ybod2 = nul(xbod)

#plt.plot(xx, DataFit, 'y-', label='Fit')
#plt.plot(x1, y1, 'b.', label='Naměřené hodnoty (Hg Mon 365 nm)')
#plt.plot(xp, fun(xp, min(y1)), 'C1--', linewidth=1, label='Saturovaná hodnota proudu')
plt.plot(xx, nul(xx), 'k--', linewidth=1, label='Nulová hladina proudu')


plt.plot(xx, DataFit, 'g--', label='Katodový proud')
plt.plot(xx, fun(xx, min(y1)), 'C1--', label='Anodový proud')
plt.plot(xx, DataFit + fun(xx, min(y1)), 'b-', label='Celkový proud (naše měření)')

plt.plot(x, y, 'r.', markersize=10, label='Hodnota $V_0$')
plt.plot(xbod, ybod, 'C1.', markersize=10, label='')
plt.plot(xbod, ybod2, 'C1.', markersize=10, label='Extrapolovaná hodnota $V_0$')
#plt.text(x, y, round(x, 3), horizontalalignment='right',verticalalignment='top', color='r')

#plt.ylim(1.53,1.58)
#plt.xlim(1,2)

leg = ax1.legend()

plt.show()
