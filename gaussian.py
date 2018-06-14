from equadratures import *
import matplotlib.pyplot as plt
import numpy as np

def blackbox(x):
    return 3*x - 2.5

a = 3.0
b = 2.5
#------------------------ parameters for input distributio -------------#
sigma_quad = 3
std_dev = np.sqrt(sigma_quad)
mu = 2.0

# Using polynomials!

xo = Parameter(order=5, distribution='Gaussian', shape_parameter_A=mu, shape_parameter_B=sigma_quad)
myBasis = Basis('Tensor')
myPoly = Polyint([xo], myBasis)
myPoly.computeCoefficients(blackbox)
myStats = myPoly.getStatistics()
print 'Effective Quadratures'
print 'mean:',myStats.mean,'variance:', myStats.variance

# Using Monte Carlo!
#N = 900000
#xi = np.random.normal(mu, std_dev,  N)
#yi = evalfunction(np.reshape(xi, (N, 1) ), blackbox)
#print 'Monte Carlo'
#print np.mean(yi), np.var(yi)

x = np.linspace(-10,30.0,100)
sigma_quad = 3
std_dev = np.sqrt(sigma_quad)
mu = 2.0

def fx(x):
    primo = 1/(std_dev*np.sqrt(2*3.14))
    secondo = np.exp((-(x-mu)**2)/sigma_quad)
    #f = (1/(std_dev*np.sqrt(2*3.14)))*np.exp(-(x-mu)^2/(sigma_quad))
    return primo*secondo
    #return x

sigma_quad = sigma_quad*a**2
mu = mu*a +b

#print sigma_quad
#print mu
y = blackbox(x)
y_i = fx(x)

def fxo(x):
    primo = 1/(std_dev*np.sqrt(2*3.14))
    secondo = np.exp((-(x-mu)**2)/sigma_quad)
    return primo*secondo

#    return (1/(std_dev*np.sqrt(2*3.14)))*np.exp(-(x-mu)^2/(std_dev**2))
#
y_o = fxo(y)


plt.figure()
plt.grid()
#plt.plot(x,y_i,'r-', label='f_x (x)')
plt.plot(y,y_o,'b-',label='f_Y(y)')
plt.legend(loc='lower high')
#------ for output-----#
axes=plt.gca()
axes.set_xlim([-20,40])
plt.show()




# Plot the distribution!
#a, b = xo.getPDF(1000)
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#count, bins, ignored = plt.hist(xi, 30, normed=True)
#plt.plot(a, b)
#adjust_spines(ax, ['left', 'bottom'])
#plt.show()

