###Project 3

#part 1:
import numpy as np
from scipy import integrate

integrate.solve_ivp()

mu = 2
def Ode(r,y,rhoval):
    rho, m = y
    if rho <= 0:
        return [0,0]
    
    x = rho**(1/3)
    gammaX = x**2/np.sqrt(3*(1+x**2))
    drho_dr = -m*rho / (gammaX * r**2)
    dm_dr = r**2 * rho
pc = np.linspace(0.1,2.5e6,10)