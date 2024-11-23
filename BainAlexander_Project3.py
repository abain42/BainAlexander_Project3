###Project 3

#part 1:
import numpy as np
from scipy import integrate

integrate.solve_ivp()

mu = 2
def Ode(r,y):
    rho, m = y
    if rho <= 0:
        return [0,0]
    
    x = rho**(1/3)
    gammaX = x**2/np.sqrt(3*(1+x**2))
    drho_dr = -m*rho / (gammaX * r**2)
    dm_dr = r**2 * rho
    return [drho_dr, dm_dr]
pc = np.linspace(0.1,2.5e6,10)
r_i = 0.00001
r_f = 10000000000
ans =[]

for rhoC in pc:
    sovler = integrate.solve_ivp(Ode,[r_i,r_f],[rhoC,0], method = "RK45",
                                 dense_output=True,events=lambda t, y:y[0])
    r_wd40 = sovler.t_events[0][0]
    m_wd40 = sovler.sol((r_wd40))[1]
    ans.append((rhoC,r_wd40,m_wd40))