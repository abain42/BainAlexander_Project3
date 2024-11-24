###Project 3

#part 1:
import numpy as np
from scipy import integrate


mu = 2
def Ode(r,y):
    rho, m = y
    if rho <= 0:
        return [0,0]
    
    x = rho**(1/3)
    gammaX = x**2/(3*np.sqrt((1+x**2)))
    drho_dr = -m*rho / (gammaX * r**2)
    dm_dr = r**2 * rho
    return [drho_dr, dm_dr]
pc = np.logspace(-1,6,10)
r_i = 0.00001
r_f = 10000000000
ans =[]

for rhoC in pc:
    sovler = integrate.solve_ivp(Ode,[r_i,r_f],[rhoC,0], method = "RK45",
                                 dense_output=True,events=lambda t, y:y[0])
    
    r_wd40 = sovler.t_events[0][0]
    m_wd40 = sovler.sol((r_wd40))[1]
    ans.append((rhoC,r_wd40,m_wd40))

### part 2
from astropy.constants import G, M_sun, R_sun, m_e,m_p,h
import matplotlib.pyplot as plt
print(M_sun)
mu_e = 2
R0 = 7.72e5 / mu_e 
M0 = 5.67e33 / (mu_e**2)
msungrams = M_sun.value*1000 
rsuncm =R_sun.value
masses = [m * M0 / msungrams for _, _, m in ans]  
radii = [r * R0 / rsuncm for _, r, _ in ans]  
print(R_sun)

plt.figure(figsize=(8, 6))
plt.plot(masses, radii, marker="o", label="Numerical Results")
plt.xlabel("Mass ($M_\odot$)")
plt.ylabel("Radius ($R_\odot$)")
plt.title("Mass-Radius Relation for White Dwarfs")
plt.legend()
plt.grid()
plt.show()
