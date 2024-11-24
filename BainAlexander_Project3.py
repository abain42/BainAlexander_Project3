###Project 3

#part 1:
#import libraries
import numpy as np
from scipy import integrate


#this is the function that we created to solve our ODE's
def Ode(r,y):
    rho, m = y
    if rho <= 0:
        #stops integration if density goes negative 
        return [0,0]
    #this is a dimensionless variable given to us 
    x = rho**(1/3)
    #EQ given to us
    gammaX = x**2/(3*np.sqrt((1+x**2)))
    drho_dr = -m*rho / (gammaX * r**2)
    dm_dr = r**2 * rho
    return [drho_dr, dm_dr]
#this sets our variables for the question, logspace was used instead of linspace because 
#it will space the data better when on the graph
pc = np.logspace(-1,6,10)
#arbitraty start and stop points, start point isnt zero do avoid division by 0
#0.00001 was chosed because it felt small enough that it would 
r_i = 0.00001
r_f = 10000000000
ans =[]

#this creates a for loop with length, the amount of rhoc's we have, which is our initial conditions

for rhoC in pc:
    #this is the solver_ivp function that we use to solve the ODE's
    sovler = integrate.solve_ivp(Ode,[r_i,r_f],[rhoC,0], method = "RK45",
                                 dense_output=True,events=lambda t, y:y[0])
    #gives us radius and mass when density = 0 and appends to ans 
    r_wd40 = sovler.t_events[0][0]
    m_wd40 = sovler.sol((r_wd40))[1]
    ans.append((rhoC,r_wd40,m_wd40))

### part 2
#imported constants and plot
from astropy.constants import G, M_sun, R_sun, m_e,m_p,h
import matplotlib.pyplot as plt
#here are constants given to us 
mu_e = 2
R0 = 7.72e8 / mu_e 
M0 = 5.67e33 / (mu_e**2)
#turning values into appropriate units, kg to g, m to cm
msungrams = M_sun.value*1000 
rsuncm =R_sun.value*100
#mass in solar masses and radius in solar radii
masses = [m * M0 / msungrams for _, _, m in ans]  
radii = [r * R0 / rsuncm for _, r, _ in ans]  

#this pplots the relationship between the two
plt.figure(figsize=(8, 6))
plt.plot(masses, radii, marker="o", label="Numerical Results")
plt.xlabel("Mass ($M_\odot$)")
plt.ylabel("Radius ($R_\odot$)")
plt.title("Mass-Radius Relation for White Dwarfs")
plt.legend()
plt.grid()
plt.show()



#Part 3: 
#this is the same as part 1, however we have only 3 rho_c and we are using different
#solver_ivp methods of integration
#this is not the most efficient way to do this, if we wanted to be efficient 
#we would combine this with part 1, (same can be said for part 4 with part 2)
#but we have purposely seperated these so it is more clear what work we did for each part
rho_c_valuespt2 = [1e2, 1e4, 1e6]
methods = ['DOP853', 'RK23', 'BDF'] 
results = {}
#this creates a for loop of length, the amount of integration methods we are trying
for i in methods:
    results[i] = []
    #this is the exact same method as part 1, we are just doing it multiple times, for each 
    #method of integration as mention before with the for loop 
    for rhoC in rho_c_valuespt2:
        sovler = integrate.solve_ivp(Ode,[r_i,r_f],[rhoC,0], method = i,
                                 dense_output=True,events=lambda t, y:y[0])
    
        r_wd40 = sovler.t_events[0][0]
        m_wd40 = sovler.sol((r_wd40))[1]
        results[i].append((rhoC,r_wd40,m_wd40))


plt.figure(figsize=(12, 8))
#the plot is essentially the same as before however we need a for loop to plot for each method
for i in methods:
    masses_new = [m * M0 / msungrams for _, _, m in results[i]]
    radii_new = [r * R0 / rsuncm for _, r, _ in results[i]]
    plt.plot(masses_new, radii_new, marker="o", label=f"Method: {i}")

plt.xlabel("Mass ($M_\odot$)")
plt.ylabel("Radius ($R_\odot$)")
plt.title("Answers for new SOLVE_IVP methods")
plt.legend()
plt.grid()
plt.show()



#part 4
#part 4 we import the file given to us from wherever it is locally on my computer
#note: this will be unique to each computer
import pandas as pd
file_location = r"C:\Users\Alex\Downloads\wd_mass_radius.csv"
#read the data and added it to a dataset
data = pd.read_csv(file_location)
#this makes it start on row 2
datar2 = data.iloc[1:]
#each of these collect each collumn startin on row 2
datamass = datar2.iloc[:,0]
datarad = datar2.iloc[:, 2]
mass_err = datar2.iloc[:, 1]
radii_err = datar2.iloc[:, 3]
#we then plot these, against our original plot in part 1
# you will notice that the error bars are quite large ,
# it says that the error is in units of solar radii, but I am wondering if it is meant to be 
# a percent error instead
plt.figure(figsize=(8, 6))
plt.plot(masses, radii, marker="o", label="Numerical Results")
plt.errorbar(datamass,datarad,xerr =mass_err,yerr = radii_err, fmt ='o',label='Data from Tremblay et al',capsize=5)
plt.xlabel("Mass ($M_\odot$)")
plt.ylabel("Radius ($R_\odot$)")
plt.title("Mass-Radius Relation for White Dwarfs")
plt.legend()
plt.grid()
plt.show()
