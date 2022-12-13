"""
Script to estimate mass of auto-tail configs. 
Empirical relations for mass estimation from 'Aircraft Design: A Conceptual Approach' by Raymer, D. used.

Run as python module.
"""

from tail_sizing import AutoTail
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

tail=AutoTail("projects/tail_MDDP_v0.config")
tail.generate_planes()
tail.run()
solutions_df=tail.results(display=False)

def calc_mass_v(HtHv,Wdg,Nz,Lt,Svt,Kz,delta_v,Av,tc):
    lb=0.0026*((1+HtHv)**0.225)*(Wdg**0.556)*(Nz**0.536)*(Lt**-0.5)*(Svt**0.5)*(Kz**0.875)*((np.cos(np.deg2rad(delta_v)))**-1)*(Av**0.35)*(tc**-0.5)
    kg=lb*2.205**-1

    return round(kg,2)

def calc_mass_h(Kuht,Fw,bh,Wdg,Nz,Sht,Lt,Ky,delta_h,Ah,Se):
    lb=0.0379*Kuht*((1+Fw/bh)**-0.25)*(Wdg**0.639)*(Nz**0.1)*(Sht**0.75)*(Lt**-1)*(Ky**0.704)*((np.cos(np.deg2rad(delta_h)))**-1)*(Ah**0.166)*(1+Se/Sht)**0.1
    kg=lb*2.205**-1

    return round(kg,2)

Av=5
delta_h=0       # deg
delta_v=0       # deg
Fw=5.5          # ft
HtHv=1
Kuht=1          
Nz=1.95         
tc=0.12     
Wdg=65506.31498 # lb
Xw=11.178       # m

mass_v=[]
mass_h=[]
for i,row in solutions_df.iterrows():
    Xt=row["Xt (Lunit)"]*3.281
    Lt=row["Lt (Lunit)"]*3.281
    Ky=0.3*Lt
    Sh=row["Sh (Lunit^2)"]*3.281**2
    Sv=row["Sv (Lunit^2)"]*3.281**2
    Se=0.3*Sh
    Kz=Lt
    Ah=row["ARh"]
    bh=np.sqrt(Sh*Ah)

    mass_h.append(calc_mass_h(Kuht,Fw,bh,Wdg,Nz,Sh,Lt,Ky,delta_h,Ah,Se))
    mass_v.append(calc_mass_v(HtHv,Wdg,Nz,Lt,Sv,Kz,delta_v,Av,tc))

solutions_df["mass_h (kg)"]=mass_h
solutions_df["mass_v (kg)"]=mass_v

print("")
print(solutions_df)

fig1,ax1=plt.subplots()
fig2,ax2=plt.subplots()

xs1=solutions_df["Sh (Lunit^2)"]
xs2=solutions_df["Sv (Lunit^2)"]
ys1=solutions_df["mass_h (kg)"]
ys2=solutions_df["mass_v (kg)"]

ax1.plot(xs1,ys1,color='k')
ax2.plot(xs2,ys2,color='k')

ax1.set_xlabel(r"$Sh$ ($m^2$)")
ax1.set_ylabel("Mass (Horizontal Stab.) (kg)")
ax2.set_xlabel(r"$Sv$ ($m^2$)")
ax2.set_ylabel("Mass (Vertical Stab.) (kg)")

plt.show()
