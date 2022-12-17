"""
Script to estimate mass of auto-tail configs. 
Empirical relations for mass estimation from 'Aircraft Design: A Conceptual Approach' by Raymer, D. used.

Run as python module.
"""

from tail_sizing import AutoTail
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def calc_mass_v(HtHv,Wdg,Nz,Lt,Svt,Kz,delta_v,Av,tc):
    lb=0.0026*((1+HtHv)**0.225)*(Wdg**0.556)*(Nz**0.536)*(Lt**-0.5)*(Svt**0.5)*(Kz**0.875)*((np.cos(np.deg2rad(delta_v)))**-1)*(Av**0.35)*(tc**-0.5)
    kg=lb*2.205**-1

    return round(kg,2)

def calc_mass_h(Kuht,Fw,bh,Wdg,Nz,Sht,Lt,Ky,delta_h,Ah,Se):
    lb=0.0379*Kuht*((1+Fw/bh)**-0.25)*(Wdg**0.639)*(Nz**0.1)*(Sht**0.75)*(Lt**-1)*(Ky**0.704)*((np.cos(np.deg2rad(delta_h)))**-1)*(Ah**0.166)*(1+Se/Sht)**0.1
    kg=lb*2.205**-1

    return round(kg,2)

def calc_mass_fuselage(Kdoor,Klg,Wdg,Nz,L,Sf,Kws,D):
    lb=0.3280*Kdoor*Klg*(Wdg*Nz)**0.5*L**0.25*Sf**0.302*(1+Kws)**0.04*(L/D)**0.1
    kg=lb*2.208**-1

    return round(kg,2)

Av=5
D=136           # ft
delta_w=0       # deg   sweep
delta_h=0       # deg
delta_v=0       # deg
Fw=5.5          # ft
HtHv=1
Kdoor=1
Klg=1.12
Kuht=1  
lambda_w=0.4    # taper ratio
Nz=4.5
Sf=2294.3       # ft^2         
tc=0.12     
Wdg=65506.31498 # lb
Xw=11.178       # m

Lcb=455             # in
Lf=2.02*39.3701     # in
Lcp=157             # in
Lg=14               # in
Lep=56              # in
Lb=(Lcb+Lf+Lcp+Lg+Lep)*0.083    # ft

tail=AutoTail("projects/tail_MDDP_v0.config")
tail.generate_planes()
tail.run()
solutions_df=tail.results(display=False)
planes=tail.get_planes()

mass_v=[]
mass_h=[]
mass_f=[]
mass_total=[]
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
    Bw=planes[i].b_w

    mass_h.append(calc_mass_h(Kuht,Fw,bh,Wdg,Nz,Sh,Lt,Ky,delta_h,Ah,Se))
    mass_v.append(calc_mass_v(HtHv,Wdg,Nz,Lt,Sv,Kz,delta_v,Av,tc))

    MAC_h=np.sqrt(row["Sh (Lunit^2)"]/Ah)
    L_emp=Lt-(Lb-planes[i].Xw_root-planes[i].mac/4)+MAC_h*(3/4)
    L=Lb+L_emp
    Kws=0.75*((1+2*lambda_w)/(1+lambda_w))*(Bw*np.tan(np.deg2rad(delta_w))/L)

    mass_f.append(calc_mass_fuselage(Kdoor,Klg,Wdg,Nz,L,Sf,Kws,D))

    mass_total.append(mass_h[i]+mass_v[i]+mass_f[i])

solutions_df["mass_h (kg)"]=mass_h
solutions_df["mass_v (kg)"]=mass_v
solutions_df["mass_f (kg)"]=mass_f
solutions_df["mass_total (kg)"]=mass_total
solutions_df["S_tail_sum (Lunit^2)"]=solutions_df["Sh (Lunit^2)"]+solutions_df["Sv (Lunit^2)"]

solutions_df.sort_values(by=["mass_total (kg)"],ascending=True,inplace=True)

print("")
print(solutions_df)

##### 3D PLOTS #####
"""
fig=plt.figure(figsize=(6,6))
ax=plt.axes(projection='3d')

xs=solutions_df["S_tail_sum (Lunit^2)"]
ys=solutions_df["Lt (Lunit)"]

zs1=solutions_df["mass_h (kg)"]+solutions_df["mass_v (kg)"]
zs2=solutions_df["mass_f (kg)"]
zs3=solutions_df["mass_total (kg)"]

ax.scatter(xs,ys,zs1,edgecolors='r',facecolors='none',depthshade=False,label="Tail")
ax.scatter(xs,ys,zs2,edgecolors='b',facecolors='none',depthshade=False,label="Fuselage")
ax.scatter(xs,ys,zs3,edgecolors='k',facecolors='none',depthshade=False,label="Combined")

ax.set_xlabel(r"St ($Lunit^2$)")
ax.set_ylabel(r"$Lt$ ($Lunit$)")
ax.set_zlabel("Mass (kg)")

plt.show()
"""
##### 2D PLOTS #####

fig1,ax1=plt.subplots()
fig2,ax2=plt.subplots()

xs1=solutions_df["S_tail_sum (Lunit^2)"]
ys1=solutions_df["mass_h (kg)"]+solutions_df["mass_v (kg)"]
ax1.scatter(xs1,ys1,edgecolors='k',facecolors='none')

ax1.set_xlabel(r"St ($Lunit^2$)")
ax1.set_ylabel("Tail Mass (kg)")


xs2=solutions_df["Lt (Lunit)"]
ys2=solutions_df["mass_f (kg)"]
ax2.scatter(xs2,ys2,edgecolors='r',facecolors='none',label="Fuselage")
ax2.scatter(xs2,ys1,edgecolors='b',facecolors='none',label="Tail")

ax2.set_xlabel(r"$Lt$ ($Lunit$)")
ax2.set_ylabel("Mass (kg)")
ax2.legend()

# plt.show()
