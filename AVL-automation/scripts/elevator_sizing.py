import numpy as np
from pandas import read_csv

from tail_sizing import AutoTail
from scripts.plain_flap_chord import Iterate

#### general parameters ####
Vstall=55.03
Cl_max=2.29
Cl0=0.273
Sw=95
Cd0=0.03    # zero lift drag

thrust=2*59400  # N
rolling_resistance=0.02

Cla_w=0.092
Cma_w=-0.022
Cla_h=0.0754
Cla_h_ground=1.1*Cla_h

i_h=-2  # incidence deg
i_w=0
mac_w=3.106 # m
AR_h=5

mass=30000  # kg
x_cg=11.86  # m
x_mg=13.39       # main gear m
z_mg=-2.761
z_D=0
#z_T=0.487   # m
z_T=0
z_cg=0.316

Ry=0.34
Lcb=455             # in
Lf=2.02*39.3701     # in
Lcp=157             # in
Lg=14               # in
Lep=56              # in
Lb=(Lcb+Lf+Lcp+Lg+Lep)*0.083    # ft
Lb_m=Lb*3.281**-1   # m

#### elevator parameters ####
ddtheta=10  # deg/s/s (rotation acceleration)
span_ratio_h=0.9
angle=30    # deg

#### initial calcs ####
Vr=1.2*Vstall

W=9.81*mass
Lto=0.5*1.225*Cl0*Sw*Vr**2  # lift at take-off N
D=0.5*1.225*Cd0*Sw*Vr**2    # drag at take-off N
friction=rolling_resistance*(Lto-W)
acceleration=(thrust-D-np.abs(friction))/mass

#### analysis stuff ####
tail=AutoTail("projects/tail_MDDP_v0.config")
tail.generate_planes()
tail.run()
_,curve_fit=tail.results(display=False)

Lt,St_h,St_v=curve_fit.curve_fit_slice()

plane=curve_fit.planes[0]
x_ac_w=plane.Xw_root+0.25*mac_w

kf_data=read_csv("scripts/Kf_plot.csv")
angle=30
chord_ratio_initial=0.3
convergence=0.005
span_ratio=0.8
sweep=0

iterate=Iterate()
chord_ratios=[]
for i,_ in enumerate(Lt):
    Xt=curve_fit.Lt_to_Xt(Lt[i])
    MAC_h=np.sqrt(St_h[i]/AR_h)

    L_emp_m=(Xt-Lb_m)+MAC_h # m
    L_fuselage=Lb_m+L_emp_m
    Iyy=((L_fuselage*3.281)**2*W*Ry**2)/(4*9.81)    # about cg?
    #print(L_fuselage,Iyy)
    Iyy=879000  # kg/m2

    x_ac_h=Xt+0.25*MAC_h

    M_ac_w=0.5*1.225*St_h[i]*Cma_w*Vr**2
    M_W=W*np.abs(x_mg-x_cg)
    M_D=D*np.abs(z_D-z_mg)
    M_T=thrust*np.abs(z_T-z_mg)
    M_Lw=Lto*np.abs(x_mg-x_ac_w)
    M_a=mass*acceleration*np.abs(z_cg-z_mg)

    #print(M_ac_w,-M_W,M_D,-M_T,M_Lw,M_a)
    #print(-M_W,M_D,-M_T,M_Lw,M_a)

    Lift_h=(-M_W+M_ac_w+M_a+M_Lw+M_D-M_T-Iyy*np.deg2rad(ddtheta))/(x_ac_h-x_mg)
    Clh=Lift_h/(0.5*1.225*St_h[i]*Vr**2)

    r=iterate(chord_ratio_initial,convergence,angle,kf_data,Clh,span_ratio,sweep)
    chord_ratios.append(round(r[-1],2))

print(chord_ratios)