"""
Script to estimate mass of auto-tail configs. 
Empirical relations for mass estimation from 'Aircraft Design: A Conceptual Approach'
    by Raymer (2018) used (section on weight estimation).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from avlautomation.tail import AutoTail

def calc_mass_v(HtHv: float, Wdg: float, Nz: float, Lt: float, Svt: float, Kz: float, delta_v: float, Av: float, tc: float) -> float:
    """_summary_

    Args:
        HtHv (float): _description_
        Wdg (float): _description_
        Nz (float): _description_
        Lt (float): _description_
        Svt (float): _description_
        Kz (float): _description_
        delta_v (float): _description_
        Av (float): _description_
        tc (float): _description_

    Returns:
        float: _description_
    """
    lb = 0.0026*((1+HtHv)**0.225)*(Wdg**0.556)*(Nz**0.536)*(Lt**-0.5)*(Svt**0.5) * \
        (Kz**0.875)*((np.cos(np.deg2rad(delta_v)))**-1)*(Av**0.35)*(tc**-0.5)
    kg = lb*2.205**-1

    return round(kg, 2)


def calc_mass_h(Kuht, Fw, bh, Wdg, Nz, Sht, Lt, Ky, delta_h, Ah, Se):
    lb = 0.0379*Kuht*((1+Fw/bh)**-0.25)*(Wdg**0.639)*(Nz**0.1)*(Sht**0.75)*(Lt**-1)*(
        Ky**0.704)*((np.cos(np.deg2rad(delta_h)))**-1)*(Ah**0.166)*(1+Se/Sht)**0.1
    kg = lb*2.205**-1

    return round(kg, 2)


def calc_mass_fuselage(Kdoor, Klg, Wdg, Nz, L, Sf, Kws, D):
    lb = 0.3280*Kdoor*Klg*(Wdg*Nz)**0.5*L**0.25 * \
        Sf**0.302*(1+Kws)**0.04*(L/D)**0.1
    kg = lb*2.208**-1

    return round(kg, 2)


Av = 1.8
D = 136           # ft
delta_w = 0       # deg   sweep
delta_h = 0       # deg
delta_v = 33       # deg
Fw = 5.67          # ft
HtHv = 1
Kdoor = 1
Klg = 1.12
Kuht = 1
lambda_w = 0.4    # taper ratio
Nz = 3
Sf = 3.51E+03 # ft^2
tc = 0.12
Wdg = 66138.6  # lb

Lb = 19.354*0.083    # ft

tail = AutoTail("../projects/tail_MDDP_v1.config")
tail.generate_planes()
tail.run()
_, curve_fit = tail.results(display=False)
Lts, St_hs, St_vs = curve_fit.curve_fit_slice()
Xts=curve_fit.Lt_to_Xt(Lts) # m

mass_v = []
mass_h = []
mass_f = []
mass_total = []
for i, _ in enumerate(Lts):
    Xt = curve_fit.Lt_to_Xt(Lts[i])*3.281  # ft

    Lt = Lts[i]*3.281  # ft
    Ky = 0.3*Lt
    Sh = St_hs[i]*3.281**2    # ft^2
    Sv = St_vs[i]*3.281**2
    Se = 0.3*Sh
    Kz = Lt

    Ah = 5
    bh = np.sqrt(Sh*Ah)
    Bw = curve_fit.planes[i].b_w

    mass_h.append(calc_mass_h(Kuht, Fw, bh, Wdg,
                  Nz, Sh, Lt, Ky, delta_h, Ah, Se))
    mass_v.append(calc_mass_v(HtHv, Wdg, Nz, Lt, Sv, Kz, delta_v, Av, tc))

    MAC_h = np.sqrt(Sh/Ah)
    L_emp_m = (Xt*3.281**-1-Lb*3.281**-1)+MAC_h  # m
    L_emp_ft = L_emp_m*3.281  # ft
    L = Lb+L_emp_ft   # ft

    Kws = 0.75*((1+2*lambda_w)/(1+lambda_w))*(Bw*np.tan(np.deg2rad(delta_w))/L)

    mass_f.append(calc_mass_fuselage(Kdoor, Klg, Wdg, Nz, L, Sf, Kws, D))

    mass_total.append(mass_h[i]+mass_v[i]+mass_f[i])

mass_total=np.array(mass_total)
mass_f=np.array(mass_f)
mass_h=np.array(mass_h)
mass_v=np.array(mass_v)

# print results
res=pd.DataFrame(columns=["Lt","Sh","Sv","mass_h","mass_v","mass_f","mass_total"])
res["Lt"]=Lts
res["Xt"]=Xts
res["Sh"]=St_hs
res["Sv"]=St_vs
res["mass_h"]=mass_h
res["mass_v"]=mass_v
res["mass_f"]=mass_f
res["mass_total"]=mass_total

print(res.iloc[res["mass_total"].idxmin()])

# plotting
fig, ax1, _ = curve_fit.plot_slice(Lts, St_hs, St_vs)
ax2 = ax1.twinx()

ax2.plot(Lts, mass_total, linestyle='--', color='k')
# ax2.plot(Lts,mass_h+mass_v,linestyle=':',color='k')
# ax2.plot(Lts,mass_f,linestyle='-.',color='k')
ax2.set_ylabel("Tail + Fuselage Mass (kg)")

plt.show()
