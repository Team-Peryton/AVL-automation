import numpy as np
from pandas import read_csv
import matplotlib.pyplot as plt

from avlautomation.tail import AutoTail
from plain_flap_chord import Iterate

#### general parameters ####
Vstall = 55.03
Cl_max = 2.29
Cl0 = 0.273
Sw = 95
Cd0 = 0.03    # zero lift drag

thrust = 2*59400  # N
rolling_resistance = 0.02

Cla_w = 0.092
Cma_w = -0.022
Cla_h = 0.0754
Cla_h_ground = 1.1*Cla_h

i_h = -2  # incidence deg
i_w = 0
mac_w = 3.106  # m
AR_h = 5

mass = 30000  # kg
x_cg = 13.603  # m
x_mg = 15.03       # main gear m
z_mg = -2.81
z_D = 0
z_T = 1.487   # m
# z_T=0
z_cg = 0.6

Iyy = 1359617.415  # kg/m2

#### elevator parameters ####
ddtheta = 10  # deg/s/s (rotation acceleration)
span_ratio = 0.9
angle = 35    # deg
sweep = 0

#### initial calcs ####
Vr = 1.2*Vstall

W = 9.81*mass
Lto = 0.5*1.225*Cl0*Sw*Vr**2  # lift at take-off N
D = 0.5*1.225*Cd0*Sw*Vr**2    # drag at take-off N
friction = rolling_resistance*(Lto-W)
acceleration = (thrust-D-np.abs(friction))/mass

#### analysis stuff ####
tail = AutoTail("../projects/tail_MDDP_v1.config")
tail.generate_planes()
tail.run()
_, curve_fit = tail.results(display=False)

Lt, St_h, St_v = curve_fit.curve_fit_slice()

plane = curve_fit.planes[0]
x_ac_w = plane.Xw_root+0.25*mac_w

kf_data = read_csv("Kf_plot.csv")
chord_ratio_initial = 0.3
convergence = 0.005

iterate = Iterate()
chord_ratios = []
for i, _ in enumerate(Lt):
    Xt = curve_fit.Lt_to_Xt(Lt[i])
    MAC_h = np.sqrt(St_h[i]/AR_h)

    x_ac_h = Xt+0.25*MAC_h

    M_ac_w = 0.5*1.225*St_h[i]*Cma_w*Vr**2
    M_W = W*np.abs(x_mg-x_cg)
    M_D = D*np.abs(z_D-z_mg)
    M_T = thrust*np.abs(z_T-z_mg)
    M_Lw = Lto*np.abs(x_mg-x_ac_w)
    M_a = mass*acceleration*np.abs(z_cg-z_mg)

    # print(M_ac_w,-M_W,M_D,-M_T,M_Lw,M_a)
    # print(-M_W,M_D,-M_T,M_Lw,M_a)

    Lift_h = np.abs((-M_W+M_ac_w+M_a+M_Lw+M_D-M_T-Iyy *
                    np.deg2rad(ddtheta))/(x_ac_h-x_mg))
    Clh = Lift_h/(0.5*1.225*St_h[i]*Vr**2)

    try:
        r = iterate(chord_ratio_initial, convergence,
                    angle, kf_data, Clh, span_ratio, sweep)
        chord_ratios.append(r[-1])
    except ValueError:
        print(
            f"\u001b[33m[Warning]\u001b[0m Cl too high: {round(Clh,3)}. No possible elevator configuration.")

        chord_ratios.append(np.NaN)
        #print(f"L={round(Lift_h,2)}")

fig, ax1, _ = curve_fit.plot_slice(Lt, St_h, St_v)
ax2 = ax1.twinx()

ax2.plot(Lt, chord_ratios, linestyle='--', color='r')
ax2.set_ylabel("cf/c")

plt.show()
