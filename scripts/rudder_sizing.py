import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
from tqdm import tqdm

from plain_flap_chord import Iterate
from avl_automation.tail import AutoTail

tail = AutoTail("../projects/tail_MDDP_v1.config")
tail.generate_planes()
tail.run()
_, curve_fit = tail.results(display=False)

# setup values
Lt, St_h, St_v = curve_fit.curve_fit_slice()

kf_data = read_csv("Kf_plot.csv")
angle = 37
chord_ratio_initial = 0.3
convergence = 0.005
span_ratio = 0.9
sweep = 25

AR_v = 5
thrust_moment_arm = 6.83  # m
thrust_per_side = 59400   # N
x_cg = 13.603
V1 = 37.597

# find chord ratios
iterate = Iterate()
chord_ratios = []
# for i,_ in tqdm(enumerate(St_v),total=len(St_v),desc="Calculating cf/c"):
for i, _ in enumerate(St_v):
    c_v = np.sqrt(St_v[i]/AR_v)
    b_v = np.sqrt(St_v[i]*AR_v)

    N = thrust_moment_arm*thrust_per_side

    x_tail = curve_fit.Lt_to_Xt(Lt[i])
    tail_moment_arm = x_tail-x_cg

    Lv = N/tail_moment_arm
    Clv = Lv/(0.5*1.225*St_v[i]*V1**2)

    try:
        r = iterate(chord_ratio_initial, convergence,
                    angle, kf_data, Clv, span_ratio, sweep)
        chord_ratios.append(r[-1])
    except ValueError:
        print(
            f"\u001b[33m[Warning]\u001b[0m Cl too high: {round(Clv,3)}. No possible elevator configuration.")
        chord_ratios.append(np.NaN)

# plot
fig, ax1, _ = curve_fit.plot_slice(Lt, St_h, St_v)

ax2 = ax1.twinx()
ax2.plot(Lt, chord_ratios, color='b', linestyle='--')

ax2.set_ylabel("cf/c")

plt.show()

"""
cf/c depends only on flap angle & tail volume coefficient
"""
