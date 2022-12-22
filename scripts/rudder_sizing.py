import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
from tqdm import tqdm

from scripts.plain_flap_chord import Iterate
from tail_sizing import AutoTail

tail = AutoTail("projects/tail_MDDP_v0.config")
tail.generate_planes()
tail.run()
_, curve_fit = tail.results(display=False)

# setup values
Lt, St_h, St_v = curve_fit.curve_fit_slice()

kf_data = read_csv("scripts/Kf_plot.csv")
angle = 35
chord_ratio_initial = 0.3
convergence = 0.005
span_ratio = 0.9
sweep = 0

AR_v = 5
thrust_moment_arm = 6.83  # m
thrust_per_side = 59400   # N
x_cg = 11.86
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

    r = iterate(chord_ratio_initial, convergence,
                angle, kf_data, Clv, span_ratio, sweep)
    chord_ratios.append(round(r[-1], 2))

# plot
fig, ax1, _ = curve_fit.plot_slice(Lt, St_h, St_v)

ax2 = ax1.twinx()
ax2.plot(Lt, chord_ratios, color='b', linestyle='--')

ax2.set_ylabel("cf/c")

plt.show()

"""
cf/c depends only on flap angle & tail volume coefficient
"""
