import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import optimize
from warnings import filterwarnings


class Iterate():
    filterwarnings(
        "ignore", message="Covariance of the parameters could not be estimated")

    def __call__(self, chord_ratio_initial: float, convergence_criteria: float,
                 angle: float, kf_data: pd.DataFrame, dCl: float, wet_area_ratio: float,
                 sweep: float):
        """
        Calculates chord ratio for a plain flap.

        Arguments:
            chord_ratio_initial: {float} -- Initial guess for cf/c
            convergence_criteria: {float} -- 
            angle: {float} -- Flap angle (degrees)
            kf_data: {pd.DataFrame} -- Angle, Kf, cf/c plot data (Raymer, 2018)
            dCl: {float} -- Required delta flap Cl (3D)
            wet_area_ratio: {float} -- Ratio of wetted area of flap to wing (Raymer, 2018)
            sweep: {float} -- Surface sweep angle (degrees)

        Returns:
            chord_ratios: {list[float]} -- List of chord ratios across iterations. [-1] is converged value
        """

        return self.iterate(chord_ratio_initial, convergence_criteria,
                            angle, kf_data, dCl, wet_area_ratio,
                            sweep)

    def func_alpha_Kf(self, data, a, b, c, d, e, f, g, h, i):
        x = data
        b = 1
        return a*(g*np.arctan((c*x)**d)*(e*x)**f)+(h*x)**i+b

    def curve_fit(self, data, x):
        np.seterr(all="ignore")
        parameters, _ = optimize.curve_fit(
            self.func_alpha_Kf, data[:, 0], data[:, 1],
            (-0.2, 1, 0.1, 1.5, 0, 1, 1, 0.001, 2)
        )

        y = self.func_alpha_Kf(x, *parameters)

        np.seterr(all="warn")
        return y, parameters

    def thickness_ratio_12(self, cl_increment):
        a = -13.738
        b = 19.397
        c = 0.5337-cl_increment

        x = (-b+np.sqrt(b**2-4*a*c))/(2*a)   # quadratic equation
        return x

    def iterate(self, chord_ratio_initial: float, convergence_criteria: float,
                angle: float, kf_data: pd.DataFrame, dCl: float, wet_area_ratio: float,
                sweep: float):

        ratio_unique = kf_data["cf/c"].unique()
        ratio_unique.sort()
        chord_ratio = [chord_ratio_initial]
        Kf = []
        cl_increment = []
        n = 0
        while True:

            i = np.searchsorted(ratio_unique, chord_ratio[n], side="left")
            if i == len(ratio_unique):
                i -= 1
            ratio_left = ratio_unique[i-1]
            ratio_right = ratio_unique[i]

            data_left = kf_data[kf_data["cf/c"] ==
                                ratio_left][["Angle", "Kf"]].to_numpy()
            data_right = kf_data[kf_data["cf/c"] ==
                                 ratio_right][["Angle", "Kf"]].to_numpy()

            Kf_left, parameters_left = self.curve_fit(data_left, angle)
            Kf_right, parameters_right = self.curve_fit(data_right, angle)

            Kf.append(
                np.interp(chord_ratio[n], (ratio_left, ratio_right), (Kf_left, Kf_right)))

            cl_increment_ = dCl / \
                (0.9*Kf[n]*wet_area_ratio*np.deg2rad(angle)
                 * np.cos(np.deg2rad(sweep)))
            if cl_increment_ > 6:
                raise ValueError("Cl out of range.")

            cl_increment.append(cl_increment_)
            cr_ = self.thickness_ratio_12(cl_increment[n])
            if cr_ != 0:
                chord_ratio.append(cr_)
            else:
                chord_ratio.append(0)

            """ plot for debugging
            angle_range=np.linspace(10,60,100)
            kf_left=func_alpha_Kf(angle_range,*parameters_left)
            kf_right=func_alpha_Kf(angle_range,*parameters_right)
            
            fig,ax=plt.subplots()
            ax.plot(angle_range,kf_left)
            ax.plot(angle_range,kf_right)
            ax.scatter(data_left[:,0],data_left[:,1])
            ax.scatter(data_right[:,0],data_right[:,1])

            ax.scatter((angle,angle,angle),(Kf_left,Kf[n],Kf_right))

            plt.show()
            """

            if np.isclose(chord_ratio[n+1], chord_ratio[n], atol=convergence_criteria):
                break
            n += 1

        return chord_ratio


if __name__ == "__main__":
    kf_data = pd.read_csv("scripts/Kf_plot.csv")
    angle = 30
    chord_chord_ratio_initial = 0.3
    convergence = 0.005
    dCl = 1.4
    wet_area_ratio = 0.9
    sweep = 0

    iterate = Iterate()

    chord_ratios = iterate(chord_chord_ratio_initial, convergence,
                           angle, kf_data, dCl, wet_area_ratio, sweep)
    print(len(chord_ratios), round(chord_ratios[-1], 2))

    # for a in range(15,50,5):
    #     chord_ratios=iterate(chord_chord_ratio_initial,convergence,a,kf_data,dCl,wet_area_ratio,sweep)
    #     print(f"{a}: {round(chord_ratios[-1],2)} in {len(chord_ratios)} iterations")
