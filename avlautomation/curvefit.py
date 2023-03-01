import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from scipy import optimize

from avlautomation.geometry import Plane


class CurveFit():
    def __init__(self, planes: list[Plane], sm_ideal: float):
        self.planes = planes
        self.sm_ideal = sm_ideal

        self.Lts = np.array([plane.Lt for plane in self.planes])
        self.Xts = np.array([plane.Xt for plane in self.planes])
        self.SMs = np.array([plane.sm for plane in self.planes])
        self.St_hs = np.array([plane.St_h for plane in self.planes])

        if self.SMs.min() > sm_ideal or self.SMs.max() < sm_ideal:
            print(
                "\u001b[33m[Warning]\u001b[0m No stable configurations found. Consider changing limits.")

            self.unstable = True
        else:
            self.unstable = False

    def func(self, data: np.ndarray, a: float, b: float, c: float, d: float) -> np.ndarray:
        """Equation for 3D surface (applicable to Lt, Sh, SM datapoints)

        Args:
            data (np.ndarray): x,y data.
            a (float): Surface control parameter.
            b (float): Surface control parameter.
            c (float): Surface control parameter.
            d (float): Surface control parameter.

        Returns:
            np.ndarray: z values.
        """
        x = data[0]
        y = data[1]
        z = a*(x**b)*(y**c)+d
        return z

    def func_inv_const_z(self, x: np.ndarray, z: float, a: float, b: float, c: float, d: float) -> np.ndarray:
        """Solve 'func' for y given x and z.

        Args:
            x (np.ndarray): x data.
            z (float): z value at which to solve func.
            a (float): Surface control parameter.
            b (float): Surface control parameter.
            c (float): Surface control parameter.
            d (float): Surface control parameter.

        Returns:
            np.ndarray: y values.
        """
        y = np.exp((1/c)*np.log((z-d)/(a*x**b)))
        return y
    
    def Lt_to_Xt(self, x):
        return np.interp(x, self.Lts, self.Xts)

    def Xt_to_Lt(self, x):
        return np.interp(x, self.Xts, self.Lts)

    def curve_fit(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        """Scipy curve fit for func.

        Args:
            x (np.ndarray): Intiial datapoints.
            y (np.ndarray): Intiial datapoints.
            z (np.ndarray): Intiial datapoints.

        Returns:
            np.ndarray: Surface control parameters.
        """

        parameters, covariance = optimize.curve_fit(self.func, [x, y], z)

        return parameters

    def curve_fit_slice(self):
        """
        Slices surface fit to AVL datapoints.

        Returns:
            Lt: {np.ndarray} -- Tail moment arm
            St_h: {np.ndarray} -- Horizontal tail area
            St_v: {np.array} -- Vertical tail area

        """

        if self.unstable == True:
            print("\u001b[31m[Error]\u001b[0m SM ideal is out of range of analysis datapoints. Stable configurations are required to slice at SM ideal.")
            exit()

        Lts = self.Lts
        SMs = self.SMs
        St_hs = self.St_hs

        parameters = self.curve_fit(St_hs, Lts, SMs)

        St_h = np.linspace(St_hs.min(), St_hs.max(), 20)
        Lt = self.func_inv_const_z(St_h, self.sm_ideal, *parameters)
        St_v = self.planes[0].Ct_v*self.planes[0].Sw*self.planes[0].b_w/Lt

        return Lt, St_h, St_v

    def curve_fit_surface(self) -> list[np.ndarray]:
        """Fits surface to Lt, St_h, SM datapoints.

        Returns:
            list[np.ndarray]: x,y,z data (Lt, St_h, SM)
        """

        St_hs = [plane.St_h for plane in self.planes]
        Lts = self.Lts
        Sms = [plane.sm for plane in self.planes]

        parameters = self.curve_fit(self.Lts, self.St_hs, self.SMs)

        St_h_range = np.linspace(min(self.St_hs), max(self.St_hs), 20)
        Lt_range = np.linspace(min(self.Lts), max(self.Lts), 20)

        x2, y2 = np.meshgrid(Lt_range, St_h_range)
        z2 = self.func(np.array((x2, y2)), *parameters)

        return x2, y2, z2

    def plot_slice(self, Lt: np.ndarray, St_h: np.ndarray, St_v: np.ndarray):
        """
        Plots slices.

        Returns:
            fig: {plt.figure}
            ax1: {plt.axes}
            ax2: {plt.axes}
        """
        fig, ax1 = plt.subplots(figsize=(7, 7))

        ax1.plot(Lt, St_h, color='r', linestyle='-',
                 label=fr"Horizontal Tail"
                 )
        ax1.plot(Lt, St_v, color='b', linestyle='-',
                 label=fr"Vertical Tail"
                 )

        ax1.set_ylabel(r"$St$ ($Lunit^2$)")
        ax1.set_xlabel(r"$Lt$ ($Lunit$)")
        if max((max(St_h), max(St_v))) > 1000:
            ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        ax1.legend()

        ax2 = ax1.secondary_xaxis(
            'top', functions=(self.Lt_to_Xt, self.Xt_to_Lt))
        ax2.set_xlabel(r"Xt ($Lunits$)")

        ax1.set_title(fr"Tail Configurations with $SM={self.sm_ideal}$")
        # fig.tight_layout()

        return fig, ax1, ax2

    def plot_surface(self, x, y, z):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        ax.plot_surface(x, y, z, cmap=cm.viridis)
        ax.scatter(self.Lts, self.St_hs, self.SMs, color='k', depthshade=False)

        ax.set_xlabel("${St_h}$ (${Lunit^2}$)")
        ax.set_ylabel("${Lt}$ (${Lunit}$)")
        ax.set_zlabel("SM")

        fig.tight_layout()

        return plt

    def plot_surface_contour(self, x, y, z):

        fig, ax = plt.subplots()

        cs = ax.contour(x, y, z, 10, colors='k')
        ax.clabel(cs, cs.levels, inline=True, colors='k')

        ax.set_xlabel("${St_h}$ (${Lunit^2}$)")
        ax.set_ylabel("${Lt}$ (${Lunit}$)")

        fig.tight_layout()

        return plt