import os
import shutil
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import copy
import pandas as pd
from scipy import optimize

from .geometry import Plane, Section
from .aero import Case, avl_cmd


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

    def curve_fit_slice(self) -> list[plt.figure, plt.axes]:
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


class AutoTail():
    def __init__(self, config_file: str):
        
        self.path = os.path.split(config_file)[0]
        
        try:
            if os.path.isdir(self.path+"/results") == True:
                shutil.rmtree(self.path+"/results")
            if os.path.isdir(self.path+"/generated planes") == True:
                shutil.rmtree(self.path+"/generated planes")
            if os.path.isdir(self.path+"/cases") == True:
                shutil.rmtree(self.path+"/cases")
        except PermissionError:
            raise PermissionError("Close all results, geometry, case files")

        os.mkdir(self.path+"/generated planes")
        os.mkdir(self.path+"/results")
        os.mkdir(self.path+"/cases")
        
        if os.path.exists(f"{self.path}/avl.exe")==False:
            print("\u001b[31m[Error]\u001b[0m avl.exe not found.")
            exit()

        self.read_config(config_file)

        return None

    def read_config(self, file: str):
        """Reads tail.config data. 

        Args:
            file (str): tail.config file path.
        """

        with open(file, 'r') as f:
            lines = [line for line in f.readlines()]
        lines = [line for line in lines if line[0] != "#" and line != "\n"]

        if lines[0].strip() != "TAIL CONFIG":
            print(
                f"\u001b[31m[Error]\u001b[0m Wrong config file type ({file}).")
            exit()
        lines = lines[1:]

        self.plane_file = "".join(
            (self.path, '/', lines[0].split(": ")[1:][0].strip()))
        self.wing_aerofoil = "".join(
            (self.path, '/', lines[1].split(": ")[1:][0]))
        self.elevator_aerofoil = "".join(
            (self.path, '/', lines[2].split(": ")[1:][0]))
        self.fin_aerofoil = "".join(
            (self.path, '/', lines[3].split(": ")[1:][0]))

        self.Xcg = lines[4].split()[1]
        self.Ycg = lines[5].split()[1]
        self.Zcg = lines[6].split()[1]
        self.mass = float(lines[7].split()[1])

        self.Xt_upper = float(lines[8].split()[1])
        self.Xt_lower = float(lines[9].split()[1])
        self.St_h_upper = float(lines[10].split()[1])
        self.St_h_lower = float(lines[11].split()[1])
        self.Ct_v = float(lines[12].split()[1])
        self.steps = 7

        self.sm_ideal = float(lines[13].split()[1])
        self.tolerance = float(lines[14].split()[1])
        self.config = float(lines[15].split()[1])
        self.b_th = lines[16].split()[1]

        self.threads = int(lines[17].split()[1])

        if self.b_th != "NA":
            self.b_th = float(self.b_th)
        if self.b_th == 0:
            print(
                "\u001b[31m[Error]\u001b[0m Input non-zero horizontal span constraint. (NA to ignore)")
            exit()
        if self.b_th != "NA" and self.config == 0:
            print(
                "\u001b[33m[Warning]\u001b[0m Span constraint will be ignored for conventional tails.")

        if self.Xt_lower == 0 or self.St_h_lower == 0:
            print("\u001b[31m[Error]\u001b[0m Input non-zero lower bound.")
            exit()

        if self.Xcg == "NA" and self.Ycg == "NA" and self.Zcg == "NA":
            self.calc_cg = True
        else:
            self.calc_cg = False

            self.Xcg = float(self.Xcg)
            self.Ycg = float(self.Ycg)
            self.Zcg = float(self.Zcg)

        if self.config != 0 and self.config != 1:
            print(
                "\u001b[31m[Error]\u001b[0m Invalid tail configuration selected.")
            exit()

        self.St_h_range = np.linspace(
            self.St_h_lower, self.St_h_upper, self.steps)
        self.Xt_range = np.linspace(self.Xt_lower, self.Xt_upper, self.steps)

    def generate_planes(self):
        """Generates planes according to user tail limits. Generates and writes AVL geometry file.

        Returns:
            List[Plane]: List of Plane generated plane objects.
        """
        self.ref_plane = Plane(name="REF")
        self.ref_plane.read(self.plane_file)
        try:
            self.ref_plane.strip_section("Elevator")
        except KeyError:
            print(
                "\u001b[33m[Warning]\u001b[0m No section 'Elevator' found. Check if geometry of generated planes looks correct.")
        try:
            self.ref_plane.strip_surface("Fin")
        except KeyError:
            print(
                "\u001b[33m[Warning]\u001b[0m No surface 'Fin' found. Check if geometry of generated planes looks correct.")

        planes = []

        self.Sw = self.ref_plane.Sw
        self.mac = self.ref_plane.mac
        self.b_w = self.ref_plane.b_w
        ARw = self.ref_plane.ARw
        ARh = ARw*2/3
        Xw_root = self.ref_plane.Xw_root
        Cw_root = self.ref_plane.Cw_root

        count = 0
        for St_h in self.St_h_range:
            for Xt in self.Xt_range:
                St_h = float(St_h)
                Xt = Xt

                name = str(count)  # Creates plane name
                plane = Plane(name=name)  # Initializes new plane

                plane.Xt = Xt
                plane.Sw = self.Sw
                plane.Xw_root = Xw_root
                plane.Cw_root = Cw_root
                plane.St_h = St_h
                plane.ARh = ARh
                plane.mac = self.mac
                plane.b_w = self.b_w
                plane.sm_ideal = self.sm_ideal
                plane.tail_config = self.config
                plane.Ct_v = self.Ct_v

                if self.calc_cg == False:
                    plane.Xcg = self.Xcg

                mod_geom = copy.copy(self.ref_plane.file_str)

                if self.b_th != "NA" and self.config == 1:  # if span constraint used:
                    chord = St_h/self.b_th  # Calculate chord based off span & area, not area & AR
                    span = self.b_th
                else:
                    # Calculates h chord based on area & AR
                    chord = np.sqrt(St_h/plane.ARh)
                    # Calculates HTP span (Lunit)
                    span = np.sqrt(St_h*plane.ARh)

                plane.b_th = span
                plane.c_t = chord

                plane.Lt = (plane.Xt+plane.c_t*0.25) - \
                    (plane.Xw_root+0.25*plane.Cw_root)
                if plane.Lt <= 0:
                    print(
                        "\u001b[31m[Error]\u001b[0m Tail moment arm <=0. Increase Xt lower bound.")
                    exit()

                plane.St_v = plane.Ct_v*plane.Sw*plane.b_w/plane.Lt  # Vertical tail sizing

                # Calculates tip height (inverted v tail) (Lunit)
                Zle = (plane.St_v)/(2*chord)
                plane.theta = np.rad2deg(np.arctan(Zle/(span/2)))

                if self.config == 0:
                    # Defines root section (object)
                    root = Section(Xt, 0, 0, chord, 10, -1,
                                   self.elevator_aerofoil)
                elif self.config == 1:
                    root = Section(Xt, 0, Zle, chord, 10, -
                                   1, self.elevator_aerofoil)

                # Defines tip section (object)
                tip = Section(Xt, span/2, 0, chord, 10, -
                              2, self.elevator_aerofoil)
                # Combines 2 sections to insert into reference plane
                mod_str = str(root)+str(tip)

                for index, line in enumerate(mod_geom):
                    if line == "MARKER\n":
                        mod_geom.pop(index)  # Removes marker
                        # Inserts modified sections
                        mod_geom.insert(index, mod_str)

                file_name = f"{plane.name}-{str(round(St_h,2))}Sh-{str(round(plane.Lt,2))}Lt"
                plane.geom_file = f"{self.path}/generated planes/{file_name}.avl"

                with open(plane.geom_file, 'w') as file:
                    file.write("".join(mod_geom))
                count += 1

                planes.append(plane)

        print("[Info] Planes generated.")
        self.planes = planes

        return planes

    def run(self):
        """Runs AVL stability analysis. Multithreaded due to high io throughput.
        """
        self.case = Case(self.path,self.Xcg, self.Ycg, self.Zcg, self.mass)
        self.case.write_stab_case()

        tasks = [(self.case, plane) for plane in self.planes]
        # Starts analysis on multiple threads
        with ThreadPoolExecutor(max_workers=self.threads) as pool:
            list(tqdm(pool.map(self.stab_analysis, tasks),
                 total=len(tasks), desc="Stability analysis"))

        tasks = [plane for plane in self.planes]
        # Starts post processing on multiple threads
        with ThreadPoolExecutor(max_workers=self.threads) as pool:
            pool.map(self.calc_SM, tasks)

    def stab_analysis(self, tasks):
        """Creates AVL input string and executes AVL analysis.

        Args:
            tasks (List): Case and Plane to run [Case, Plane].

        Returns:
            None
        """
        case, plane = tasks

        cmd_str = "load {0}\n".format(plane.geom_file)  # Load plane
        cmd_str += "case {0}\n".format(case)  # Load case
        cmd_str += "oper\n x\n"  # Run analysis
        cmd_str += "st\n"  # View stability derivatives

        plane.results_file = f"{self.path}/results/"+plane.name+".txt"
        cmd_str += plane.results_file+"\n"  # Saves results

        avl_cmd(cmd_str, self.path)

    def calc_SM(self, tasks):
        """Calculates static margin for each plane.

        Args:
            tasks (Plane)
        """
        plane = tasks
        if self.calc_cg == False:
            plane.calc_SM()
        else:
            plane.calc_Xcg_ideal()

    def results(self, display=True):
        """Collates results.

        Args:
            display (bool, optional): _description_. Defaults to True.

        Returns:
            pd.DataFrame: solutions dataframe (redundant).
            CurveFit: CurveFit object for stable tail configs.
        """

        if self.calc_cg == False:

            columns = ["Plane ID", "Static Margin", "Xnp (Lunit)", "Xt (Lunit)", "Lt (Lunit)",
                       "Span (Lunit)", "Chord (Lunit)", "Angle (deg)", "Sh (Lunit^2)", "Sv (Lunit^2)", "ARh"]
            solutions = []
            for plane in self.planes:
                if np.isclose(plane.sm, plane.sm_ideal, atol=self.tolerance) == True:
                    solutions.append([
                        plane.name.split("-")[0],
                        plane.sm,
                        plane.np,
                        plane.Xt,
                        plane.Lt,
                        plane.b_th,
                        plane.c_t,
                        plane.theta,
                        plane.St_h,
                        plane.St_v,
                        plane.ARh,
                    ])

            solutions_df = pd.DataFrame(solutions, columns=columns)
            solutions_df = solutions_df.round(2)

            if self.config == 0:
                solutions_df = solutions_df[[
                    "Plane ID", "Static Margin", "Xnp (Lunit)", "Xt (Lunit)", "Lt (Lunit)", "Sh (Lunit^2)", "Sv (Lunit^2)", "ARh"]]

            curve_fit = CurveFit(self.planes, self.sm_ideal)
            if curve_fit.unstable==False:
                print(
                    "Consider refining limits around possible configurations.\n")
            print(self.planes[0].mac)
            print(f"np: {self.Xcg-(self.planes[0].mac*self.sm_ideal)}")

            if display == True:
                ##### Generated planes SM results (3D plot) #####
                x2, y2, z2 = curve_fit.curve_fit_surface()
                curve_fit.plot_surface(x2, y2, z2)

                if curve_fit.unstable==False:
                    Lt, St_h, St_v = curve_fit.curve_fit_slice()
                    curve_fit.plot_slice(Lt, St_h, St_v)


                plt.show()

            return solutions_df, curve_fit

        elif self.calc_cg == True:

            columns = ["Plane ID", "Xcg (Lunit)",
                       "np (Lunit)", "Static Margin"]
            solutions = []
            for plane in self.planes:
                solutions.append(
                    [plane.name.split("-")[0], plane.Xcg, plane.np, self.sm_ideal])

            solutions_df = pd.DataFrame(solutions, columns=columns)
            solutions_df = solutions_df.round(2)

            if display == True:
                fig = plt.figure()
                ax = fig.add_subplot(projection='3d')

                x = [plane.St_h for plane in self.planes]
                y = [plane.Lt for plane in self.planes]
                z = [plane.np for plane in self.planes]

                ax.scatter(x, y, z)

                ax.set_xlabel("St_h (Lunit^2)")
                ax.set_ylabel("Xt")
                ax.set_zlabel(f"Xcg for SM={self.planes[0].sm_ideal}")

                print("\n", solutions_df)

                plt.show()

            return solutions_df


if __name__ == "__main__":
    tail = AutoTail("projects/tail_MDDP_v0.config")
    tail.generate_planes()
    tail.run()
    tail.results()
