import copy
import os
from pathlib import Path
import shutil
from concurrent.futures import ThreadPoolExecutor
from uuid import uuid4

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from tqdm import tqdm

from avlautomation.aero import Case, avl_cmd
from avlautomation.config import TailConfig
from avlautomation.curvefit import CurveFit
from avlautomation.geometry import Plane, Section


class AutoTail():
    def __init__(self, config: TailConfig):
        self.config = config
        try:
            path = self.config.path.joinpath
            shutil.rmtree(path("results")) if os.path.isdir(path("results")) else None
            shutil.rmtree(path("generated planes")) if os.path.isdir(path("generated planes")) else None
            shutil.rmtree(path("cases")) if os.path.isdir(path("cases")) else None
        except PermissionError:
            raise PermissionError("Close all results, geometry, case files")

        list(map(os.mkdir, map(path, ["generated planes", "results", "cases"])))
        
        if os.path.exists(f"{self.config.path}/avl.exe")==False:
            raise Exception("\u001b[31m[Error]\u001b[0m avl.exe not found.")

        self.calc_cg = True if self.config.Cg_xyz == (0,0,0) else False
        self.St_h_range = np.linspace(self.config.St_h_lower, self.config.St_h_upper, self.config.steps)
        self.Xt_range = np.linspace(self.config.Xt_lower, self.config.Xt_upper, self.config.steps)


    def generate_plane(self, St_h, Xt):
        St_h = float(St_h)
        Xt = Xt

        plane = Plane(name=str(uuid4()))  # Initializes new plane

        plane.Xt = Xt
        plane.Sw = self.ref_plane.Sw
        plane.Xw_root = self.ref_plane.Xw_root
        plane.Cw_root = self.ref_plane.Cw_root
        plane.St_h = St_h
        plane.ARh = self.ref_plane.ARw*2/3
        plane.mac = self.ref_plane.mac
        plane.b_w = self.ref_plane.b_w
        plane.tail_config = self.config.tail_config_vtail
        plane.Ct_v = self.config.Ct_v

        if self.calc_cg == False:
            plane.Xcg = self.config.Cg_xyz[0]

        mod_geom = copy.copy(self.ref_plane.file_str)

        if self.config.horz_tail_span != 0 and self.config.tail_config_vtail:  # if span constraint used:
            chord = St_h/self.config.horz_tail_span  # Calculate chord based off span & area, not area & AR
            span = self.config.horz_tail_span
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
            raise Exception("\u001b[31m[Error]\u001b[0m Tail moment arm <=0. Increase Xt lower bound.")

        plane.St_v = plane.Ct_v*plane.Sw*plane.b_w/plane.Lt  # Vertical tail sizing

        # Calculates tip height (inverted v tail) (Lunit)
        Zle = (plane.St_v)/(2*chord)
        plane.theta = np.rad2deg(np.arctan(Zle/(span/2)))

        if self.config.tail_config_vtail:
            root = Section(Xt, 0, Zle, chord, 10, -1, self.config.elevator_aerofoil)
        else:
            root = Section(Xt, 0, 0, chord, 10, -1, self.config.elevator_aerofoil)


        # Defines tip section (object)
        tip = Section(Xt, span/2, 0, chord, 10, -
                        2, self.config.elevator_aerofoil)
        # Combines 2 sections to insert into reference plane
        mod_str = str(root)+str(tip)

        for index, line in enumerate(mod_geom):
            if line == "MARKER\n":
                mod_geom.pop(index)  # Removes marker
                # Inserts modified sections
                mod_geom.insert(index, mod_str)

        file_name = f"{plane.name}-{str(round(St_h,2))}Sh-{str(round(plane.Lt,2))}Lt"
        plane.geom_file = Path(f"{self.config.path}/generated planes/{file_name}.avl")

        with open(plane.geom_file, 'w') as file:
            file.write("".join(mod_geom))

        return plane


    def generate_planes(self):
        """Generates planes according to user tail limits. Generates and writes AVL geometry file.

        Returns:
            List[Plane]: List of Plane generated plane objects.
        """
        self.ref_plane = Plane(name="REF")
        self.ref_plane.read(self.config.input_plane)
        try:
            self.ref_plane.strip_section("Elevator")
        except KeyError:
            raise Exception("\u001b[33m[Warning]\u001b[0m No section 'Elevator' found. Check if geometry of generated planes looks correct.")
        try:
            self.ref_plane.strip_surface("Fin")
        except KeyError:
            raise Exception("\u001b[33m[Warning]\u001b[0m No surface 'Fin' found. Check if geometry of generated planes looks correct.")

        self.planes = [self.generate_plane(St_h, Xt) for Xt in self.Xt_range for St_h in self.St_h_range]

        print("[Info] Planes generated.")


    def run(self):
        """Runs AVL stability analysis. Multithreaded due to high io throughput.
        """
        self.case = Case(self.config.path,self.config.Cg_xyz[0], self.config.Cg_xyz[1], self.config.Cg_xyz[2], self.config.mass)
        self.case.write_stab_case()

        tasks = [(self.case, plane) for plane in self.planes]
        # Starts analysis on multiple threads
        with ThreadPoolExecutor(max_workers=self.config.threads) as pool:
            list(
                tqdm(pool.map(
                    self.stab_analysis, 
                    tasks
                    ),
                 total=len(tasks), 
                 desc="Stability analysis")
            )

        tasks = [plane for plane in self.planes]
        # Starts post processing on multiple threads
        with ThreadPoolExecutor(max_workers=self.config.threads) as pool:
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

        plane.results_file = f"{self.config.path}/results/"+plane.name+".txt"
        cmd_str += plane.results_file+"\n"  # Saves results
        
        avl_cmd(cmd_str, self.config.path)


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
                if np.isclose(plane.sm, self.config.sm_ideal, atol=self.config.tolerance) == True:
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

            print("\n", solutions_df)

            if self.config.tail_config_vtail == False:
                solutions_df = solutions_df[[
                    "Plane ID", 
                    "Static Margin", 
                    "Xnp (Lunit)", 
                    "Xt (Lunit)", 
                    "Lt (Lunit)", 
                    "Sh (Lunit^2)", 
                    "Sv (Lunit^2)", 
                    "ARh"
                ]]

            curve_fit = CurveFit(self.planes, self.config.sm_ideal)
            if curve_fit.unstable==False:
                print("Consider refining limits around possible configurations.\n")

            print(f"np: {self.config.Cg_xyz[0]-(self.planes[0].mac*self.config.sm_ideal)} Lunit")

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
                    [plane.name.split("-")[0], plane.Xcg, plane.np, self.config.sm_ideal])

            solutions_df = pd.DataFrame(solutions, columns=columns)
            solutions_df = solutions_df.round(2)

            if display == True:
                fig = plt.figure()
                ax:Axes3D = fig.add_subplot(projection='3d')

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
