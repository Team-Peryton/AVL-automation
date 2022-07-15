from cv2 import split
from avl import Aero
from geometry import Plane,Section
from concurrent.futures import ThreadPoolExecutor
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os
import shutil
from tqdm import tqdm
import copy

class Dihedral():
    def __init__(self,dihedral_config_file:str,aero_config_file:str):
        path=os.path.abspath(os.getcwd())
        try:
            if os.path.isdir(path+"/results")==True:
                shutil.rmtree(path+"/results")
            if os.path.isdir(path+"/generated planes")==True:
                shutil.rmtree(path+"/generated planes")
            if os.path.isdir(path+"/cases")==True:
                shutil.rmtree(path+"/cases")
        except PermissionError:
            print("! Close all results/geometry/case files !")
            exit()
        os.mkdir(path+"/generated planes")
        os.mkdir(path+"/results")
        os.mkdir(path+"/cases")

        self.read_config(dihedral_config_file)
        self.aero_config_file=aero_config_file

        return None

    def read_config(self,config_file):
        str_to_bool=lambda x:True if (x=="Y") else False

        with open(config_file,'r') as f:
            lines=f.readlines()
        lines=[line for line in lines if line[0]!="#" and line!="\n"]
    
        self.plane_file         = lines[0].split(": ")[1:][0].strip()
        self.wing_aerofoil      = lines[1].split(": ")[1:][0]
        self.elevator_aerofoil  = lines[2].split(": ")[1:][0]
        self.fin_aerofoil       = lines[3].split(": ")[1:][0]
        self.angle_min          = float(lines[4].split()[1])
        self.angle_max          = float(lines[5].split()[1])
        self.increment          = int(lines[6].split()[1])
        self.span_loc           = float(lines[7].split()[1])
        self.threads            = float(lines[8].split()[1])
        self.show_geom_plt      = str_to_bool(lines[9].split(": ")[1][0])

        return None

    def generate_planes(self):
        self.ref_plane=Plane(name="REF")
        self.ref_plane.read(self.plane_file)
        self.ref_plane.strip_section("Main Wing")
        self.ref_plane.strip_surface("Fin")

        planes=[]

        mac=self.ref_plane.mac
        span=self.ref_plane.span
        hspan=span/2

        count=0
        theta_range=np.linspace(
            self.angle_min,
            self.angle_max,
            int(1+(self.angle_max-self.angle_min)/self.increment)
        )
        for theta in theta_range:
            name=str(count)

            plane=Plane(name)
            plane.dihedral_angle=theta
            plane.dihedral_split=self.span_loc
            
            split_loc=hspan*self.span_loc/100
            plane.dihedral_splitY=split_loc
            plane.span=span

            mod_geom=copy.copy(self.ref_plane.file_str)

            Zle=round((hspan-split_loc)*np.sin(np.radians(theta)),3)    #   Calculates tip Z due to dihedral angle
            plane.tipZ=Zle
            Yle=round((hspan-split_loc)*np.cos(np.radians(theta))+split_loc,3)  #   Calcualtes tip Y due to dihedral angle
            plane.tipY=Yle

            root=Section(self.ref_plane.Xle,0,0,mac,int(split_loc*0.02),-2,self.elevator_aerofoil)
            split=Section(self.ref_plane.Xle,split_loc,0,mac,int((np.sqrt(Yle**2+Zle**2)-split_loc)*0.02),-1,self.elevator_aerofoil)
            tip=Section(self.ref_plane.Xle,Yle,Zle,mac,0,0,self.elevator_aerofoil)   #   Creates tip section based off tip geometry
            mod_str=root.string()  #   Gets section string in avl format
            mod_str+=split.string()
            mod_str+=tip.string()

            for index,line in enumerate(mod_geom):
                if line=="MARKER\n":    #   Finds marker
                    mod_geom.pop(index) #   Removes marker
                    mod_geom.insert(index,mod_str)  #   Inserts modified sections
                
            file_name=name=f"{plane.name}-{theta}deg-{self.span_loc}%"
            plane.geom_file=f"generated planes/{file_name}.avl"
            with open(plane.geom_file,'w') as file:
                file.write("".join(mod_geom))
            count+=1

            planes.append(plane)
        
        self.planes=planes

        return planes

    def run(self):
        aero=Aero(self.aero_config_file)    #   initialises aero analysis, reads config file.

        tasks=[(plane,aero) for plane in self.planes]

        print("Polar analysis...")
        with ThreadPoolExecutor(max_workers=self.threads) as pool:
            list(tqdm(pool.map(self.aero_analysis,tasks),total=len(tasks)))

        return None

    def aero_analysis(self,tasks):
        plane,aero=tasks
        aero.run(plane)

        return None

    def plot_polars(self):
        fig,(ax1,ax2,ax3)=plt.subplots(ncols=3,figsize=(12,4))
        
        dihedral_angles=[plane.dihedral_angle for plane in self.planes]
        Cl_0=self.planes[0].polars['Cl'].iloc[-1]
        Cd_0=self.planes[0].polars['Cd'].iloc[-1]
        Cl_delta=[100*(plane.polars['Cl'].iloc[-1]-Cl_0)/Cl_0 for plane in self.planes]
        Cd_delta=[100*(plane.polars['Cd'].iloc[-1]-Cd_0)/Cd_0 for plane in self.planes]

        ax1.plot(dihedral_angles,Cl_delta,label="\u0394Cl",color='r')
        ax1.plot(dihedral_angles,Cd_delta,label="\u0394Cd",color='b')
        ax1.set_xlabel("Dihedral Angles (\u00B0)")
        ax1.set_ylabel(f"\u0394 (%) @ {self.planes[0].polars['Alpha (deg)'].iloc[-1]}\u00B0")
        ax1.legend(loc='upper left')
        ax1.set_title("Aero Coeffients")

        Clb_0=self.planes[0].polars['Clb'].iloc[0]
        Clp_0=self.planes[0].polars['Clp'].iloc[0]
        Clb_delta=[100*(plane.polars['Clb'].iloc[0]-Clb_0)/Clb_0 for plane in self.planes]
        Clp_delta=[100*(plane.polars['Clp'].iloc[0]-Clp_0)/Clp_0 for plane in self.planes]
        
        ax2.plot(dihedral_angles,Clb_delta,label="Dihedral Effect Derivative",color='r')
        ax2.plot(dihedral_angles,Clp_delta,label="Roll rate derivative",color='b')
        ax2.set_xlabel("Dihedral Angles (\u00B0)")
        ax2.set_ylabel("\u0394 (%)")
        ax2.legend()
        ax2.set_title("Stability Derivatives")

        spiral=[plane.polars['spiral'].iloc[0] for plane in self.planes]

        ax3.plot(dihedral_angles,spiral,color='k')   
        ax3.set_xlabel("Dihedral Angles (\u00B0)")
        ax3.set_title("Spiral Stability (>1 = stable)")

        fig.tight_layout()

        return plt

    def eigenvalues(self): 
        dihedral_angles=[plane.dihedral_angle for plane in self.planes]
        roll_0=self.planes[0].modes_modes["roll"][0]
        dutch_0=self.planes[0].modes_modes["dutch"][0]
        roll_delta=[100*(plane.modes_modes["roll"][0]-roll_0)/roll_0 for plane in self.planes]
        dutch_delta=[100*(plane.modes_modes["dutch"][0]-dutch_0)/dutch_0 for plane in self.planes]

        plt.figure(figsize=(4,4))
        plt.title("Eigenmode Damping")
        plt.xlabel(f"Dihedral Angle (\u00B0)\nSplit Location={[plane.dihedral_split for plane in self.planes][0]}% of Span")
        plt.ylabel("\u0394 Damping (%)")

        plt.plot(dihedral_angles,roll_delta,color='r',label="Roll")
        plt.plot(dihedral_angles,dutch_delta,color='b',label="Dutch Roll")
        plt.legend()
        plt.tight_layout()
        
        return plt

    def geom_plot(self):
        plt.figure()
        plt.title("Spanwise Geometry Plot")
        plt.xlabel("Y (mm)")
        plt.ylabel("Z (mm)")  
        plt.xlim(0,max([plane.tipY for plane in self.planes]))
        plt.ylim(0,max([plane.tipY for plane in self.planes]))
        for plane in self.planes:
            plt.plot([0,plane.dihedral_splitY,plane.tipY],[0,0,plane.tipZ])
        
        return plt


if __name__=="__main__":
    dihedral=Dihedral('dihedral.config','aero.config')
    dihedral.generate_planes()
    dihedral.run()
