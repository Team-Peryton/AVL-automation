import os
import shutil
import numpy as np
from matplotlib import pyplot as plt
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import copy

from geometry import Plane,Section
from avl import Case,avl_cmd

class AutoTail():
    def __init__(self,config_file:str): 
        path=os.path.abspath(os.getcwd())
        if os.path.isdir(path+"/results")==True:
            shutil.rmtree(path+"/results")
        if os.path.isdir(path+"/generated planes")==True:
            shutil.rmtree(path+"/generated planes")
        if os.path.isdir(path+"/cases")==True:
            shutil.rmtree(path+"/cases")
        os.mkdir(path+"/generated planes")
        os.mkdir(path+"/results")
        os.mkdir(path+"/cases")

        self.read_config(config_file)

        return None

    def read_config(self,file):
        with open(file,'r') as f:
            lines=[line for line in f.readlines()]
        lines=[line for line in lines if line[0]!="#" and line!="\n"]

        self.plane_file         = lines[0].split(": ")[1:][0].strip()
        self.wing_aerofoil      = lines[1].split(": ")[1:][0]
        self.elevator_aerofoil  = lines[2].split(": ")[1:][0]
        self.fin_aerofoil       = lines[3].split(": ")[1:][0]

        self.Xcg                = lines[4].split()[1]
        self.Ycg                = lines[5].split()[1]
        self.Zcg                = lines[6].split()[1]
        self.mass               = float(lines[7].split()[1])

        self.Lt_upper           = float(lines[8].split()[1])
        self.Lt_lower           = float(lines[9].split()[1])
        self.St_h_upper         = float(lines[10].split()[1])
        self.St_h_lower         = float(lines[11].split()[1])
        self.b_th               = lines[12].split()[1]
        self.St_v               = float(lines[13].split()[1])
        self.steps              = int(lines[14].split()[1])

        self.sm_ideal           = float(lines[15].split()[1])
        self.tolerance          = float(lines[16].split()[1])
        self.config             = float(lines[17].split()[1])
        self.threads            = int(lines[18].split()[1])
        
        if self.b_th!="-":
            self.b_th=float(self.b_th)

        if self.Lt_lower==0 or self.St_h_lower==0:
            print("Input non-zero lower bound.")
            exit()
        
        if self.Xcg=="-" and self.Ycg=="-" and self.Zcg=="-":
            self.calc_cg=True
        else:
            self.calc_cg=False

            self.Xcg=float(self.Xcg)
            self.Ycg=float(self.Ycg)
            self.Zcg=float(self.Zcg)

        return None

    def generate_planes(self):
        self.ref_plane=Plane(name="REF")
        self.ref_plane.read(self.plane_file)
        self.ref_plane.strip_section("Elevator")
        self.ref_plane.strip_surface("Fin")

        planes=[]

        St_h_range=np.linspace(self.St_h_lower,self.St_h_upper,self.steps)
        Lt_range=np.linspace(self.Lt_lower,self.Lt_upper,self.steps)

        mac=self.ref_plane.mac
        span=self.ref_plane.span
        ARw=self.ref_plane.ARw
        ARt=ARw*2/3

        count=0
        for St_h in St_h_range:
            for Lt in Lt_range:
                St_h=round(float(St_h),2)
                Lt=round(Lt,2)

                name=str(count)  #   Creates plane name
                plane=Plane(name)   #   Initializes new plane
                
                plane.Lt=round(Lt,0)
                plane.St_h=St_h
                plane.mac=mac
                plane.sm_ideal=self.sm_ideal

                if self.calc_cg==False:
                    plane.Xcg=self.Xcg
                plane.tail_config=self.config

                mod_geom=copy.copy(self.ref_plane.file_str)

                if self.b_th!="-":   #   if span constraint used:
                    chord=round((St_h*1000**2)/self.b_th,3)    #   Calculate chord based off span & area, not area & AR
                    span=self.b_th
                else:
                    chord=round(np.sqrt(St_h/ARt)*1000,3)     #   Calculates h chord based on area & AR
                    span=round(np.sqrt(St_h*ARt)*1000,3)      #   Calculates HTP span (mm)
                try:
                    Zle=round((self.St_v*(1000**2))/(2*chord),3)         #   Calculates tip height (inverted v tail) (mm)
                except:
                    exit()

                plane.b_th=round(span,0)
                plane.b_tv=round(Zle,0)
                plane.c_t=round(chord,0)

                if self.config==0 and self.b_th=="-":                
                    root=Section(Lt,0,0,chord,10,-1,self.elevator_aerofoil)    #   Defines root section (object)
                elif self.config==1:          
                    root=Section(Lt,0,Zle,chord,10,-1,self.elevator_aerofoil)
                else:
                    print("\nInvalid self.configuration selected.")
                    exit()

                tip=Section(Lt,span/2,0,chord,10,-2,self.elevator_aerofoil)    #   Defines tip section (object)
                mod_str=root.string()+tip.string()  #   Combines 2 sections to insert into reference plane

                for index,line in enumerate(mod_geom):
                    if line=="MARKER\n":
                        mod_geom.pop(index)    #   Removes marker
                        mod_geom.insert(index,mod_str) #   Inserts modified sections

                file_name=f"{plane.name}-{str(St_h)}{St_h}-{str(Lt)}Lt"
                plane.geom_file=f"generated planes/{file_name}.avl"
                
                with open(plane.geom_file,'w') as file:
                    file.write("".join(mod_geom))
                count+=1

                planes.append(plane)

        print("Planes generated...")
        self.planes=planes

        return(planes)

    def run(self):
        self.case=Case(self.Xcg,self.Ycg,self.Zcg,self.mass)
        self.case.write_stab_case()

        print("\nStability analysis...")
        tasks=[(self.case,plane) for plane in self.planes]
        with ThreadPoolExecutor(max_workers=self.threads) as pool: #   Starts analysis on multiple threads
            list(tqdm(pool.map(self.stab_analysis,tasks),total=len(tasks)))

        tasks=[plane for plane in self.planes]
        with ThreadPoolExecutor(max_workers=self.threads) as pool: #   Starts post processing on multiple threads
            pool.map(self.calc_SM,tasks)

    def stab_analysis(self,tasks):
        case,plane=tasks

        cmd_str="load {0}\n".format(plane.geom_file)    #   Load plane
        cmd_str+="case {0}\n".format(case)  #   Load case
        cmd_str+="oper\n x\n"   #   Run analysis
        cmd_str+="st\n" #   View stability derivatives
        
        plane.results_file="results/"+plane.name+".txt"
        cmd_str+=plane.results_file+"\n"    #   Saves results
        
        avl_cmd(cmd_str)

        return None

    def calc_SM(self,tasks):
        plane=tasks
        if self.calc_cg==False:
            plane.calc_SM()
        else:
            plane.calc_Xcg_ideal()

        return None

    def results(self):
        """
        Plots results
        """
        fig=plt.figure()
        ax=fig.add_subplot(projection='3d')

        x=[plane.St_h for plane in self.planes]
        y=[plane.Lt for plane in self.planes]

        if self.planes[0].tail_config==1:
                zz="z (mm):"
        else:
            zz=""

        if self.calc_cg==False:
            z=[plane.sm for plane in self.planes]

            ax.scatter(x,y,z,c=z)
            ax.set_xlabel("St_h (m^2)")
            ax.set_ylabel("Lt (m)")
            ax.set_zlabel("SM")
    
            solutions=[f"\nPossible configurations:\nPlane ID:\tSM:\tnp\tLt (mm):\tb (mm):\tc (mm):\t{zz}\n"]
            for plane in self.planes:
                if np.isclose(plane.sm,plane.sm_ideal,rtol=self.tolerance)==True:
                    solutions.append(plane.name.split("-")[0])
                    solutions.append(f"\t\t{str(plane.sm)}\t{str(plane.np)}\t{str(plane.Lt)}\t\t{str(plane.b_th)}\t{str(plane.c_t)}\t{str(plane.b_tv if zz!='' else '')}\n")
            
            if len(solutions)==1:
                print("\nNo ideal configurations possible. Consider changing limits.")
            else:
                solutions.append("\nConsider refining limits around possible configurations.\n")
                print("".join(solutions))

        else:
            z=[plane.np for plane in self.planes]

            ax.scatter(x,y,z,c=z)
            ax.set_xlabel("St_h (m^2)")
            ax.set_ylabel("Lt (m)")
            ax.set_zlabel(f"Xcg for SM={self.planes[0].sm_ideal}")

            solutions=["\nPossible configurations:\nPlane ID:\tXcg:\tnp (mm)\tLt (mm):\tb (mm):\tc (mm):\t{zz}\n"]
            for plane in self.planes:
                solutions.append(plane.name.split("-")[0])
                solutions.append(f"\t\t{str(plane.Xcg)}\t{str(plane.np)}\t{str(plane.Lt)}\t\t{str(plane.b_th)}\t{str(plane.c_t)}\t{str(plane.b_tv if zz!='' else '')}\n")

            print("".join(solutions))

        plt.show()

        return None

    def plot_plane(self,id):
        """
        Plots plane geometry with dimensions.
        """
        

        return None

if __name__=="__main__":
    tail=AutoTail("tail.config")
    tail.generate_planes()
    tail.run()
    tail.results()
