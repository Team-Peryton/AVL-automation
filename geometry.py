import subprocess as sp
from matplotlib.pyplot import clabel
from numpy import linspace
import os
import shutil
from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor
import pandas as pd

class Case():
    def __init__(self,Xcg,Ycg,Zcg,Ixx,Iyy,Izz,mass,velocity,alpha,case_file=None,results_file=None,Cl=None,Cd=None,eigen=False,id=False,Clb=None,Clp=None,spiral=None):
        self.Xcg=Xcg
        self.Ycg=Ycg
        self.Zcg=Zcg
        self.Ixx=Ixx
        self.Iyy=Iyy
        self.Izz=Izz
        self.mass=mass
        self.velocity=velocity
        self.alpha=alpha
        self.case_file=case_file
        self.Cl=Cl
        self.Cd=Cd
        self.results_file=results_file
        self.eigen=eigen
        self.id=id
        self.Clb=Clb
        self.Clp=Clp
        self.spiral=spiral

########################    INPUT & RUN    ##############################
class Aero():
    def __init__(self,input_file,polars=None):
        self.input_file=input_file
        self.inputs=self.load_inputs()

        path=os.path.abspath(os.getcwd())
        if os.path.isdir(path+"/cases")==True:
            shutil.rmtree(path+"/cases")
        os.mkdir(path+"/cases")
        if os.path.isdir(path+"/results")==True:
            shutil.rmtree(path+"/results")
        os.mkdir(path+"/results")

    def initialize_cases(self):
        """
        Creates case strings for each alpha and writes to file. 
        """       
        #   Creates alpha range
        alpha_range=linspace(
                            self.inputs["alpha0"],
                            self.inputs["alpha1"],
                            int(1+(self.inputs["alpha1"]-self.inputs["alpha0"])/self.inputs["increment"])
                            )
        #   Initiates class objects
        cases=[Case(self.inputs["Xcg"],self.inputs["Ycg"],self.inputs["Zcg"],self.inputs["Ixx"],self.inputs["Iyy"],self.inputs["Izz"],self.inputs["mass"],self.inputs["velocity"],alpha) for alpha in alpha_range]
        
        #   Writes case files & adds filepath to case obj
        tasks=[case for case in cases]
        with ThreadPoolExecutor(max_workers=self.inputs["threads"]) as pool:
            pool.map(self.case_create,tasks)

        return cases

    ### Loads inputs... 
    def load_inputs(self):
        with open(self.input_file,'r') as file:
            lines=file.readlines()
        
        try:
            inputs={"mass":float(lines[1].split()[1]),
                    "Xcg":float(lines[2].split()[1]),
                    "Ycg":float(lines[3].split()[1]),
                    "Zcg":float(lines[4].split()[1]),
                    "Ixx":float(lines[5].split()[1]),
                    "Iyy":float(lines[6].split()[1]),
                    "Izz":float(lines[7].split()[1]),
                    "velocity":float(lines[9].split()[1]),
                    "alpha0":float(lines[11].split()[1]),
                    "alpha1":float(lines[12].split()[1]),
                    "increment":float(lines[13].split()[1]),
                    "threads":int(lines[15].split()[1]),
                    "units":lines[16].split()[1]
                    }
        except IndexError:
            print("Parameters must have a value assigned. (AERO_CONFIG.txt)")
            exit()

        return(inputs)
    ########################    ANALYSIS    ##############################

    ### Opens AVL
    def AVL(self):
        return sp.Popen(['avl.exe'],
                    stdin=sp.PIPE,
                    stdout=sp.PIPE,
                    stderr=sp.PIPE)

    ### Write command to AVL
    def issueCmd(self,cmd: str):
        self.AVL().communicate(input=cmd.encode())

    ### Creates case file according to AVL format.
    def case_create(self,tasks)->str:
        case=tasks

        case_str="\n---------------------------------------------\n"
        case_str+="Run case  1:\n\n"
        case_str+="alpha -> alpha = {0}\n".format(case.alpha)
        case_str+="X_cg={0} Lunit\n".format(case.Xcg)
        case_str+="Y_cg={0} Lunit\n".format(case.Ycg)
        case_str+="Z_cg={0} Lunit\n".format(case.Zcg)
        case_str+="mass={0} kg\n".format(case.mass)
        case_str+="Ixx={0} kg-m^2\n".format(case.Ixx)
        case_str+="Iyy={0} kg-m^2\n".format(case.Iyy)
        case_str+="Izz={0} kg-m^2\n".format(case.Izz)
        case_str+="velocity={0} m/s\n".format(case.velocity)
        case_str+="density=1.225 kg-m^3\n"
        case_str+="grav.acc.=0.98 m/s^2\n"
        
        path="cases/"+str(case.alpha)+"deg.txt"

        with open(path,'w') as file:    #   Saves case file
            file.write(case_str)
        case.case_file=path
        
    ### Runs analysis through AVL interface options & saves stability derivatives.
    def analysis(self,plane,case)->str:
        """
        Runs aero analysis.
        """     
        run="load {0}\n".format(plane.geom_file)    #   Load plane
        run+="case {0}\n".format(case.case_file)  #   Load case
        run+="mass {0}\n".format(self.inputs["units"])
        run+="oper\n o\n v\n\n x\n"   #   Run analysis

        if case.eigen==True:
            run+="\nmode\n N\n W\n"   #   Run eigenvalue analysis
        else:
            run+="st\n" #   View stability derivatives
        case.results_file="".join(["results/",plane.name.split("-")[0],"-",str(case.alpha),"deg.aero" if case.eigen==False else "deg.eig"])       
        run+=case.results_file+"\n"    #   Saves results
        
        self.issueCmd(run)

        pass

    ########################    RESULTS    ##############################

    @staticmethod
    def read_aero(case):
        with open(case.results_file,'r') as file:
            lines=file.readlines()

            Cl=float(lines[23].split()[2])
            Cd=float(lines[24].split()[2])
            Clb=float(lines[38].split()[8])
            Clp=float(lines[46].split()[5])
            try:
                spiral=float(lines[52].split()[6])
            except IndexError:
                print("! No vertical stabilisation !")
                spiral=0
                pass

        return Cl,Cd,Clb,Clp,spiral

    @staticmethod
    def read_eigen(plane,case):
        with open(case.results_file,'r') as file:
            lines=file.readlines()

            try:
                modes={"dutch":tuple(map(float,lines[3].split()[1:])),
                        #"-dutch":tuple(map(float,lines[4].split()[1:])),
                        "roll":tuple(map(float,lines[5].split()[1:])),
                        #"short":tuple(map(float,lines[6].split()[1:])),
                        #"-short":tuple(map(float,lines[7].split()[1:])),
                        #"lateral":tuple(map(float,lines[8].split()[1:])),
                        #"phugoid":tuple(map(float,lines[9].split()[1:])),
                        #"-phugoid":tuple(map(float,lines[10].split()[1:])),
                        }
            except IndexError as e:
                print(f"Eigenmode analysis/read failed: Case {case.results_file}")
                print(f"\n{e}")
                exit()
        modes_df=pd.DataFrame.from_dict(modes,orient='index',columns=["Damping Ratio","Frequency"])
        #print(modes_df)
        plane.eigen_modes=modes

        pass
