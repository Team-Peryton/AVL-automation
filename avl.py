import subprocess as sp
import numpy as np
from concurrent.futures import ThreadPoolExecutor

def avl_cmd(cmd_str:str)->None:
    avl_subprocess=sp.Popen(
        ['avl.exe'],
        stdin=sp.PIPE,
        stdout=sp.PIPE,
        stderr=sp.PIPE
    )

    avl_subprocess.communicate(input=cmd_str.encode())

    return None

class Case():
    def __init__(self,Xcg,Ycg,Zcg,mass,Ixx=None,Iyy=None,Izz=None,velocity=None,alpha=None,case_file=None,results_file=None,Cl=None,Cd=None,eigen=False,id=False,Clb=None,Clp=None,spiral=None):
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

    def write_aero_case(self):
        case_str =  "\n---------------------------------------------\n"
        case_str += f"Run case  1:\n\n"
        case_str += f"alpha -> alpha = {self.alpha}\n"
        case_str += f"X_cg={self.Xcg} Lunit\n"
        case_str += f"Y_cg={self.Ycg} Lunit\n"
        case_str += f"Z_cg={self.Zcg} Lunit\n"
        case_str += f"mass={self.mass} kg\n"
        case_str += f"Ixx={self.Ixx} kg-m^2\n"
        case_str += f"Iyy={self.Iyy} kg-m^2\n"
        case_str += f"Izz={self.Izz} kg-m^2\n"
        case_str += f"velocity={self.velocity} m/s\n"
        case_str += "density=1.225 kg-m^3\n"
        case_str += "grav.acc.=0.98 m/s^2\n"

        path=f"cases/{str(self.alpha)}deg.case"

        with open(path,'w') as f:
            f.write(case_str)

        self.case_file=path

        return None

    def write_stab_case(self)->None:
        """
        Creates case string & writes to file
        """
        case_str =  "\n---------------------------------------------\n"
        case_str += "Run case  1:\n\n"
        case_str += f"X_cg={self.Xcg} Lunit\n"
        case_str += f"Y_cg={self.Ycg} Lunit\n"
        case_str += f"Z_cg={self.Zcg} Lunit\n"
        case_str += f"mass={self.mass} kg\n"

        path="cases/tail.case"
        with open(path,'w') as file:    #   Saves case file
            file.write(case_str)

        self.case_file=path

        return None 

class Aero():
    def __init__(self,config_file:str):
        self.read_config(config_file)

        alpha_range=np.linspace(
            self.alpha0,
            self.alpha1,
            int(1+(self.alpha1-self.alpha0)/self.increment)
        )

        self.cases=[]
        for alpha in alpha_range:
            self.cases.append(Case(
                self.Xcg,
                self.Ycg,
                self.Zcg,
                self.Ixx,
                self.Iyy,
                self.mass,
                self.velocity,
                alpha
            ))
            
        #   Writes case files & adds filepath to case obj
        with ThreadPoolExecutor(max_workers=self.inputs["threads"]) as pool:
            pool.map(self.create_self.cases,self.cases)

        
    def read_config(self,file:str)->None:
        """
        Reads aero config file.
        """
        with open(file,'r') as f:
            lines=f.readlines()

        try:
            self.mass       = float(lines[1].split()[1])
            self.Xcg        = float(lines[2].split()[1])
            self.Ycg        = float(lines[3].split()[1])
            self.Zcg        = float(lines[4].split()[1])
            self.Ixx        = float(lines[5].split()[1])
            self.Iyy        = float(lines[6].split()[1])
            self.Izz        = float(lines[7].split()[1])
            self.velocity   = float(lines[9].split()[1])
            self.alpha0     = float(lines[11].split()[1])
            self.alpha1     = float(lines[12].split()[1])
            self.increment  = float(lines[13].split()[1])
            self.threads    = int(lines[15].split()[1])
            self.units      = lines[16].split()[1]
        except IndexError:
            print("Parameters must have a value assigned. (AERO_CONFIG.txt)")
            exit()

        return None

    def create_cases(self,case):
        case.write_file()

        return None


    def run(self,plane):
        

        return results_file

    def run_(self,tasks):
        case,plane=tasks

        cmd_str=f"load{plane.geom_file}\n"
        cmd_str+=f"case {self.case_file}\n"
        cmd_str+=f"mass {self.config['units']}\n"
        cmd_str+="oper\n0\nv\n\nx\n"

        results_file=f"results/{plane.name.split('-')[0]}-{str(case.alpha)}deg"
        if case.eigen==True:
            cmd_str+="\nmode\nN\nW\n"
            case.results_file+=".eig"
        else:
            run+="st\n"
            case.results_file+=".aero"
        
        avl_cmd(cmd_str)

        return None


class Stability():
    def __init__(self):
        pass