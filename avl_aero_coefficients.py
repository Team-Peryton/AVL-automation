import subprocess as sp
from numpy import linspace
import os
import shutil
from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor

class Case():
    def __init__(self,Xcg,Ycg,Zcg,mass,velocity,alpha,case_file=None,results_file=None,Cl=None,Cd=None):
        self.Xcg=Xcg
        self.Ycg=Ycg
        self.Zcg=Zcg
        self.mass=mass
        self.velocity=velocity
        self.alpha=alpha
        self.case_file=case_file
        self.results_file=results_file
        self.Cl=Cl
        self.Cd=Cd

########################    INPUT & RUN    ##############################
class Aero():
    def __init__(self,input_file,plane_file,cases=None,inputs=None,polars=None):
        self.input_file=input_file
        self.plane_file=plane_file
        self.cases=cases
        self.inputs=inputs
        self.polars=polars

        path=os.path.abspath(os.getcwd())
        if os.path.isdir(path+"/cases")==True:
            shutil.rmtree(path+"/cases")
        os.mkdir(path+"/cases")
        if os.path.isdir(path+"/results")==True:
            shutil.rmtree(path+"/results")
        os.mkdir(path+"/results")

    ### Main run function.
    def run(self):
        
        inputs=self.load_inputs()

        #   Creates alpha range
        alpha_range=linspace(
                            inputs["alpha0"],
                            inputs["alpha1"],
                            int(1+(inputs["alpha1"]-inputs["alpha0"])/inputs["increment"])
                            )
        #   Initiates class objects
        self.cases=[Case(inputs["Xcg"],inputs["Ycg"],inputs["Zcg"],inputs["mass"],inputs["velocity"],alpha) for alpha in alpha_range]
        
        #   Creates cases
        tasks=[case for case in self.cases]
        with ThreadPoolExecutor(max_workers=inputs["threads"]) as pool:
            pool.map(self.case_create,tasks)
        #with ProcessPoolExecutor(max_workers=8) as pool:
        #    pool.map(self.case_create,tasks)

        #   Runs analysis
        tasks=[(case,self.plane_file) for case in self.cases]
        with ThreadPoolExecutor(max_workers=inputs["threads"]) as pool:
            pool.map(self.run_analysis,tasks)
        #with ProcessPoolExecutor(max_workers=8) as pool:
        #    pool.map(self.run_analysis,tasks)
        self.results()

    ### Loads inputs... 
    def load_inputs(self):
        with open(self.input_file,'r') as file:
            lines=file.readlines()
        
        inputs={"Xcg":float(lines[1].split()[1]),
                "Ycg":float(lines[2].split()[1]),
                "Zcg":float(lines[3].split()[1]),
                "mass":float(lines[4].split()[1]),
                "velocity":float(lines[5].split()[1]),
                "alpha0":float(lines[7].split()[1]),
                "alpha1":float(lines[8].split()[1]),
                "increment":float(lines[9].split()[1]),
                "threads":int(lines[11].split()[1])}

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
        case_str+="velocity={0} Lunit/Tunit\n".format(case.velocity)
        
        path="cases/case "+str(case.alpha)+"deg.txt"

        with open(path,'w') as file:    #   Saves case file
            file.write(case_str)
        case.case_file=path
        
    ### Runs analysis through AVL interface options & saves stability derivatives.
    def run_analysis(self,tasks)->str:
        case,plane=tasks

        run="load {0}\n".format(plane)    #   Load plane
        run+="case {0}\n".format(case.case_file)  #   Load case
        run+="oper\n o\n v\n\n x\n"   #   Run analysis
        run+="ft\n" #   View stability derivatives
        
        case.results_file="results/"+str(case.alpha)+"deg.txt"

        run+=case.results_file+"\n"    #   Saves results
        
        self.issueCmd(run)

    ########################    RESULTS    ##############################

    def results(self):
        for case in self.cases:
            with open(case.results_file,'r') as file:
                lines=file.readlines()

                case.Cl=float(lines[23].split()[2])
                case.Cd=float(lines[24].split()[2])

        polars=[(case.alpha,case.Cl,case.Cd) for case in self.cases]
        self.polars=tuple(polars)