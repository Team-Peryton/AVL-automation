import subprocess as sp

from matplotlib.colors import same_color
from geometry import Plane,Section
import math
import numpy
import os
import shutil
from concurrent.futures import ThreadPoolExecutor
import matplotlib.pyplot as plt
from tqdm import tqdm

########################    INPUT & RUN    ##############################

### Loads inputs... 
def load_inputs(input_file):
    with open(input_file,'r') as file:
        lines=[line for line in file.readlines()]
    
    inputs={"input_plane":lines[1].split()[1],
            "wing_aerofoil":lines[3].split(": ")[1:][0],
            "elevator_aerofoil":lines[4].split(": ")[1:][0],
            "fin_aerofoil":lines[5].split(": ")[1:][0],
            "Xcg":lines[8].split()[1],
            "Ycg":lines[9].split()[1],
            "Zcg":lines[10].split()[1],
            "mass":float(lines[11].split()[1]),
            "Lt_upper":float(lines[14].split()[1]),
            "Lt_lower":float(lines[15].split()[1]),
            "St_h_upper":float(lines[16].split()[1]),
            "St_h_lower":float(lines[17].split()[1]),
            "St_v":float(lines[18].split()[1]),
            "steps":int(lines[19].split()[1]),
            "sm_ideal":float(lines[21].split()[1]),
            "tolerance":float(lines[22].split()[1]),
            "config":float(lines[24].split()[1]),
            "threads":int(lines[26].split()[1]),
            "calc_cg":False
            }

    if inputs["Lt_lower"]==0 or inputs["St_h_lower"]==0:
        print("Input non-zero lower bound.")
        exit()

    if inputs["Xcg"]=="-" and inputs["Ycg"]=="-" and inputs["Zcg"]=="-":
        inputs["calc_cg"]=True
    else:
        inputs["Xcg"]=float(inputs["Xcg"])
        inputs["Ycg"]=float(inputs["Ycg"])
        inputs["Zcg"]=float(inputs["Zcg"])

    return inputs

### Main run function.
def run(input_file):
    inputs=load_inputs(input_file)
    
    #   Generate planes
    mac,plane_geom, ARt=load_plane(inputs["input_plane"])    #   Loads input plane
    ref_plane=make_ref_plane(plane_geom,inputs["config"])    #   Strips input plane of elevator sections
    planes=generate_planes(ref_plane,inputs["St_h_lower"],
                            inputs["St_h_upper"],
                            inputs["Lt_lower"],
                            inputs["Lt_upper"],
                            inputs["steps"],ARt,mac,
                            inputs["elevator_aerofoil"],
                            inputs["Xcg"],
                            inputs["St_v"],
                            inputs["config"],
                            inputs["calc_cg"],
                            inputs["sm_ideal"]
                            )
    
    #   Run analysis   
    case=case_create(inputs["Xcg"],inputs["Ycg"],inputs["Zcg"],inputs["mass"])  #   Generates case file

    print("\nStability analysis...")
    tasks=[(case,plane) for plane in planes]    #   Generator for running all plane configs
    with ThreadPoolExecutor(max_workers=inputs["threads"]) as pool: #   Starts analysis on multiple threads
        list(tqdm(pool.map(run_analysis,tasks),total=len(tasks)))

    tasks=[plane for plane in planes]
    with ThreadPoolExecutor(max_workers=inputs["threads"]) as pool: #   Starts post processing on multiple threads
        pool.map(calc_SM,tasks)

    #   Plot SM against St,Lt & outputs possibel configurations.
    results(planes,inputs["tolerance"],inputs["calc_cg"])

########################    GEOMETRY    ##############################

### Reads input plane file for modification later
def load_plane(plane_file:str):
    """
    Loads input plane file & saves reference dims.
    """
    plane_geom=list()
    with open(plane_file,'r') as file:
        for line in file.readlines():
            if line=="\n" or line[0]=="#":  #   Ignores blank lines & comments
                continue
            else:
                plane_geom.append(line)
    ref_dims=[n for n in plane_geom][3]
    mac=float(ref_dims.split()[1])
    span=float(ref_dims.split()[2])
    ARw=span/mac
    ARt=ARw*2/3

    return mac, plane_geom, ARt

### Creates reference plane str exlucing elevator
def make_ref_plane(plane_geom:list,tail_config:str)->tuple:
    """
    Removes tail sections from input plane for modification
    """
    ref_plane=Plane("reference")
    ref_plane.tail_config=tail_config
    ref_plane_geom=ref_plane.make_reference(plane_geom)

    return tuple(ref_plane_geom)

### Creates modified plane object
def generate_planes(ref_plane:list,St_h_lower,St_h_upper,Lt_lower,Lt_upper,steps,ARt,mac,aerofoil,Xcg,St_v,config,calc_cg,sm_ideal)->object:
    planes=list()
    count=0
    a=0
    for St_h in numpy.linspace(St_h_lower,St_h_upper,steps):
        for Lt in numpy.linspace(Lt_lower,Lt_upper,steps):
            St_h=round(St_h,2)
            Lt=round(Lt,2)

            name=str(count)+"-"+str(St_h)+"St_h-"+str(Lt)+"Lt"  #   Creates plane name
            plane=Plane(name)   #   Initializes new plane
            plane.Lt=round(Lt,0)
            plane.St_h=St_h
            plane.mac=mac
            plane.sm_ideal=sm_ideal
            if calc_cg==False:
                plane.Xcg=Xcg
            plane.tail_config=config

            mod_geom=list(ref_plane)
            chord=round(math.sqrt(St_h/ARt)*1000,3)     #   Calculates HTP chord (mm)
            Zle=round((St_v*(1000**2))/(2*chord),3)         #   Calculates tip height (inverted v tail) (mm)
            span=round(math.sqrt(St_h*ARt)*1000,3)      #   Calculates HTP span (mm)
            
            plane.b_th=round(span,0)
            plane.b_tv=round(Zle,0)
            plane.c_t=round(chord,0)

            if config==0:
                root=Section(Lt,0,0,chord,10,-1,aerofoil)    #   Defines root section (object)
            elif config==1:
                root=Section(Lt,0,Zle,chord,10,-1,aerofoil)
            else:
                print("\nInvalid configuration selected.")
                exit()

            tip=Section(Lt,span/2,0,chord,10,-2,aerofoil)    #   Defines tip section (object)
            mod_str=root.create_input()+tip.create_input()  #   Combines 2 sections to insert into reference plane

            for index,line in enumerate(mod_geom):
                if line=="YES PLEASE\n":
                    mod_geom.pop(index)    #   Removes marker
                    mod_geom.insert(index,mod_str) #   Inserts modified sections

            plane.geom_file="generated planes/"+plane.name+".avl"
            with open(plane.geom_file,'w') as file:
                file.write("".join(mod_geom))
            count+=1

            planes.append(plane)

    print("Planes generated...")
    return(planes)

########################    ANALYSIS    ##############################

### Opens AVL
def AVL():
    """
    Initiates AVL.exe in subprocess.
    """
    return sp.Popen(['avl.exe'],
                stdin=sp.PIPE,
                stdout=sp.PIPE,
                stderr=sp.PIPE)

### Write command to AVL
def issueCmd(cmd: str):
    AVL().communicate(input=cmd.encode())

### Creates case file according to AVL format.
def case_create(Xcg,Ycg,Zcg,mass)->str:
    """
    Creates case string & writes to file
    """
    case_str="\n---------------------------------------------\n"
    case_str+="Run case  1:\n\n"
    case_str+="X_cg={0} Lunit\n".format(Xcg)
    case_str+="Y_cg={0} Lunit\n".format(Ycg)
    case_str+="Z_cg={0} Lunit\n".format(Zcg)
    case_str+="mass={0} kg\n".format(mass)

    path="cases/tail_case.txt"
    with open(path,'w') as file:    #   Saves case file
        file.write(case_str)
    return path 
    
### Runs analysis through AVL interface options & saves stability derivatives.
def run_analysis(tasks)->str:
    """
    Writes commands string to AVL
    """
    case,plane=tasks

    run="load {0}\n".format(plane.geom_file)    #   Load plane
    run+="case {0}\n".format(case)  #   Load case
    run+="oper\n x\n"   #   Run analysis
    run+="st\n" #   View stability derivatives
    
    plane.results_file="results/"+plane.name+".txt"
    run+=plane.results_file+"\n"    #   Saves results
    
    issueCmd(run)

########################    POST PROCESSING    ##############################

### Finds neutral point in results file & calculates SM.
def calc_SM(tasks):
    plane=tasks

    plane.calc_SM() #   Calculates static margin

def results(planes,tolerance,calc_cg):
    """
    Plots results
    """
    fig=plt.figure()
    ax=fig.add_subplot(projection='3d')

    x=[plane.St_h for plane in planes]
    y=[plane.Lt for plane in planes]

    if calc_cg==False:
        z=[plane.sm for plane in planes]

        ax.scatter(x,y,z,c=z)
        ax.set_xlabel("St_h (m^2)")
        ax.set_ylabel("Lt (m)")
        ax.set_zlabel("SM")

        if planes[0].tail_config==0:
            solutions=["\nPossible configurations:\nPlane ID:\tSM:\tLt (mm):\tSt_h (m^2):\n"]
            for plane in planes:
                if math.isclose(plane.sm,plane.sm_ideal,rel_tol=tolerance)==True:
                    solutions.append(plane.name.split("-")[0])
                    solutions.append(f"\t\t{str(plane.sm)}\t{str(plane.Lt)}\t\t{str(plane.St_h)}\n")
        elif planes[0].tail_config==1:
            solutions=["\nPossible configurations:\nPlane ID:\tSM:\tLt (mm)\tb (mm)\tc (mm)\tz (mm)\n"]
            for plane in planes:
                if math.isclose(plane.sm,plane.sm_ideal,rel_tol=tolerance)==True:
                    solutions.append(plane.name.split("-")[0])
                    solutions.append(f"\t\t{str(plane.sm)}\t{str(plane.Lt)}\t{str(plane.b_th)}\t{str(plane.c_t)}\t{str(plane.b_tv)}\n")
        if len(solutions)==1:
            print("\nNo ideal configurations possible. Consider changing limits.")
        else:
            solutions.append("\nConsider refining limits around possible configurations.\n")
            print("".join(solutions))

    elif calc_cg==True:
        z=[plane.np for plane in planes]

        ax.scatter(x,y,z,c=z)
        ax.set_xlabel("St_h (m^2)")
        ax.set_ylabel("Lt (m)")
        ax.set_zlabel(f"Xcg for SM={planes[0].sm_ideal}")

        if planes[0].tail_config==0:
            solutions=["\nPossible configurations:\nPlane ID:\tXcg:\tLt (mm):\tSt_h (m^2):\n"]
            for plane in planes:
                solutions.append(plane.name.split("-")[0])
                solutions.append(f"\t\t{str(plane.Xcg)}\t{str(plane.Lt)}\t\t{str(plane.St_h)}\n")
        elif planes[0].tail_config==1:
            solutions=["\nPossible configurations:\nPlane ID:\tXcg:\tLt (mm)\tb (mm)\tc (mm)\tz (mm)\n"]
            for plane in planes:
                solutions.append(plane.name.split("-")[0])
                solutions.append(f"\t\t{str(plane.Xcg)}\t{str(plane.Lt)}\t{str(plane.b_th)}\t{str(plane.c_t)}\t{str(plane.b_tv)}\n")
        print("".join(solutions))

    plt.show()
    pass


if __name__=="__main__":
    os.system('cls')
    input_file="TAIL_CONFIG.txt"

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

    run(input_file)
