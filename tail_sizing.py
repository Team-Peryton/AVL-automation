import subprocess as sp
from geometry import Plane,Section
import math
import numpy
import os
import shutil
from concurrent.futures import ThreadPoolExecutor
import matplotlib.pyplot as plt

########################    INPUT & RUN    ##############################

### Main run function.
def run(input_file):
    inputs=load_inputs(input_file)
    
    #   Generate planes
    mac,plane_geom, ARt=load_plane(inputs["input_plane"])    #   Loads input plane
    ref_plane=make_ref_plane(plane_geom,inputs["config"])    #   Strips input plane of elevator sections
    planes=generate_planes(ref_plane,inputs["St_lower"],inputs["St_upper"],inputs["Lt_lower"],inputs["Lt_upper"],inputs["steps"],ARt,mac,inputs["elevator_aerofoil"],inputs["Xcg"]) #   Generates configured planes
    
    #   Run analysis
    analysis_manager(inputs,planes)

    #   Plot SM against St,Lt & outputs possibel configurations.
    results(planes,inputs["SM_ideal"],inputs["tolerance"])

### Loads inputs... 
def load_inputs(input_file):
    with open(input_file,'r') as file:
        lines=[line for line in file.readlines()]
    
    inputs={"input_plane":lines[1].split()[1],
            "wing_aerofoil":lines[3].split(": ")[1:][0],
            "elevator_aerofoil":lines[4].split(": ")[1:][0],
            "fin_aerofoil":lines[5].split(": ")[1:][0],
            "Xcg":float(lines[8].split()[1]),
            "Ycg":float(lines[9].split()[1]),
            "Zcg":float(lines[10].split()[1]),
            "mass":float(lines[11].split()[1]),
            "Lt_upper":float(lines[14].split()[1]),
            "Lt_lower":float(lines[15].split()[1]),
            "St_upper":float(lines[16].split()[1]),
            "St_lower":float(lines[17].split()[1]),
            "steps":int(lines[18].split()[1]),
            "SM_ideal":float(lines[20].split()[1]),
            "tolerance":float(lines[21].split()[1]),
            "config":lines[23].split()[1],
            "threads":int(lines[25].split()[1])}

    if inputs["Lt_lower"]==0 or inputs["St_lower"]==0:
        print("Input non-zero lower bound.")
        exit()

    return inputs

########################    ANALYSIS    ##############################

### Analysis manager
def analysis_manager(inputs,planes):
    #   Perform analysis
    case=case_create(inputs["Xcg"],inputs["Ycg"],inputs["Zcg"],inputs["mass"])  #   Generates case file

    tasks=[(case,plane) for plane in planes]    #   Generator for running all plane configs
    with ThreadPoolExecutor(max_workers=inputs["threads"]) as pool: #   Starts analysis on multiple threads
        pool.map(run_analysis,tasks)
    tasks=[plane for plane in planes]
    with ThreadPoolExecutor(max_workers=inputs["threads"]) as pool: #   Starts post processing on multiple threads
        pool.map(calc_SM,tasks)
    
    print("Analysis complete...")

### Opens AVL
def AVL():
    return sp.Popen(['avl.exe'],
                stdin=sp.PIPE,
                stdout=sp.PIPE,
                stderr=sp.PIPE)

### Write command to AVL
def issueCmd(cmd: str):
    AVL().communicate(input=cmd.encode())

### Creates case file according to AVL format.
def case_create(Xcg,Ycg,Zcg,mass)->str:
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

def results(planes,sm_ideal,tolerance):
    fig=plt.figure()
    ax=fig.add_subplot(projection='3d')

    x=[plane.St for plane in planes]
    y=[plane.Lt for plane in planes]
    z=[plane.sm for plane in planes]

    ax.scatter(x,y,z,c=z)
    ax.set_xlabel("St (m^2)")
    ax.set_ylabel("Lt (m)")
    ax.set_zlabel("SM")

    solutions=["\nPossible configurations:\nSM:      Lt:    St:\n"]
    for plane in planes:
        if math.isclose(plane.sm,sm_ideal,rel_tol=tolerance)==True:
            solutions.append(str(plane.sm)+"     "+str(plane.Lt)+"    "+str(plane.St)+"\n")
    if len(solutions)==1:
        print("No ideal configurations possible. Consider changing limits.")
    else:
        solutions.append("\nConsider refining limits around possible configurations.")
        print("".join(solutions))

    return plt.show()

########################    GEOMETRY    ##############################

### Reads input plane file for modification later
def load_plane(plane_file:str):
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
    ref_plane=Plane("reference")
    ref_plane.tail_config=tail_config
    ref_plane_geom=ref_plane.make_reference(plane_geom)

    return tuple(ref_plane_geom)

### Creates modified plane object
def generate_planes(ref_plane:list,St_lower,St_upper,Lt_lower,Lt_upper,steps,ARt,mac,aerofoil,Xcg)->object:
    planes=list()
    count=0
    a=0
    for St in numpy.linspace(St_lower,St_upper,steps):
        for Lt in numpy.linspace(Lt_lower,Lt_upper,steps):
            St=round(St,2)
            Lt=round(Lt,2)

            name=str(count)+"-"+str(St)+"St-"+str(Lt)+"Lt"  #   Creates plane name
            plane=Plane(name)   #   Initializes new plane
            plane.Lt=Lt
            plane.St=St
            plane.mac=mac
            plane.Xcg=Xcg

            mod_geom=list(ref_plane)

            chord=round(math.sqrt(St/ARt)*1000,3)    #   Calculates HTP chord (mm)
            span=round(math.sqrt(St*ARt)*1000,3) #   Calculates HTP span (mm)

            root=Section(Lt,0,0,chord,10,0,aerofoil)    #   Defines root section (object)
            tip=Section(Lt,span/2,0,chord,10,0,aerofoil)    #   Defines tip section (object)
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

if __name__=="__main__":
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
