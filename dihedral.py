from typing import Counter
from aero import Aero
from geometry import Plane,Section
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor,as_completed
from multiprocessing import freeze_support
import matplotlib.pyplot as plt
import pandas as pd
import numpy
import math
import os
import shutil
import time
from tqdm import tqdm

def load_inputs(input_file):
    with open(input_file,'r') as f:
        lines=f.readlines()
    
    inputs={"input_plane":lines[1].split()[1],
            "aero_config":lines[2].split()[1],
            "wing_aerofoil":lines[4].split(": ")[1:][0],
            "elevator_aerofoil":lines[5].split(": ")[1:][0],
            "fin_aerofoil":lines[6].split(": ")[1:][0],
            "angle_min":float(lines[9].split()[1]),
            "angle_max":float(lines[10].split()[1]),
            "increment":int(lines[11].split()[1]),
            "span_loc":float(lines[13].split()[1]),
            "threads":float(lines[15].split()[1]),
            "show_geom_plt":lines[16].split(": ")[1][0]
            }

    return inputs

def load_plane(input_plane:str):
    with open(input_plane,'r') as f:
        plane_geom=[line for line in f.readlines()  if line!="\n" and line[0]!="#"]
    
    ref_dims=[line for line in plane_geom][3]
    mac=float(ref_dims.split()[1])
    span=float(ref_dims.split()[2])

    return plane_geom, mac, span

def run(input_file):
    inputs=load_inputs(input_file)
    plane_geom,mac,span=load_plane(inputs["input_plane"])
    ref_plane_geom,ref_plane=make_ref_plane(plane_geom,mac,span,inputs["span_loc"],inputs["wing_aerofoil"])

    analysis=Aero(inputs["aero_config"])
    planes=generate_planes(ref_plane_geom,inputs["angle_min"],inputs["angle_max"],inputs["increment"],inputs["span_loc"],span,ref_plane,mac,inputs["wing_aerofoil"],analysis)
   
    tasks=[]

    for plane in planes: #list
        for case in plane.cases: #list
            tasks.append([plane,case,analysis]) # one list

    print("Polar analysis...")
    with ThreadPoolExecutor(max_workers=inputs["threads"]) as pool:
        list(tqdm(pool.map(run_analysis,tasks),total=len(tasks)))

    print("Reading polar results...")
    polars(planes)
    
    plot_polars(planes)

    tasks=[]
    for plane in planes:
        for case in plane.cases:
            if case.alpha==0:
                case.eigen=True
                tasks.append([plane,case,analysis])

    print("\nEigenmode analysis...")
    with ThreadPoolExecutor(max_workers=inputs["threads"]) as pool:
        list(tqdm(pool.map(run_analysis,tasks),total=len(tasks)))
    
    print("Reading eigenmode results...\n")
    eigenvalues(planes,analysis)
    if inputs["show_geom_plt"]=="Y":
        geom_plot(planes)

    plt.show()

    pass

def make_ref_plane(plane_geom:list,mac,span,span_loc,wing_aerofoil)->tuple:
    ref_plane=Plane("reference",mac=mac)    #   Creates reference plane object
    ref_plane.dihedral_splitY=(span/2)*(span_loc/100)   #   Convert from %
    ref_plane_geom=ref_plane.make_dihedral_ref(plane_geom)

    return tuple(ref_plane_geom), ref_plane

def generate_planes(ref_plane_geom:list,angle_min,angle_max,increment,span_loc,span,ref_plane,mac,aerofoil,analysis):
    """
    Takes plane modification info & writes to string. String inserted to reference plane geometry & saved to a file.
    """    
    planes=[]
    count=0
    hspan=span/2
  
    #   Generates range of angles from min, max, and increment
    for angle in numpy.linspace(angle_min,
                                angle_max,
                                int(1+(angle_max-angle_min)/increment)):

        name="".join([str(count),"-",str(angle),"deg-",str(span_loc),"%"])
        
        plane=Plane(name)   #   Creates plane object with name
        plane.dihedral_angle=angle  #   Assigns dihedral angle
        plane.dihedral_split=span_loc   #   Assigns dihedral split percentage plcation
        split_loc=hspan*span_loc/100    #   Convert from %
        plane.dihedral_splitY=split_loc
        plane.span=span

        mod_geom=list(ref_plane_geom)

        Zle=round((hspan-split_loc)*math.sin(math.radians(angle)),3)    #   Calculates tip Z due to dihedral angle
        plane.tipZ=Zle
        Yle=round((hspan-split_loc)*math.cos(math.radians(angle))+split_loc,3)  #   Calcualtes tip Y due to dihedral angle
        plane.tipY=Yle
        
        root=Section(ref_plane.Xle,0,0,mac,int(split_loc*0.02),-2,aerofoil)
        split=Section(ref_plane.Xle,split_loc,0,mac,int((math.sqrt(Yle**2+Zle**2)-split_loc)*0.02),-1,aerofoil)
        tip=Section(ref_plane.Xle,Yle,Zle,mac,0,0,aerofoil)   #   Creates tip section based off tip geometry
        mod_str=root.create_input()  #   Gets section string in avl format
        mod_str+=split.create_input()
        mod_str+=tip.create_input()

        for index,line in enumerate(mod_geom):
            if line=="YES PLEASE\n":    #   Finds marker
                mod_geom.pop(index) #   Removes marker
                mod_geom.insert(index,mod_str)  #   Inserts modified sections
                
        plane.geom_file="generated planes/"+plane.name+".avl"
        with open(plane.geom_file,'w') as file:
            file.write("".join(mod_geom))
        count+=1

        planes.append(plane)
        plane.results_file=list()
        plane.cases=[case for case in analysis.initialize_cases()]  #   Assigns analysis cases

    print("Planes generated...")
    return(planes)

def run_analysis(tasks):
    time.sleep(0.001)
    plane,case,analysis=tasks
    
    analysis.analysis(plane,case)

    pass

def polars(planes):
    for plane in planes:
        polars=list()
        for case in plane.cases:
            case.Cl,case.Cd,case.Clb,case.Clp,case.spiral=Aero.read_aero(case)
            polars.append((case.alpha,case.Cl,case.Cd,case.Clb,case.Clp,case.spiral))
        plane.polars=pd.DataFrame(polars,columns=["Alpha (deg)","Cl","Cd","Clb","Clp","spiral"])

    pass

def plot_polars(planes):
    fig,(ax1,ax2,ax3)=plt.subplots(ncols=3,figsize=(12,4))
    
    dihedral_angles=[plane.dihedral_angle for plane in planes]
    Cl_0=planes[0].polars['Cl'].iloc[-1]
    Cd_0=planes[0].polars['Cd'].iloc[-1]
    Cl_delta=[100*(plane.polars['Cl'].iloc[-1]-Cl_0)/Cl_0 for plane in planes]
    Cd_delta=[100*(plane.polars['Cd'].iloc[-1]-Cd_0)/Cd_0 for plane in planes]

    ax1.plot(dihedral_angles,Cl_delta,label="Cl_delta",color='r')
    ax1.plot(dihedral_angles,Cd_delta,label="Cd_delta",color='b')
    ax1.set_xlabel("Dihedral Angles (deg)")
    ax1.set_ylabel(f"% Difference @ {planes[0].polars['Alpha (deg)'].iloc[-1]} deg")
    ax1.legend(loc='upper left')
    ax1.set_title("Aero Coeffients")

    Clb_0=planes[0].polars['Clb'].iloc[0]
    Clp_0=planes[0].polars['Clp'].iloc[0]
    Clb_delta=[100*(plane.polars['Clb'].iloc[0]-Clb_0)/Clb_0 for plane in planes]
    Clp_delta=[100*(plane.polars['Clp'].iloc[0]-Clp_0)/Clp_0 for plane in planes]
    
    ax2.plot(dihedral_angles,Clb_delta,label="Dihedral Effect Derivative",color='r')
    ax2.plot(dihedral_angles,Clp_delta,label="Roll rate derivative",color='b')
    ax2.set_xlabel("Dihedral Angles (deg)")
    ax2.set_ylabel("% Difference")
    ax2.legend()
    ax2.set_title("Stability Derivatives")

    spiral=[plane.polars['spiral'].iloc[0] for plane in planes]

    ax3.plot(dihedral_angles,spiral,color='k')   
    ax3.set_xlabel("Dihedral Angles (deg)")
    ax3.set_title("Spiral Stability (>1 = stable)")

    fig.tight_layout()

    return plt

def eigenvalues(planes,analysis): 
    for plane in planes:
        for case in plane.cases:
            if case.eigen==True:
                analysis.read_eigen(plane,case)

    dihedral_angles=[plane.dihedral_angle for plane in planes]
    roll_0=planes[0].eigen_modes["roll"][0]
    dutch_0=planes[0].eigen_modes["dutch"][0]
    roll_delta=[100*(plane.eigen_modes["roll"][0]-roll_0)/roll_0 for plane in planes]
    dutch_delta=[100*(plane.eigen_modes["dutch"][0]-dutch_0)/dutch_0 for plane in planes]

    df=pd.DataFrame(list(zip(dihedral_angles,roll_delta,dutch_delta)), columns=["Angle (deg)","Roll Damping (/s)","Dutch Damping (/s)"])
    """
    try:
        df.to_excel("dihedral results.xlsx")
    except PermissionError:
        print("! Could not write dihedral results. Close excel. !")
        pass
    """
    plt.figure(figsize=(4,4))
    plt.title("Eigenmode Damping")
    plt.xlabel(f"Dihedral Angle (deg)\nSplit Location={[plane.dihedral_split for plane in planes][0]}% of Span")
    plt.ylabel("Damping % Diff")

    plt.plot(dihedral_angles,roll_delta,color='r',label="Roll")
    plt.plot(dihedral_angles,dutch_delta,color='b',label="Dutch Roll")
    plt.legend()
    plt.tight_layout()
    
    return plt

def geom_plot(planes):
    plt.figure()
    plt.title("Spanwise Geometry Plot")
    plt.xlabel("Y (mm)")
    plt.ylabel("Z (mm)")  
    plt.xlim(0,max([plane.tipY for plane in planes]))
    plt.ylim(0,max([plane.tipY for plane in planes]))
    for plane in planes:
        plt.plot([0,plane.dihedral_splitY,plane.tipY],[0,0,plane.tipZ])
    
    return plt

if __name__=='__main__':
    os.system('cls')
    freeze_support()

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

    input_file="DIHEDRAL_CONFIG.txt"

    plt.close("all")

    run(input_file)
