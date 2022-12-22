import argparse
import os

from .aero import Aero
from .geometry import Plane
from .dihedral import Dihedral
from .tail import AutoTail

# if os.path.exists('avl.exe')==False:
#     print("\u001b[31m[Error]\u001b[0m avl.exe not found.")
#     exit()

parser=argparse.ArgumentParser(description="AVL Automation.")

parser.add_argument('run_type',choices=['aero','tail','dihedral'],help='Type of analysis to run.')
parser.add_argument('-p','--plane',action='store',help="Plane .avl file for aero analysis.")
parser.add_argument('-c','--config',nargs='+',action='store',help="Config file for analysis.")

args=parser.parse_args()

if args.run_type=='aero':
    if args.plane is None:
        parser.error("Aero requires --plane.")

    if args.config is None:
        parser.error("Aero requires --config.")

    if len(args.config)>1:
        parser.error("Aero requires only 1 config file.")

    if os.path.exists(args.config[0])==False:
        print(f"\u001b[31m[Error]\u001b[0m {args.config[0]} not found.")
        exit()

    plane=Plane(geom_file=args.plane)

    aero=Aero(args.config[0])
    aero.run(plane)

    if aero.polars==True:
        print('\nPolars:\n',plane.polars)
    if aero.modes==True:
        print('\nEigenmodes:\n',plane.modes,'\n')

if args.run_type=='tail':
    if args.config is None:
        parser.error("Aero requires --config.")
        
    if len(args.config)>1:
        parser.error("Aero requires only 1 config file.")

    if os.path.exists(args.config[0])==False:
        print(f"\u001b[31m[Error]\u001b[0m {args.config[0]} not found.")
        exit()
    
    tail=AutoTail(args.config[0])
    tail.generate_planes()
    tail.run()
    tail.results()

if args.run_type=='dihedral':
    if args.config is None:
        parser.error("Aero requires --config.")
    
    if len(args.config)!=2:
        parser.error("Aero requires 2 config files: dihedral, aero.")

    for config in args.config:
        if os.path.exists(config)==False:
            print(f"\u001b[31m[Error]\u001b[0m {config} not found.")
            exit()

    
    dihedral=Dihedral(args.config[0],args.config[1])
    dihedral.generate_planes()
    dihedral.run()
    dihedral.plot()