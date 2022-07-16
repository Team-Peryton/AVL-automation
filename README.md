# AVL-automation
Surrey Team Peryton Heron IMechE UAS 2022 AVL automation program. Enabled automatic tail sizing based on input aircraft geometry, and investigation of dihedral angle effects on stability and aerodynamic performance. 

A lose wrapping of the vortex lattice method AVL provides a means of scripting these tools.

## Usage
1. Install python (tested on 3.9.5)
2. Install packages in requirements.txt
   >pip install -r requirements.txt
3. Run avl-automation.py from the command line. For help:
   >py avl-automation.py -h

## Tail Sizing
- Suggests tail configurations based on input geometry to achieve a desired static margin.
- Works for convensional & V-tail configurations.
- Can be run in reverse i.e. input tail dimensions & output optimal CG location
- Anslysis considers only the horizontal tail plane (if convensional tail is selected). Vertical tail plane should be sized independently based on rudder and yaw stability requirements.



## Dihedral
- Generates wing configurations with varying dihedral angle with a set spanwise dihedral location. 
- Calculates roll and dutch roll damping modes for each configuration.
- Runs aero analysis, generating Cl and Cd polars, and relevant stability derivatives.

- AVL solution becomes very unstable as dihedral span location approaches tail sections in y.

#### aero.py useage:
example useage to come

#### Required files in same directory:
- avl.exe
- Config files
- Input plane .avl file
- aerofoil .dat files
