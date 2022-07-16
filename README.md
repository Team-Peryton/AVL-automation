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

The main 3 variables to consider when sizing the horizontal tail for longitudinal static stability are: tail moment arm, tail plane area, and CG position. In short, increasing the moment arm and area increase the longitudinal stability of the aircraft for a given CG because it moves the neutral point away from the CG. The neutral point is the point where, if the CG was placed on it, $C_M/\alpha=0$. 


<p align="center">
  ![image](https://user-images.githubusercontent.com/79290428/179371409-41668bca-5fb2-44d5-ba4a-9e70294bdff1.png)
</p>

The tail moment arm is likely restricted in some way because of structural constraints so there is a balance between the tail area and moment arm.

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
