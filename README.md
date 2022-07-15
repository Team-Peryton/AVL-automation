# AVL-Wrapper
CODE CLEANUP IN PROGRESS

#### Tail Sizing
- Suggests tail configurations based on input geometry to achieve a desired static margin.
- Works for convensional & A tail configurations.
- Can be run in reverse i.e. input tail dimensions & output optimal CG location

#### Dihedral
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
