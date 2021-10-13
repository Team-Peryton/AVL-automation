# AVL-Wrapper
Wraps avl. Code is very messy. You have been warned.

#### Tail Sizing
- Suggests tail configurations based on input geometry to achieve a desired static margin.
- Works for convensional & inverted V configurations.
- Can be run in reverse i.e. input tail dimensions & output optimal CG location

#### Dihedral
- Generates wing configurations with varying dihedral angle and set spanwise dihedral location. 
- Calculates roll and dutch roll damping modes for each configuration.
- Runs aero analysis, generating Cl and Cd polars.


#### avl_aero_coefficients.py useage:

```
analysis=Aero("AERO_CONFIG.txt")
cases=analysis.initialize_cases()
#generate plane objects (geometry.py)
for plane in planes:
  plane.cases=cases
  for case in plane.cases:
    analysis.analysis(plane,plane.case)
analysis.read_aero(case)
```

#### Required files in same directory:
- avl.exe
- Config files
- Input plane .avl file
- aerofoil .dat files
