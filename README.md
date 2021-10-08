# AVL-Wrapper
Wraps avl.

### avl_aero_coefficients.py useage:

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
