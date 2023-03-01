from dataclasses import dataclass
from pathlib import Path

@dataclass
class Solver:
    """
    Cross Solver Configuration
    """
    threads: int
    steps: int
    path: Path


@dataclass
class Plane(Solver):
    """
    Core aircraft attributes widely used across solvers
    """
    # files
    input_plane: Path
    wing_aerofoil: Path
    elevator_aerofoil: Path
    fin_aerofoil: Path

    # aircraft data
    Cg_xyz: tuple[float, float, float] # set all to 0 if unknown
    mass: float # set to 0 if unknown


@dataclass
class TailConfig(Plane):
    """
    Tail solver specifics
    """
    # Solver inputs
    Xt_upper: float     # lever arm distance 
    Xt_lower: float     # lever arm distance
    St_h_upper: float   # Equivilent horizontal tail area
    St_h_lower: float   # Equivilent horizontail tail area
    Ct_v: float         # Equivilent vertical tail area

    sm_ideal: float     # target static margin
    tolerance: float    # tolerance to the static margin solution

    tail_config_vtail: bool

    horz_tail_span: float # set to 0 if unknown


@dataclass
class TailConfigSpanFind(TailConfig):
    """
    Used when tail span is known but the CG position is to be found
    """
    horz_tail_span: float


@dataclass
class DihedralConfig(Plane):
    """
    
    """
    angle_min: int          # min dihedral solution angle
    angle_max: int          # max dihedral solution angle
    increment: int          # solve at this increment in dihedral angle
    span_loc: int           # % distance for cranked dihedral
    show_geo_plot: bool     # show the geometry plot

