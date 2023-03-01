from avlautomation.config import TailConfig
from pathlib import Path
from avlautomation.tail import AutoTail

for i in [0,0.1,0.2,0.3]:
    tail_config = TailConfig(
        threads=8,
        steps=7,
        path=Path("./"),
        input_plane=Path("Aria3.avl"),
        wing_aerofoil=Path("NACA_2412.dat"),
        elevator_aerofoil=Path("NACA_0012H.dat"),
        fin_aerofoil=Path("NACA_0012H.dat"),
        Cg_xyz=(520,0,0),
        mass=10,
        Xt_upper=1800,
        Xt_lower=800,
        St_h_upper=100000,
        St_h_lower=400000,
        Ct_v=0.06,
        sm_ideal=i,
        tolerance=0.05,
        tail_config_vtail=True,
        horz_tail_span=0
    )

    at = AutoTail(tail_config)
    at.generate_planes()
    at.run()
    at.results()