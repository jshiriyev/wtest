from dataclasses import dataclass

@dataclass
class rrock:
    """
    perm    : permeability in mD
    poro    : porosity
    height  : formation height, ft
    comp    : rock compressibility, 1/psi
    """
    perm    : float = None
    poro    : float = None
    comp    : float = None
    height  : float = None
    radius  : float = None