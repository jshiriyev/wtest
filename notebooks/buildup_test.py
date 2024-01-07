import numpy as np

from borepy.wtest import reservoir,well,fluid
from borepy.wtest import buildup

res = reservoir(permeability=200,porosity=0.15,thickness=30)

oil = fluid(viscosity=1.5,formation_volume_factor=1.2,density=60)

prod = well(rate=[800,0],time=[16,20],radius=0.25,skin=0)

Pbu = buildup(prod)

Pbu.set_parameters(res,oil)

Pbu.set_compressibility(total=25e-6)

Pbu.initialize(scale="log",size=1000)

Pbu.view()
# Pbu.horner()

