from ase.build import bulk
from ase.build import surface
from ase.spacegroup import crystal
from pynep.calculate import NEP
from calorine import GPUNEP
from ase.optimize import QuasiNewton
from ase.constraints import ExpCellFilter
from ase.optimize import FIRE
from ase.optimize import LBFGS
import numpy as np
import os

def Surf(lattice, nep, miller):
    a = lattice
    m = miller
    atoms = crystal('W',[0,0,0] ,spacegroup=229, cellpar=[a,a,a,90,90,90])
    slab = surface(atoms, m, 10, vacuum=10)
    slab = slab * (10,10,1)
    cubic = bulk('W', 'bcc', a, cubic = True) * (10, 10, 10)
    
    relax(slab, nep, qn = True)  
    relax(cubic, nep, qn = True)  
    
    box = slab.get_cell()
    S = (box[0][0] * box[1][1] - box[0][1] * box[1][0]) 
    slab_energy = slab.get_potential_energy()
    cubic_energy = cubic.get_potential_energy()
    surf_energy = - (cubic_energy-slab_energy) / (2 * S)
    f = open('surf.out', 'a')
    print(miller, 'Surface_Energy:', surf_energy * 1000, 'meV/A^2', file = f)
    f.close()

def relax(atoms, nep, f_max = 0.005 ,cell = True, qn = False, lbfgs = False, fire = False):             
    calc = NEP(nep)
    atoms.calc = calc     
    if cell:
        ucf = ExpCellFilter(atoms, scalar_pressure = 0.0, hydrostatic_strain = True)
    else:
        ucf = atoms
    if qn:
        dyn = QuasiNewton(ucf)
        dyn.run(fmax=f_max)    
    if lbfgs:
        dyn = LBFGS(ucf)
        dyn.run(fmax=f_max)
    if fire:
        dyn = FIRE(ucf)
        dyn.run(fmax=f_max)      

def main():
    Lattice = 3.1852
    Nep = ('train/nep.txt')
    Miller=[(1,1,0),(1,0,0),(1,1,1),(2,1,1),(2,1,0),(2,2,1),(3,1,1),(3,1,0),(3,2,1),(3,2,0)]
    for i in Miller:
        Surf(Lattice, Nep, i)

if __name__ == "__main__":
    main()


