from ase.build import bulk
from pynep.calculate import NEP
from calorine import GPUNEP
from ase.optimize import QuasiNewton
from ase.io import read, write
from ase.constraints import ExpCellFilter
from ase.optimize import FIRE
from ase.optimize import LBFGS
import numpy as np
import os

def Cubic(lattice, supercell):
    a = lattice
    atoms = bulk('W', 'bcc', a, cubic = True)
    atoms = atoms * (10,10,10)
    return atoms
   
def relax(atoms, nep, f_max = 0.05 ,cell = True, qn = False, lbfgs = False, fire = False, gpumd = False):            
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
  nep = ('../Train/GAP/GAP1/nep.txt')  
  atoms = bulk('W', 'bcc', 3.185, cubic = True) * (10, 10, 10)
  relax(atoms, nep, fire = True)
  
  Lattice= (atoms.get_cell()[0][0] + atoms.get_cell()[1][1] + atoms.get_cell()[2][2]) / 10 / 3
  Coh = atoms.get_potential_energy() / len(atoms)

  f = open('lc.out','a')
  print('Lattice_Constant:', Lattice, file = f)
  print('Coherent_energy:', Coh, file = f)
  f.close()

if __name__ == "__main__":
  main()
