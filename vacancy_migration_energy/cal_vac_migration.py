from ase.build import bulk
from ase import Atom
from calorine import GPUNEP, CPUNEP
from pynep.calculate import NEP
from ase.optimize import QuasiNewton
from ase.io import read, write
from ase.constraints import ExpCellFilter
from ase.optimize import FIRE
from ase.optimize import LBFGS
from ase import io
from ase.neb import NEB
import numpy as np
import os

def Point(lattice, supercell=(1,1,1), point = 'none'):
    a = lattice
    atoms = bulk('W', 'bcc', a, cubic = True)
    atoms = atoms * supercell 
    
    if point == 'vacancy1':
        for atom in atoms:
            if  atom.x == lattice * 2 and atom.y == lattice * 2 and atom.z == lattice * 2:
                index = atom.index             
        del atoms[index] 
                               
    if point == 'vacancy2':
        for atom in atoms:
            if  atom.x == lattice * 2.5 and atom.y == lattice * 2.5 and atom.z == lattice * 2.5:
                index = atom.index             
        del atoms[index]       
 
    return atoms
          
def main():
    Lattice = 3.1852
    Nep = ('train/nep.txt')  
    Supercell = (4,4,4)
    initial = Point(Lattice, Supercell, 'vacancy1')
    final = Point(Lattice, Supercell, 'vacancy2')
    # Make a band consisting of 5 images:
    images = [initial]
    images += [initial.copy() for i in range(13)]
    images += [final]
    print(len(images))
    neb = NEB(images)
    # Interpolate linearly the potisions of the three middle images:
    neb.interpolate()
    # Set calculators:
    for image in images[1:15]:
        image.calc = NEP(Nep)
    # Optimize:
    optimizer = FIRE(neb)
    optimizer.run(fmax=0.01)
    energy = []
    for image in images:
        image.calc = NEP(Nep)
        energy.append(image.get_potential_energy())
    f = open('vac_migration.out','a')
    print(energy, file = f)
    print(energy[7]-energy[0], file = f)
    f.close()
   
if __name__ == "__main__":
    main()


