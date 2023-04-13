from ase.build import bulk
from ase import Atom
from calorine import GPUNEP
from pynep.calculate import NEP
from ase.optimize import QuasiNewton
from ase.io import read, write
from ase.constraints import ExpCellFilter
from ase.optimize import FIRE
from ase.optimize import LBFGS
import numpy as np
import os

def Vac(lattice, supercell=(1,1,1), point = 'none'):
    a = lattice
    atoms = bulk('W', 'bcc', a, cubic = True)
    atoms = atoms * supercell 
    
    if point == '0':
        for atom in atoms:
            if  atom.x == lattice * 2 and atom.y == lattice * 2 and atom.z == lattice * 2:
                index = atom.index             
        del atoms[index] 
        
    if point == '1':
        for atom in atoms:
            if  atom.x == lattice * 0 and atom.y == lattice * 0 and atom.z == lattice * 0:
                index = atom.index
        del atoms[index] 
        
        for atom in atoms:
            if  atom.x == lattice * 0.5 and atom.y == lattice * 0.5 and atom.z == lattice * 0.5:
                index = atom.index       
        del atoms[index]
            
    if point == '2':
        for atom in atoms:
            if  atom.x == lattice * 2 and atom.y == lattice * 2 and atom.z == lattice * 2:
                index = atom.index
        del atoms[index] 
        
        for atom in atoms:
            if  atom.x == lattice * 3 and atom.y == lattice * 2 and atom.z == lattice * 2:
                index = atom.index       
        del atoms[index]
            
    if point == '3':
        for atom in atoms:
            if  atom.x == lattice * 2 and atom.y == lattice * 2 and atom.z == lattice * 2:
                index = atom.index
        del atoms[index] 
        
        for atom in atoms:
            if  atom.x == lattice * 3 and atom.y == lattice * 2 and atom.z == lattice * 3:
                index = atom.index       
        del atoms[index]
            
    if point == '4':
        for atom in atoms:
            if  atom.x == lattice * 2 and atom.y == lattice * 2 and atom.z == lattice * 2:
                index = atom.index
        del atoms[index] 
        
        for atom in atoms:
            if  atom.x == lattice * 3.5 and atom.y == lattice * 2.5 and atom.z == lattice * 2.5:
                index = atom.index       
        del atoms[index]
            
    if point == '5':
        for atom in atoms:
            if  atom.x == lattice * 2 and atom.y == lattice * 2 and atom.z == lattice * 2:
                index = atom.index
        del atoms[index] 
        
        for atom in atoms:
            if  atom.x == lattice * 3 and atom.y == lattice * 3 and atom.z == lattice * 3:
                index = atom.index       
        del atoms[index]
            
    return atoms
   
  
def relax(atoms, nep, f_max = 0.001 ,cell = True, qn = False, lbfgs = False, fire = False):   
    calc = NEP(nep)
    atoms.calc = calc      
    if cell:
        ucf = ExpCellFilter(atoms, scalar_pressure = 0.0, hydrostatic_strain = True)
    else:
        ucf = atoms
    if qn:
        dyn = QuasiNewton(ucf)
        dyn.run(fmax=f_max, steps = 500)    
    if lbfgs:
        dyn = LBFGS(ucf)
        dyn.run(fmax=f_max, steps = 50)
    if fire:
        dyn = FIRE(ucf)
        dyn.run(fmax=f_max, steps = 500)   
    energy = atoms.get_potential_energy()
    return energy 
        
def main():
    Lattice = 3.185
    Supercell = (4,5,6)
    a = np.array(Supercell)
    Nep = ('nep.txt')  
    E = []
    for i in ['none','1','2','3','4','5','0']:
        atoms = Vac(Lattice, Supercell, i) 
        e = relax(atoms, Nep, f_max = 0.05 , cell = True ,qn = True) 
        E.append(e)
        
    k = (-1 + 2 * (a[0]*a[1]*a[2])) / (2 * (a[0]*a[1]*a[2]))   
    E_vac   = -k * E[0] + E[6]
    k = (-2 + 2 * (a[0]*a[1]*a[2])) / (2 * (a[0]*a[1]*a[2])) 
    E_1 = 2 * E_vac - (-k * E[0] + E[1])
    E_2 = 2 * E_vac - (-k * E[0] + E[2])
    E_3 = 2 * E_vac - (-k * E[0] + E[3])
    E_4 = 2 * E_vac - (-k * E[0] + E[4])
    E_5 = 2 * E_vac - (-k * E[0] + E[5])
    
    print ('Ef_<vac>: ',   E_vac)
    print ('E_1NN: ',      E_1)
    print ('E_2NN: ',      E_2)
    print ('E_3NN: ',      E_3)
    print ('E_4NN: ',      E_4)
    print ('E_5NN: ',      E_5)
    
    f = open('vacancy.out','a')    
    print ('Ef_<vac>: ',   E_vac, file = f)
    print ('E_1NN: ',      E_1, file = f)
    print ('E_2NN: ',      E_2, file = f)
    print ('E_3NN: ',      E_3, file = f)
    print ('E_4NN: ',      E_4, file = f)
    print ('E_5NN: ',      E_5, file = f)
    f.close()
    
if __name__ == "__main__":
    main()


