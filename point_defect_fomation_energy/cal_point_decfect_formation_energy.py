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

def Point(lattice, supercell=(1,1,1), point = 'none'):
    a = lattice
    atoms = bulk('W', 'bcc', a, cubic = True)
    atoms = atoms * supercell 
    
    if point == 'vac':
        for atom in atoms:
            if  atom.x == lattice * 2 and atom.y == lattice * 2 and atom.z == lattice * 2:
                index = atom.index             
        del atoms[index] 
                               
    if point == '111':
        for atom in atoms:
            if  atom.x == lattice * 2 and atom.y == lattice * 2 and atom.z == lattice * 2:
                index = atom.index             
        del atoms[index]
        i1 = Atom('W',(2.25 * a, 2.25 * a ,2.25 * a))
        i2 = Atom('W',(1.75 * a, 1.75 * a ,1.75 * a))
        atoms.append(i1)
        atoms.append(i2)        
        
    if point == '110':
        for atom in atoms:
            if  atom.x == lattice * 2 and atom.y == lattice * 2 and atom.z == lattice * 2:
                index = atom.index
        del atoms[index]
        i1 = Atom('W',(2.25 * a, 2.25 * a , 2 * a))
        i2 = Atom('W',(1.75 * a, 1.75 * a , 2 * a))
        atoms.append(i1)
        atoms.append(i2)     
        
    if point == '100':
        for atom in atoms:
            if  atom.x == lattice * 2 and atom.y == lattice * 2 and atom.z == lattice * 2:
                index = atom.index             
        del atoms[index]
        i1 = Atom('W',(2.25 * a, 2 * a , 2 * a))
        i2 = Atom('W',(1.75 * a, 2 * a , 2 * a))
        atoms.append(i1)
        atoms.append(i2)     
        
    if point == 'octa':
        oc = Atom('W',(2 * a, 2 * a ,1.5 * a))
        atoms.append(oc)
        
    if point == 'tetra':
        te = Atom('W',(1.75 * a, 1.5 * a, 2 * a))
        atoms.append(te)
        
    return atoms
        
def relax(atoms, nep, f_max = 0.001 ,cell = True, qn = False, lbfgs = False, fire = False, gpumd = False, md_dir='./'):   
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
    Lattice = 3.1852
    Supercell = (4,5,6)
    a = np.array(Supercell)
    k = (1 + 2 * (a[0]*a[1]*a[2])) / (2 * (a[0]*a[1]*a[2]))
    Nep = ('nep.txt')  
    E = []
    for i in ['none','111','110','100','octa','tetra','vac']:
        atoms = Point(Lattice, Supercell, i) 
        e = relax(atoms, Nep, f_max = 0.0001 , qn = True, md_dir = i) 
        E.append(e)
        
    E_111   = -k * E[0] + E[1]
    E_110   = -k * E[0] + E[2]
    E_100   = -k * E[0] + E[3]
    E_octa  = -k * E[0] + E[4]
    E_tetra = -k * E[0] + E[5]
    k = (-1 + 2 * (a[0]*a[1]*a[2])) / (2 * (a[0]*a[1]*a[2]))   
    E_vac   = -k * E[0] + E[6]
    f = open('defect.out','a')    
    print ('Ef_<111>: ',   E_111, file = f)
    print ('Ef_<110>: ',   E_110, file = f)
    print ('Ef_<100>: ',   E_100, file = f)
    print ('Ef_<octa>: ',  E_octa, file = f)
    print ('Ef_<tetra>: ', E_tetra, file =f)
    print ('Ef_<vac>: ',   E_vac, file =f)
    f.close()
    
if __name__ == "__main__":
    main()


