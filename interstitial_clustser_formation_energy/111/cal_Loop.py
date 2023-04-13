from ase.build import bulk
from ase import Atom
from calorine import GPUNEP
import numpy as np
import os

def Loop(lattice, supercell=(30,30,30), defect = '100', Rcut = 2, b = 3):
    a = lattice
    atoms = bulk('W', 'bcc', a, cubic = True)
    atoms = atoms * supercell 
    x = (supercell[0] / 2) * a
    y = (supercell[0] / 2) * a
    z = (supercell[0] / 2) * a
    index = []
    if defect == '100':
        for atom in atoms:
            R = ((atom.x - x)**2 +(atom.y - y)**2 + (atom.z - z)**2)**(1/2)
            if  R < Rcut and atom.z < z + 1/2 * b and atom.z > z - 1/2 * b:
                index.append(atom.index)  
                atom.z = atom.z + 1/3 * b
        for i in index:            
            interz = Atom('W',(atoms[i].x,atoms[i].y, atoms[i].z -2/3 * b))
            atoms.append(interz)              
  
    if defect == '111':
        for atom in atoms:
            R = ((atom.x - x)**2 +(atom.y - y)**2 + (atom.z - z)**2)**(1/2)
            if  R < Rcut and atom.x < atom.y + 1.4142 / 2 * b and  atom.x > atom.y - 1.4142 / 2 * b:
                index.append(atom.index)  
                atom.x = atom.x + 1.4142 /2 * 1/3 * b 
                atom.y = atom.y - 1.4142 /2 * 1/3 * b 
        for i in index:            
            interz = Atom('W',(atoms[i].x - 1.4142 /2 * 2/3 * b, atoms[i].y + 1.4142 /2 * 2/3 * b, atoms[i].z))
            atoms.append(interz)              
    return atoms
        
def relax(atoms, nep, md_dir='./'):                      
    os.chdir(md_dir)
    calc = GPUNEP('../'+nep, directory = './', maximum_neighbors = 200)
    atoms.calc = calc 
    parameters = [('ensemble', 'nve'), ('minimize', ('sd', '1.0e-4','1000')) ,('time_step','0'), ('dump_thermo',1),('dump_position',1),('run',1)]
    calc.run_custom_md(parameters)
    f = open('./thermo.out')
    data = f.readlines()
    energy = float(data[-1].split()[2])
    f.close()
    os.chdir('../')
    return energy

def main():
    Lattice = 3.1852
    Nep = ('nep.txt')  
    Supercell = (50,50,50)
    E = []
    for r in range(20):
        atoms = Loop(Lattice, Supercell, defect = '111', Rcut = r, b = 3) 
        interz = len(atoms) - Supercell[0] * Supercell[1] * Supercell[2] * 2
        path = './'+str(interz)
        isExists = os.path.exists(path)
        if not isExists:						
            os.makedirs(path)	
            print("%s success"%interz)
        else:
            print("%s dir exist"%interz)	
        e = relax(atoms, Nep, path) 
        E.append([interz, e])   
        
    f = open('111.out', 'a')
    print(len(E))
    for i in range(len(E)):
        print(E[i][0],E[i][1],file = f)
    f.close()             
   
if __name__ == "__main__":
    main()


