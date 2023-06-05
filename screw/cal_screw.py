from ase.neb import NEB
from ase import Atoms
from pynep.calculate import NEP
from ase.optimize import FIRE
import matplotlib.pyplot as plt
import numpy as np

def dump_xyz(filename, atoms, comment = ''):
    with open(filename, 'a') as f:
        Out_string = ""
        Out_string += str(int(len(atoms))) + "\n"
        Out_string += "pbc=\"" + " ".join(["T" if pbc_value else "F" for pbc_value in atoms.get_pbc()]) + "\" "
        Out_string += "Lattice=\"" + " ".join(list(map(str, atoms.get_cell().reshape(-1)))) + "\" "
        Out_string += "Properties=species:S:1:pos:R:3:mass:R:1"
        Out_string += comment + "\n"
        s = atoms.get_chemical_symbols()
        p = atoms.get_positions()
        m = atoms.get_masses()
        for j in range(int(len(atoms))):              
            Out_string += '{:2} {:>15.8e} {:>15.8e} {:>15.8e} {:>15.8e}\n'.format(s[j], *p[j], m[j])
        f.write(Out_string)

def read_xyz(filename):
    with open(filename) as f:
        lines = f.readlines()
        frames = []
        while lines:
            symbols = []
            positions = []
            natoms = int(lines.pop(0))
            comment = lines.pop(0)  # Comment line; ignored
            if "pbc=\"" in comment:
                pbc_str = comment.split("pbc=\"")[1].split("\"")[0].strip()
                pbc = [True if pbc_value == "T" else False for pbc_value in pbc_str.split()]
            else:
                pbc = [True, True, True]
            lattice_str = comment.split("Lattice=\"")[1].split("\" Properties=")[0].strip()
            lattice = [list(map(float, row.split())) for row in lattice_str.split(" ")]
            cell = [lattice[0] + lattice[1] + lattice[2], lattice[3] + lattice[4] + lattice[5], lattice[6] + lattice[7] + lattice[8]]
            for _ in range(natoms):
                line = lines.pop(0)
                symbol, x, y, z = line.split()[:4]
                symbol = symbol.lower().capitalize()
                symbols.append(symbol)
                positions.append([float(x), float(y), float(z)])
            frames.append(Atoms(symbols=symbols, positions=positions, cell = cell, pbc = pbc))
    return frames

MV4S = [0,0,0,0,0.00272027,0.018451292,0.042325469,0.055881881,0.059555449,0.055453378,0.035437329,0.008974482,0,0,0,0]
MV2B = [0,0,0.001211684,0.014145402,0.033933558,0.047929706,0.051603274,0.045556091,0.027274916,0.007043813,0,0,0]
DFT_csanyi = [0,0.048115872,0.10106564,0.071874498,0.03106564]
DFT_Rodney = [0,0.007497994,0.030386776,0.057786872,0.079858771,0.087961804,0.078850907,0.056397047,0.028762638,0.00496068,0]
DFT_stoller = [0,0.005376344,0.023318889,0.046451613,0.068284385,0.078709677,0.067839833,0.045806452,0.021808698,0.003020382,0]

initial = read_xyz('initial_135at_periodic.xyz')[0]
final = read_xyz('final_135at_periodic_dipole_move.xyz')[0]
images = [initial] + [initial.copy() for i in range(9)] + [final]
for i in images:
    i.calc = NEP('../train/nep.txt')
neb = NEB(images)
neb.interpolate()    
optimizer = FIRE(neb)
optimizer.run(fmax=0.02)
energies = [image.get_potential_energy() / 2  for image in images]
energies -= min(energies) 
for image in images:
    dump_xyz('screw.xyz', image, comment=f' config_type = screw')  

plt.plot(np.linspace(0, 1, len(energies)), energies, marker='o', label='UNEP')
plt.plot(np.linspace(0, 1, len(MV4S)), MV4S, marker='o', label='MV4S')
#plt.plot(np.linspace(0, 1, len(MV2B)), MV2B, marker='o')
plt.plot(np.linspace(0, 1, len(DFT_csanyi)), DFT_csanyi, marker='o', label='DFT-Csanyi')
plt.plot(np.linspace(0, 1, len(DFT_Rodney)), DFT_Rodney, marker='o', label='DFT-Rodney')
plt.plot(np.linspace(0, 1, len(DFT_stoller)), DFT_stoller, marker='o', label='DFT-Stoller')
plt.legend()
plt.savefig('screw_dipole_move.png')
