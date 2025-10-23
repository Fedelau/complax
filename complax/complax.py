import sys
import numpy as np
import math
import os
import time
import itertools
import subprocess
from colorama import Fore, Style
import shutil
import argparse, textwrap
from ase.io import read, write
from ase import Atoms
import threading
from functools import partial
from multiprocessing import Pool
import tqdm
from tabulate import tabulate
import pkg_resources


def print_banner():
    os.system("clear")

    # --- layout ---
    left_margin = 11  
    titolo = [
        " ██████╗ ██████╗ ███╗   ███╗██████╗ ██╗      █████╗ ██╗  ██╗",
        "██╔════╝██╔═══██╗████╗ ████║██╔══██╗██║     ██╔══██╗╚██╗██╔╝",
        "██║     ██║   ██║██╔████╔██║██████╔╝██║     ███████║ ╚███╔╝ ",
        "██║     ██║   ██║██║╚██╔╝██║██╔═══╝ ██║     ██╔══██║ ██╔██╗ ",
        "╚██████╗╚██████╔╝██║ ╚═╝ ██║██║     ███████╗██║  ██║██╔╝ ██╗",
        " ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚═╝     ╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝"
    ]

    raggio_esterno = 4 
    raggio_interno = 1.5
    fetta = []
    for y in range(-raggio_esterno, raggio_esterno + 1):
        line = []
        for x in range(-raggio_esterno * 2, raggio_esterno * 2 + 1):
            dist = math.sqrt((x / 2) ** 2 + y ** 2)
            if raggio_interno < dist < raggio_esterno:
                ang = math.degrees(math.atan2(y, x / 2)) % 360
                if any(abs(ang - d) < 9 for d in range(0, 360, 45)):
                    line.append("*")   
                else:
                    line.append("@")  
            else:
                line.append(" ")
        fetta.append("".join(line))

    titolo_width = max(len(r) for r in titolo)
    fetta_width = len(fetta[0])
    pad_left_for_fetta = left_margin + (titolo_width - fetta_width) // 2

    print("\n" * 1)
    for r in fetta:
        print(" " * pad_left_for_fetta + r)
    print()

    for r in titolo:
        print(" " * left_margin + r)

    print()
    print(" " * (left_margin) + " ____________________________________________________________ ")
    print(" " * (left_margin) + "|                                                            |")
    print(" " * (left_margin) + "|        COMPLAX — Solvation & Optimization Tool             |")
    print(" " * (left_margin) + "|  Developed by Federica Lauria, University of Turin (2025)  |")
    print(" " * (left_margin) + "|____________________________________________________________|")
    print("\n")
     

def openfile(file):
    with open(file, 'r') as xyz_file:
        lines = xyz_file.readlines()[2:]  # salto le prime due righe

    # mantieni solo le righe con almeno 4 elementi
    lines = [line for line in lines if len(line.split()) >= 4]

    atomic_symbols = [line.split()[0] for line in lines]
    atomic_coordinates = np.array([line.split()[1:4] for line in lines], dtype=float)

    return atomic_symbols, atomic_coordinates
    
def normalize(v):
    return v / np.linalg.norm(v)

def check_overlap(mol_list, mol, cutoff=1.2):
    """Controlla che mol non si sovrapponga con nessuna molecola in mol_list"""
    pos2 = mol.get_positions()
    for mol1 in mol_list:
        pos1 = mol1.get_positions()
        for a in pos1:
            for b in pos2:
                if np.linalg.norm(a - b) < cutoff:
                    return True
    return False

def random_rotation_matrix(axis):
    """Matrice di rotazione casuale attorno a un asse"""
    axis = normalize(axis)
    theta = np.random.rand() * 2 * np.pi
    c, s = np.cos(theta), np.sin(theta)
    x, y, z = axis
    R = np.array([
        [c+(1-c)*x*x,     (1-c)*x*y - s*z, (1-c)*x*z + s*y],
        [(1-c)*y*x + s*z, c+(1-c)*y*y,     (1-c)*y*z - s*x],
        [(1-c)*z*x - s*y, (1-c)*z*y + s*x, c+(1-c)*z*z    ]
    ])
    return R

spinner_used = False

def spinner_func(molecola, stop_event):
    for c in itertools.cycle(['|', '/', '-', '\\']):
        if stop_event.is_set():
            break
        sys.stderr.write(f'\rProcessing {molecola}... {c}')
        sys.stderr.flush()
        time.sleep(0.1)
    sys.stderr.write('\r') 
    sys.stderr.flush()

def task(molecola, alpb, gbsa, con, lev, chrg, uhf, proc):
    global spinner_used

    use_spinner = False
    if not spinner_used:
        spinner_used = True
        use_spinner = True
        
    cmd_parts = [
        f"xtb {molecola}.xyz",    
        f"--opt {lev}",            
        f"--input xtb{con}.inp",   
        f"--namespace {molecola}", 
        f"--chrg {chrg}",
        f"--uhf {uhf}",          
        f"-P {proc}"               
    ]
    
    if alpb:
        cmd_parts.insert(3, f"--alpb {alpb}")
    elif gbsa:
        cmd_parts.insert(3, f"--gbsa {gbsa}")

    cmd = " ".join(cmd_parts)
        
    #cmd = f"xtb {molecola}.xyz --opt {lev} --input xtb{con}.inp {alpb} --namespace {molecola} --chrg {chrg} -P {proc}"

    if use_spinner:
        stop_event = threading.Event()
        t = threading.Thread(target=spinner_func, args=(molecola, stop_event))
        t.start()

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if use_spinner:
        stop_event.set()
        t.join()

    with open(f"{molecola}.out", "w") as f:
        f.write(result.stdout)

    if "abnormal" in result.stderr.lower():
        print(Fore.RED + result.stderr + Style.RESET_ALL, molecola)

        
class SpacedHelpFormatter(argparse.HelpFormatter):
    def add_argument(self, action):
        super().add_argument(action)
        self._add_item(lambda *args: "", [])
        
class BannerArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        print_banner() 
        super().print_help(file)

def main():
 
    parser = BannerArgumentParser(
        usage='%(prog)s <file1.xyz> <file2.xyz> [options]',
        description=textwrap.dedent('''
            COMPLAX: A tool that places solvent molecules around a solute molecule at a specified distance,
            ensuring no overlaps, and performs a geometry optimization.
        '''),

        epilog="For further information, please contact the programme author.",
        formatter_class=lambda prog: SpacedHelpFormatter(prog, max_help_position=40, width=95) #SmartFormatter(prog, max_help_position=35, width=90)
    )
    
    try:
        version = pkg_resources.get_distribution('complax').version
    except pkg_resources.DistributionNotFound:
        version = "unknown" # Fallback se non è installato come pacchetto
    
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'%(prog)s {version}',
        help='Show the version number and exit.'
    )
    
    parser.add_argument(
        'file1',
        type=str,
        nargs='?',
        help=textwrap.dedent('''
                             The solute molecule file (e.g., molecule.xyz). 
                             This file represents the molecule coordinates that 
                             will coordinate with the solvent molecule."
                              ''')
    )
    
    parser.add_argument(
        'file2',
        type=str,
        nargs='?',
        help=textwrap.dedent('''
                             The solvent molecule file (e.g., solvent.xyz). 
                             This file contains the solvent coordinates with 
                             which the solute will interact." 
                             ''')
    )

    parser.add_argument(
        '--alpb',
        type=str,
        metavar='SOLVENT',
        choices=[
        'acetone', 'acetonitrile', 'aniline', 'benzaldehyde', 'benzene',
        'ch2cl2', 'chcl3', 'cs2', 'dioxane', 'dmf', 'dmso', 'ether',
        'ethylacetate', 'furane', 'hexandecane', 'hexane', 'methanol',
        'nitromethane', 'octanol', 'woctanol', 'phenol', 'toluene', 'thf', 'water'
    ],
        help=textwrap.dedent('''Analytical linearized Poisson-Boltzmann (ALPB) model,
available solvents are acetone, acetonitrile, aniline, benzaldehyde,
benzene, ch2cl2, chcl3, cs2, dioxane, dmf, dmso, ether, ethylacetate, furane,
hexandecane, hexane, methanol, nitromethane, octanol, woctanol, phenol, toluene,
thf, water. 
                             ''')                       
    )
    
    parser.add_argument(
        '--gbsa',
        type=str,
        metavar='SOLVENT',
        choices=['acetone', 'acetonitrile', 'benzene', 'CH2Cl2', 'CHCl3', 'CS2', 'DMF', 'DMSO', 'ether', 'H2O', 'methanol', 'n-hexan', 'THF', 'toluene'
    ],
        help=textwrap.dedent('''Generalized Born model with a simple switching function (GBSA), 
                             available solvents are acetone, acetonitrile, benzene (only GFN1-xTB), 
                             CH2Cl2, CHCl3, CS2, DMF (only GFN2-xTB), DMSO, ether, H2O, methanol, 
                             n-hexane (only GFN2-xTB), THF and toluene.
                             ''')
    )

    parser.add_argument(
        '-a',
        type=int,
        nargs=2,
        metavar=('MOLECULE_ATOM', 'SOLVENT_ATOM'),
        help=textwrap.dedent(''' 
                             Atom numbers of molecule and solvent.
                             Format: MOLECULE_ATOM SOLVENT_ATOM
                             The number has to be specified using 1-based indexing. 
                             ''')
    )
    
    parser.add_argument(
        '-c',
        type=int,
        metavar="NO. OF COPIES",
        default=1,
        help="Number of file2 copies to be placed around the selected atom of file1. Default is 1."        
    )
    
    parser.add_argument(
        '-t',
        type=float,
        metavar="TARGET DISTANCE",
        default=2.0,
        help=textwrap.dedent(''' 
                             Target distace from 'MOLECULE ATOM', in Ångstrom. Default is 2.0 Ångstrom
                             ''') 
    ) 

    parser.add_argument(
        '-p',
        type=int,
        metavar="INT",
        default=1,
        help=textwrap.dedent(''' 
                             Numbers of parallel processes. Default=1
                             Number of parallel processes. Default=1.
                             During the initial optimization, the program will use the specified number 
                             of parallel processes. For subsequent optimizations (one for each solvent 
                             configuration), it will automatically launch as many parallel calculations 
                             as the number of solvent molecules selected.
                                ''')
    )
    
    parser.add_argument(
        '--lev',
        type=str,
        metavar="LEVEL",
        default="--gfn2",
        help=textwrap.dedent(''' 
                             Level of theory for the optimization. Default is --gnf2 (geometry optimization).
                             Other options include --gfn0, --gfn1, --gfn2, --gfnff.
                             ''')
    )
    
    parser.add_argument(
        '--chrg',
        type=int,
        metavar="INT",
        default=0,
        help=textwrap.dedent('''Molecular charge. Default is 0.
                             ''')
    )
    
    parser.add_argument(
        '-u','--uhf',
        type=int,
        metavar="INT",
        default=1,
        help=textwrap.dedent('''
                                Number of unpaired electrons. Default is 1.
                                ''')
    )
    
    parser.add_argument(
        '--maxtries',
        type=int,
        metavar="INT",
        default=1000,
        help=textwrap.dedent(''' 
                             Maximum number of attempts to place each solvent molecule without overlaps.
If the desired number of solvent molecules cannot be positioned, try increasing this value.
However, if placement remains difficult, it is likely due to steric hindrance between the molecules. Default is 1000.
                             ''')
    )
    
    parser.add_argument(
        '--solvfx',
        action='store_true',
        help=textwrap.dedent('''
                             If specified, do a evaluation of the effect of the solvation in term of potential energy among the different systems with an increasing number of solvent molecules.
                            ''')
    )
    
    args = parser.parse_args()

    if not args.file1 or not args.file2:
        parser.print_usage()  
        print("Error: You must specify <file1> and <file2>. Use -h or --help for more information.")
        sys.exit(1)
    
    molecule_atom, solvent_atom = args.a if args.a else (None, None)
    
    # INIZIO COMPLAX

    print_banner()
    
    atomic_coordinates = openfile(args.file1)
    
    atomic_coordinates_s = openfile(args.file2)

    last_index_reag = len(atomic_coordinates)
    index_I_solv = last_index_reag+1
    last_index_I_solv = len(atomic_coordinates_s)+len(atomic_coordinates)

    atom_1 = 1
    atom_2 = len(atomic_coordinates)
    atom_3 = index_I_solv
    atom_4 = last_index_I_solv
 
    with open('xtb.inp', 'w') as f:
        f.write('$fix')
        f.write('\n')
        f.write(f'atoms: {atom_1}-{atom_2}')
        f.write('\n')
        f.write('$end')
        f.write('\n')
        f.write('$constrain')
        f.write('\n')
        f.write('force constant=5.0')
        f.write('\n')
        f.write(f'atoms: {atom_3}-{atom_4}')
        f.write('\n')
        f.write('$end')

    dir = 'outplax'
    if os.path.exists(dir):
        shutil.rmtree(dir)
    os.makedirs("outplax")

    os.system(f"cp {args.file1} outplax")
    os.system(f"cp {args.file2} outplax")
    os.system("cp xtb.inp outplax")
    os.chdir("outplax")

    level = args.lev
    constrain = ''

    molA = read(args.file1)
    molB = read(args.file2)

    iA = int(molecule_atom-1)
    iB = int(solvent_atom-1)
    
    print(f"{molA[iA].symbol}{molecule_atom} was selected for {args.file1}")
    print(f"{molB[iB].symbol}{solvent_atom} was selected for {args.file2}")
    print("")
    
    dist = int(args.t)
    
    n_copies = args.c
    
    try:
        posA = molA[iA].position
    except IndexError:
        print(Fore.RED + f"❌ Error: the atom index {iA+1} is out of range for {args.file1}!" + Style.RESET_ALL)
        print(Fore.YELLOW + f"Hint: check the number of atoms in {args.file1} (the atom indices must be specified using 1-based indexing. They are numbered from 1 to {len(molA)})." + Style.RESET_ALL)
        sys.exit(1)
    
    try:
        posB = molB[iB].position
    except IndexError:
        print(Fore.RED + f"❌ Error: the atom index {iB+1} is out of range for {args.file2}!" + Style.RESET_ALL)
        print(Fore.YELLOW + f"Hint: check the number of atoms in {args.file2} (the atom indices must be specified using 1-based indexing. They are numbered from 1 to {len(molA)})." + Style.RESET_ALL)
        sys.exit(1)
    
    placed_molecules = [molA]
    
    max_tries = args.maxtries
    
    for copy in range(n_copies):
        success = False
        for attempt in range(max_tries):
            trialB = molB.copy()

            trialB.translate(-posB)

            direction = np.random.randn(3)
            direction = normalize(direction)

            new_posB = posA - dist * direction

            trialB.translate(new_posB)

            axis = normalize(posA - new_posB)
            R = random_rotation_matrix(axis)
            trialB.positions = (trialB.positions - new_posB) @ R.T + new_posB

            if not check_overlap(placed_molecules, trialB, cutoff=1.2):
                placed_molecules.append(trialB)
                print(f"🧭 {args.file2} no.{copy+1} successfully placed after {attempt+1} attempts")
                success = True
                break

        if not success:
            print(f"❌ Unable to place {args.file2} instance {copy+1}")
       
    mol_list = []
    out_list = []
    
    
    for n in range(1, len(placed_molecules)):  
        combined = Atoms()
        for mol in placed_molecules[:n+1]: 
            combined += mol
        write(f"complax_input_{n}solvent.xyz", combined)
        print(f"💾 Saved system with {n} solvent molecule(s) to complax_input_{n}solvent.xyz")
    
    molAnoext = os.path.splitext(args.file1)[0]
    out_list.append(molAnoext)
    molBnoext = os.path.splitext(args.file2)[0]
    out_list.append(molBnoext)
         
    for n in range(1, n_copies+1):
        filename = f"complax_input_{n}solvent.xyz"
        if os.path.exists(filename):
            mol_list.append(filename)   
            out_name = os.path.splitext(filename)[0]  
            out_list.append(out_name)
                
    combined = Atoms()
    for mol in placed_molecules:
        combined += mol
    write(f"complax_input.xyz", combined)
    
    
    print(f"💾 Final system has been saved to complax_input.xyz")
    
    print("")
    print("-------------------------------------------------------------------")
    print(f"---------- Geometry optimization with {Fore.YELLOW}{args.c} solvent{Fore.RESET} molecule ----------")
    print("-------------------------------------------------------------------")
    print("")

    task(molecola="complax_input", alpb=args.alpb, gbsa=args.gbsa, con=constrain, lev=level, chrg=args.chrg, uhf=args.uhf, proc=int(args.p) )#, out=f"opt{noext}")
    
    print(f"✅ Optimized geometry has been saved to {Fore.GREEN}complax_input.xtbopt.xyz{Fore.RESET}")
    
    task_with_args = partial(task,
                            alpb=args.alpb,
                            gbsa=args.gbsa,
                            con=constrain,
                            lev=args.lev,
                            chrg=args.chrg,
                            uhf=args.uhf,
                            proc=1)

    with Pool(args.p) as pool:
        for _ in tqdm.tqdm(
            pool.imap_unordered(task_with_args, out_list),
            total=len(mol_list),
            desc=f"{Fore.YELLOW}Calculating ...{Style.RESET_ALL}",
            bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [Elapsed Time:{elapsed} ETA:{remaining}]'
        ):
            pass
        
    i, j = args.a
    
    distances = []
    alldistances = []
    
    for n in range(1, n_copies + 1):
        outfile = f"complax_input_{n}solvent.xtbopt.xyz"
        if os.path.exists(outfile):
            outmol = read(outfile)
            indexS = i - 1 
            N_solute = len(molA)
            N_solvent = len(molB)
            
            for copy_num in range(1, n + 1):
                solvent_index = N_solute + (copy_num - 1) * N_solvent + (j - 1)
                dist = outmol.get_distance(indexS, solvent_index)
                if dist > args.t + 1.0:
                    print(Fore.RED + f"⚠️ Warning: In {outfile}, solvent molecule no.{copy_num} is too far from the selected atom (distance = {dist:.2f} Å)." + Style.RESET_ALL)
                distances.append(dist)
                
    alldistances.append(distances)
    
            
    if args.solvfx:
        print("")
        print("                     ***************************")
        print("                     * Effect of the Solvation *")
        print("                     ***************************")
        print("")
        
        toten = []
        
        for n in range(1, n_copies+1):
            try:
                with open(f"complax_input_{n}solvent.xtbopt.xyz") as solv_en:
                    en = solv_en.readlines()[1].split()[1]
                    toten.append(en)
            except FileNotFoundError:
                print(f"⚠️ complax_input_{n}solvent.xtbopt.xyz not found. Skipping this configuration.")
                continue
        
        
        
        molB_file = open(f"{molBnoext}.xtbopt.xyz")
        molB_en = molB_file.readlines()[1].split()
        molB_en = float(molB_en[1])
        toten.insert(0,molB_en)
        molA_file = open(f"{molAnoext}.xtbopt.xyz")
        molA_en = molA_file.readlines()[1].split()
        molA_en = float(molA_en[1])    
        toten.insert(1,molA_en)
        
        toten = [float(str(i).strip()) for i in toten]
        [float(i) for i in toten]

        n = args.c  
        
        headers = [str(i) for i in range(1, n+1)]
        
        E_solvent = toten[0]
        E_reag = toten[1]
        
        somma_reag_solv = [E_reag + i * E_solvent for i in range(1, n+1)]

        calcolo_insieme = []
        for i in range(1, n+1):
            if 1 + i < len(toten):
                calcolo_insieme.append(toten[1 + i])
            else:
                calcolo_insieme.append("N/A")

        diff = []
        for theo, calc in zip(somma_reag_solv, calcolo_insieme):
            if isinstance(calc, float):
                diff.append((calc - theo)*627.51)
            else:
                diff.append("N/A")

        table = [
            ["Sum of reag + solv"] + somma_reag_solv,
            ["Tot Energy (Complex)"] + calcolo_insieme,
            ["ΔE (kcal/mol)"] + diff
        ]

        print(tabulate(table, headers=[""] + headers, tablefmt="grid", floatfmt=".6f"))
        
    print("")
    print("All the results have been saved in the 'outplax' folder.")
            

if __name__ == "__main__":
    main()
    exit()