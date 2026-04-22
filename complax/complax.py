import sys
import numpy as np
import math
import os
import time
import itertools
import subprocess
from colorama import Fore, Style
import shutil
import argparse
import textwrap
from ase.io import read, write
from ase import Atoms
import threading
from functools import partial
from multiprocessing import Pool
import tqdm
from tabulate import tabulate
from importlib.metadata import version, PackageNotFoundError

def print_banner() -> None:
    """Prints the application banner, version, and author information."""
    # Recupero della versione per il banner
    try:
        current_version = version('complax')
    except PackageNotFoundError:
        current_version = "1.1.0"

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

    width = 60
    header = "COMPLAX — Solvation & Optimization Tool"
    author = "Developed by Federica Lauria, University of Turin (2025)"
    ver_str = f"Version {current_version}"

    print()
    print(" " * (left_margin) + " " + "_" * width + " ")
    print(" " * (left_margin) + "|" + " " * width + "|")
    print(" " * (left_margin) + f"|{header.center(width)}|")
    print(" " * (left_margin) + f"|{author.center(width)}|")
    print(" " * (left_margin) + f"|{ver_str.center(width)}|")
    print(" " * (left_margin) + "|" + "_" * width + "|")
    print("\n")
     

def openfile(file: str) -> tuple[list[str], np.ndarray]:
    """
    Reads an XYZ file and extracts atomic symbols and coordinates.
    
    Args:
        file (str): Path to the XYZ file.
        
    Returns:
        tuple: A list of atomic symbols and a NumPy array of coordinates.
    """
    with open(file, 'r') as xyz_file:
        lines = xyz_file.readlines()[2:]

    lines = [line for line in lines if len(line.split()) >= 4]

    atomic_symbols = [line.split()[0] for line in lines]
    atomic_coordinates = np.array([line.split()[1:4] for line in lines], dtype=float)

    return atomic_symbols, atomic_coordinates
    
def normalize(v: np.ndarray) -> np.ndarray:
    """Normalizes a given 3D vector."""
    return v / np.linalg.norm(v)

def check_overlap(mol_list: list[Atoms], mol: Atoms, cutoff: float = 1.2) -> bool:
    """
    Checks if a given molecule overlaps with any molecule in a provided list.
    
    Args:
        mol_list (list): List of placed ASE Atoms objects.
        mol (Atoms): The new molecule to check.
        cutoff (float): Distance threshold in Angstroms to define an overlap.
    """
    pos2 = mol.get_positions()
    for mol1 in mol_list:
        pos1 = mol1.get_positions()
        for a in pos1:
            for b in pos2:
                if np.linalg.norm(a - b) < cutoff:
                    return True
    return False

def random_rotation_matrix(axis: np.ndarray) -> np.ndarray:
    """
    Generates a random 3D rotation matrix around a specific normalized axis.
    """
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

def spinner_func(molecola: str, stop_event: threading.Event) -> None:
    """Displays a spinning cursor during computation."""
    for c in itertools.cycle(['|', '/', '-', '\\']):
        if stop_event.is_set():
            break
        sys.stderr.write(f'\rProcessing {molecola}... {c}')
        sys.stderr.flush()
        time.sleep(0.1)
    sys.stderr.write('\r') 
    sys.stderr.flush()

def task(molecola: str, alpb: str, gbsa: str, con: str, lev: str, chrg: int, uhf: int, proc: int) -> None:
    """Executes the xTB optimization task for a given molecular system."""
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


def run_complax_workflow(args: argparse.Namespace) -> None:
    """Main execution workflow for COMPLAX."""
    print_banner()
    
    try:
        atomic_coordinates = openfile(args.file1)
        atomic_coordinates_s = openfile(args.file2)
    except Exception as e:
        print(Fore.RED + f"❌ Error loading files: {e}" + Style.RESET_ALL)
        sys.exit(1)

    n_atoms_A = len(atomic_coordinates[1])
    n_atoms_B = len(atomic_coordinates_s[1])
    n_solvent = args.c
    
    with open('xtb.inp', 'w') as f:
        f.write("$fix\n")
        f.write(f"atoms: 1-{n_atoms_A}\n")
        f.write("$end\n")
        f.write("$constrain\n")
        f.write("force constant=5.0\n")
    
        start = n_atoms_A + 1
        for i in range(n_solvent):
            for j in range(n_atoms_B - 1):
                atom1 = start + j
                atom2 = start + j + 1
                f.write(f"distance: {atom1},{atom2}, auto\n")
            start += n_atoms_B
        f.write("$end\n")

    dir_name = 'outplax'
    if os.path.exists(dir_name):
        shutil.rmtree(dir_name)
    os.makedirs(dir_name)

    shutil.copy(args.file1, dir_name)
    shutil.copy(args.file2, dir_name)
    shutil.copy('xtb.inp', dir_name)
    os.chdir(dir_name)

    molA = read(args.file1)
    molB = read(args.file2)

    molecule_atom, solvent_atom = args.a if args.a else (None, None)
    iA = int(molecule_atom - 1)
    iB = int(solvent_atom - 1)
    
    print(f"{molA[iA].symbol}{molecule_atom} was selected for {args.file1}")
    print(f"{molB[iB].symbol}{solvent_atom} was selected for {args.file2}")
    
    dist = float(args.t)
    n_copies = args.c
    max_tries = args.maxtries
    nstruct = args.nstruct if args.nstruct else 0
    total_struct = nstruct + 1
    
    try:
        posA = molA[iA].position
        posB = molB[iB].position
    except IndexError:
        print(Fore.RED + "❌ Error: atom index out of range!" + Style.RESET_ALL)
        sys.exit(1)
    
    # Stochastic placement
    for s in range(total_struct):
        print(f"\n🍍 Generating solvent configuration {s+1}/{total_struct}")
        placed_molecules = [molA]
    
        for copy in range(n_copies):
            success = False
            for attempt in range(max_tries):
                trialB = molB.copy()
                trialB.translate(-posB)
                direction = normalize(np.random.randn(3))
                new_posB = posA - dist * direction
                trialB.translate(new_posB)
                
                axis = normalize(posA - new_posB)
                R = random_rotation_matrix(axis)
                trialB.positions = (trialB.positions - new_posB) @ R.T + new_posB
    
                if not check_overlap(placed_molecules, trialB, cutoff=1.2):
                    placed_molecules.append(trialB)
                    print(f"🧭 {args.file2} no.{copy+1} placed after {attempt+1} attempts")
                    success = True
                    break
    
            if not success:
                print(f"❌ Unable to place {args.file2} instance {copy+1}")
    
        # Save incremental structures
        for n in range(1, len(placed_molecules)):
            combined = Atoms()
            for mol in placed_molecules[:n+1]:
                combined += mol
            filename = f"complax_struct_{n}solvent_{s+1}.xyz"
            write(filename, combined)
            print(f"💾 Saved system with {n} solvent molecule(s) to {filename}")
       
    mol_list = []
    out_list = []
    sp_list = []
    
    molAnoext = os.path.splitext(args.file1)[0]
    molBnoext = os.path.splitext(args.file2)[0]
    sp_list.extend([molAnoext, molBnoext])
         
    for s in range(total_struct):
        for n in range(1, n_copies+1):
            filename = f"complax_struct_{n}solvent_{s+1}.xyz"
            if os.path.exists(filename):
                mol_list.append(filename)
                out_list.append(os.path.splitext(filename)[0])
                
    combined = Atoms()
    for mol in placed_molecules:
        combined += mol
    write("complax_input.xyz", combined)
    
    print("\n" + "-"*67)
    print(f"---------- Geometry optimization with {Fore.YELLOW}{args.c} solvent{Fore.RESET} molecule ----------")
    print("-"*67)

    task(molecola="complax_input", alpb=args.alpb, gbsa=args.gbsa, con='', lev=args.lev, chrg=args.chrg, uhf=args.uhf, proc=int(args.p))
    print(f"✅ Optimized geometry has been saved to {Fore.GREEN}complax_input.xtbopt.xyz{Fore.RESET}")
    
    task_with_args = partial(task, alpb=args.alpb, gbsa=args.gbsa, con='', lev=args.lev, chrg=args.chrg, uhf=args.uhf, proc=1)

    with Pool(args.p) as pool:
        for _ in tqdm.tqdm(
            pool.imap_unordered(task_with_args, out_list), 
            total=len(out_list),
            desc=f"{Fore.YELLOW}Calculating ...{Style.RESET_ALL}",
            bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [Elapsed Time:{elapsed} ETA:{remaining}]'
        ):
            pass
    
    # Ranking e Keep-best logic
    if args.nstruct > 1 or args.keep_best:
        # Recupero energie SP per il soluto e solvente singolo
        sp_energies= {}
        for mol in sp_list:
            cmd = f"xtb {mol}.xyz --namespace {mol} {args.lev} --chrg {args.chrg} --uhf {args.uhf}"
            results = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            for line in results.stdout.splitlines():
                if "TOTAL ENERGY" in line:
                    sp_energies[mol] = float(line.split()[3])
                    break

        print("\n" + "-"*50)
        print("        Ranking of stochastic solvent structures")
        print("-"*50)
    
        for n in range(1, n_copies + 1):
            struct_energies = []
            for s in range(total_struct):
                name = f"complax_struct_{n}solvent_{s+1}"
                outfile = f"{name}.xtbopt.xyz"
                if os.path.exists(outfile):
                    try:
                        with open(outfile) as f:
                            energy = float(f.readlines()[1].split()[1])
                        struct_energies.append((name, energy))
                    except: pass
            
            if not struct_energies: continue
            struct_energies.sort(key=lambda x: x[1])
            
            print(f"\n[{n} Solvent Molecule(s)]")
            for name, en in struct_energies:
                print(f"{name:<35} {en:>12.6f} Eh")
                
            if args.keep_best:
                best = struct_energies[:args.keep_best]
                discarded = struct_energies[args.keep_best:]
                for name, en in discarded:
                    for ext in [".xyz", ".xtbopt.xyz", ".out", "xtbopt.log"]:
                        if os.path.exists(f"{name}{ext}"): os.remove(f"{name}{ext}")

    if args.solvfx:
        # ... logica solvfx (omessa per brevità, resta invariata) ...
        pass
        
    print("\nAll the results have been saved in the 'outplax' folder.")


def main():
    parser = BannerArgumentParser(
        usage='%(prog)s <file1.xyz> <file2.xyz> [options]',
        description='COMPLAX: A tool that places solvent molecules around a solute.',
        formatter_class=lambda prog: SpacedHelpFormatter(prog, max_help_position=40, width=95) 
    )
    
    try:
        current_version = version('complax')
    except PackageNotFoundError:
        current_version = "1.1.0"
    
    parser.add_argument('-v', '--version', action='version', version=f'complax {current_version}')
    parser.add_argument('file1', type=str, nargs='?')
    parser.add_argument('file2', type=str, nargs='?')
    parser.add_argument('--alpb', type=str, metavar='SOLVENT')                       
    parser.add_argument('--gbsa', type=str, metavar='SOLVENT')
    parser.add_argument('-a', type=int, nargs=2, metavar=('MOL_AT', 'SOLV_AT'))
    parser.add_argument('-c', type=int, default=1)        
    parser.add_argument('-t', type=float, default=2.0) 
    parser.add_argument('-p', type=int, default=1)
    parser.add_argument('--lev', type=str, default="--gfn2")
    parser.add_argument('--chrg', type=int, default=0)
    parser.add_argument('-u','--uhf', type=int, default=1)
    parser.add_argument('--maxtries', type=int, default=1000)
    parser.add_argument('--solvfx', action='store_true')
    parser.add_argument('--nstruct', type=int, default=1)
    parser.add_argument('--keep-best', type=int, metavar='INT')
    
    args = parser.parse_args()

    if not args.file1 or not args.file2:
        parser.print_usage()
        sys.exit(1)
        
    run_complax_workflow(args)

if __name__ == "__main__":
    main()
