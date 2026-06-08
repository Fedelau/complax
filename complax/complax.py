from __future__ import annotations
import sys
import numpy as np
import math
import os
import time
import itertools
import subprocess
from colorama import Fore, Style, init
import shutil
import argparse
from ase.io import read, write
from ase import Atoms
import threading
from functools import partial
from multiprocessing import Pool
import tqdm
from tabulate import tabulate
from importlib.metadata import version, PackageNotFoundError
import glob

init(autoreset=True)

def print_banner() -> None:
    """Prints the application banner, version, and author information."""
    
    try:
        current_version = version('complax')
    except PackageNotFoundError:
        current_version = "1.1.1"

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
    author = "Developed by Federica Lauria and Andrea Maranzana"
    affiliation = "University of Turin (2026)"
    ver_str = f"Version {current_version}"

    print()
    print(" " * left_margin + " " + "_" * width + " ")
    print(" " * left_margin + "|" + " " * width + "|")
    print(" " * left_margin + f"|{header.center(width)}|")
    print(" " * left_margin + f"|{author.center(width)}|")
    print(" " * left_margin + f"|{affiliation.center(width)}|")
    print(" " * left_margin + f"|{ver_str.center(width)}|")
    print(" " * left_margin + "|" + "_" * width + "|")
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

def check_overlap(mol_list: list[Atoms], mol: Atoms, cutoff: float = 1.8) -> bool:
    """
    Checks if a given molecule overlaps with any molecule in a provided list.
    Vectorized steric clash check using numpy broadcasting.
    """
    pos2 = mol.get_positions()
    cutoff_sq = cutoff ** 2  
    
    for mol1 in mol_list:
        pos1 = mol1.get_positions()
        
        diff = pos1[:, np.newaxis, :] - pos2[np.newaxis, :, :]
        dist_sq = np.sum(diff ** 2, axis=-1)
        
        if np.any(dist_sq < cutoff_sq):
            return True
    return False

def rotation_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
    """
    Generates a 3D rotation matrix around a specific normalized axis by a given angle theta.
    """
    axis = normalize(axis)
    c, s = np.cos(theta), np.sin(theta)
    x, y, z = axis
    R = np.array([
        [c+(1-c)*x*x,     (1-c)*x*y - s*z, (1-c)*x*z + s*y],
        [(1-c)*y*x + s*z, c+(1-c)*y*y,     (1-c)*y*z - s*x],
        [(1-c)*z*x - s*y, (1-c)*z*y + s*x, c+(1-c)*z*z    ]
    ])
    return R

def fibonacci_sphere(samples: int) -> list[np.ndarray]:
    """
    Generates uniformly distributed points on the surface of a unit sphere 
    using the Fibonacci spiral method.
    """
    points = []
    phi = math.pi * (3.0 - math.sqrt(5.0))  
    
    for i in range(samples):
        y = 1.0 - (i / float(samples - 1)) * 2.0  
        radius = math.sqrt(1.0 - y * y) 
        
        theta = phi * i 
        
        x = math.cos(theta) * radius
        z = math.sin(theta) * radius
        
        points.append(np.array([x, y, z]))
        
    return points

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

def task(molecola: str, alpb: str, gbsa: str, gbe: str, con: str, lev: str, chrg: int, uhf: int, proc: int) -> None:
    """Executes the xTB optimization task for a given molecular system."""
    global spinner_used

    use_spinner = False
    if not spinner_used:
        spinner_used = True
        use_spinner = True
        
    if use_spinner:
        stop_event = threading.Event()
        t = threading.Thread(target=spinner_func, args=(molecola, stop_event))
        t.start()

    # --- PRE-OPTIMIZATION ---
    preopt_cmd = f"xtb {molecola}.xyz --opt --cycles 50 --gfn1 --input xtb_preopt.inp --namespace {molecola}_preopt --chrg {chrg} --uhf {uhf} -P {proc}"
    subprocess.run(preopt_cmd, shell=True, capture_output=True, text=True)
    
    preopt_file = f"{molecola}_preopt.xtbopt.xyz"
    if os.path.exists(preopt_file) and os.path.getsize(preopt_file) > 0:
        start_geom = preopt_file
    else:
        start_geom = f"{molecola}.xyz"

    # --- REAL OPTIMIZATION ---
    cmd_parts = [
        f"xtb {start_geom}",    
        f"--opt --{lev.lstrip('-')}",            
        f"--input xtb_main.inp",   
        f"--namespace {molecola}", 
        f"--chrg {chrg}",
        f"--uhf {uhf}",          
        f"-P {proc}"               
    ]
    
    if alpb:
        cmd_parts.insert(3, f"--alpb {alpb}")
    elif gbsa:
        cmd_parts.insert(3, f"--gbsa {gbsa}")
    elif gbe:
        cmd_parts.insert(3, f"--gbe {gbe}")

    cmd = " ".join(cmd_parts)

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
    

    molecule_atom = args.a[0] if args.a else None
    solvent_atom = args.a[1] if args.a else None

    with open('xtb_preopt.inp', 'w') as f:
        f.write("$fix\n")
        f.write(f"atoms: 1-{n_atoms_A}\n")
        f.write("$end\n")

        if molecule_atom and solvent_atom:
            f.write("$constrain\n")
            f.write("force constant=5.0\n")
            for i in range(n_solvent):
                target_solv_atom = n_atoms_A + (i * n_atoms_B) + solvent_atom
                f.write(f"distance: {molecule_atom}, {target_solv_atom}, auto\n")
        start = n_atoms_A + 1
        for i in range(n_solvent):
            for j in range(n_atoms_B - 1):
                atom1 = start + j
                atom2 = start + j + 1
                f.write(f"distance: {atom1},{atom2}, auto\n")
            start += n_atoms_B
        f.write("$end\n")

    with open('xtb_main.inp', 'w') as f:
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
    shutil.copy('xtb_preopt.inp', dir_name)
    shutil.copy('xtb_main.inp', dir_name)
    os.chdir(dir_name)

    molA = read(args.file1)
    molB = read(args.file2)

    iA = int(molecule_atom - 1)
    iB = int(solvent_atom - 1)
    
    print(f"{molA[iA].symbol}{molecule_atom} was selected for {args.file1}")
    print(f"{molB[iB].symbol}{solvent_atom} was selected for {args.file2}")
    
    dist = float(args.t)
    n_copies = args.c
    grid_points = args.res 
    total_struct = args.nstruct 
    
    try:
        posA = molA[iA].position
        posB = molB[iB].position
    except IndexError:
        print(Fore.RED + "❌ Error: atom index out of range!" + Style.RESET_ALL)
        sys.exit(1)
    
    if args.seed != -1:
        np.random.seed(args.seed)
    else:
        np.random.seed(None)
    
    # Placement
    for s in range(total_struct):
        print(f"\n🍍 Generating solvent configuration {s+1}/{total_struct}")
        placed_molecules = [molA]
        for copy in range(n_copies):
            success = False

            fib_points = fibonacci_sphere(grid_points)
            np.random.shuffle(fib_points)

            rot_steps = 36  # rotazioni per punto spaziale

            for attempt in range(grid_points):
                direction = fib_points[attempt]
                new_posB = posA - dist * direction

                for _ in range(rot_steps):
                    trialB = molB.copy()
                    trialB.translate(-posB)

                    random_axis = normalize(np.random.randn(3))
                    random_angle = np.random.uniform(0, 2 * np.pi)
                    R = rotation_matrix(random_axis, random_angle)
                    trialB.positions = trialB.positions @ R.T

                    trialB.translate(new_posB)

                    if not check_overlap(placed_molecules, trialB, cutoff=args.cutoff):
                        placed_molecules.append(trialB)
                        print(f"🧭 {args.file2} no.{copy+1} placed after {attempt+1} spatial points")
                        success = True
                        break

                if success:
                    break

            if not success:
                print(Fore.RED + f"❌ Unable to place {args.file2} instance {copy+1}. Steric crowding is too high." + Style.RESET_ALL)
    
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
                
    # combined = Atoms()
    # for mol in placed_molecules:
    #     combined += mol
    # write("complax_input.xyz", combined)
    
    total_to_optimize = len(out_list)
    print("\n" + "-"*67)
    print(f"---------- Geometry optimization of {Fore.YELLOW}{total_to_optimize}{Fore.RESET} structure(s) ----------------")
    print("-"*67)

    # task(molecola="complax_input", alpb=args.alpb, gbsa=args.gbsa, gbe=args.gbe, con='', lev=args.lev, chrg=args.chrg, uhf=args.uhf, proc=int(args.p))
    # print(f"✅ Optimized geometry has been saved to {Fore.GREEN}complax_input.xtbopt.xyz{Fore.RESET}")
    
    task_with_args = partial(task, alpb=args.alpb, gbsa=args.gbsa, gbe=args.gbe, con='', lev=args.lev, chrg=args.chrg, uhf=args.uhf, proc=1)

    with Pool(args.p) as pool:
        for _ in tqdm.tqdm(
            pool.imap_unordered(task_with_args, out_list), 
            total=len(out_list),
            desc=f"{Fore.YELLOW}Calculating ...{Style.RESET_ALL}",
            bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [Elapsed Time:{elapsed} ETA:{remaining}]'
        ):
            pass

    sp_energies = {}
    for mol in sp_list:
        solv_cmd = ""
        if args.alpb: solv_cmd = f"--alpb {args.alpb}"
        elif args.gbsa: solv_cmd = f"--gbsa {args.gbsa}"
        elif args.gbe: solv_cmd = f"--gbe {args.gbe}"

        cmd = f"xtb {mol}.xyz --namespace {mol} --{args.lev.lstrip('-')} {solv_cmd} --chrg {args.chrg} --uhf {args.uhf}"
        results = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        for line in results.stdout.splitlines():
            if "TOTAL ENERGY" in line:
                sp_energies[mol] = float(line.split()[3])
                break

    # --- Ranking e Keep-best logic ---
    if args.nstruct > 1 or args.keep_best:
        ranking_output = "\n" + "-"*50 + "\n"
        ranking_output += "        Ranking of stochastic solvent structures\n"
        ranking_output += "-"*50 + "\n"
    
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

            ranking_output += f"\n[{n} Solvent Molecule(s)]\n"
            for name, en in struct_energies:
                ranking_output += f"{name:<35} {en:>12.6f} Eh\n"
                
            if args.keep_best:
                best = struct_energies[:args.keep_best]
                discarded = struct_energies[args.keep_best:]
                for name, en in discarded:
                    for file_to_delete in glob.glob(f"{name}*"):
                        if os.path.exists(file_to_delete):
                            os.remove(file_to_delete)

        print(ranking_output)

        with open("complax_summary.txt", "a", encoding="utf-8") as f_out:
            f_out.write(ranking_output)

    if args.solvfx:
        solvfx_output = "\n                     ***************************\n"
        solvfx_output += "                     * Effect of the Solvation *\n"
        solvfx_output += "                     ***************************\n\n"

        best_energies_per_n = []
        for n in range(1, n_copies + 1):
            min_en = None
            for s in range(total_struct):
                try:
                    with open(f"complax_struct_{n}solvent_{s+1}.xtbopt.xyz") as solv_en:
                        en = float(solv_en.readlines()[1].split()[1])
                        if min_en is None or en < min_en:
                            min_en = en
                except FileNotFoundError:
                    continue
            best_energies_per_n.append(min_en if min_en is not None else "N/A")

        E_solvent = sp_energies.get(molBnoext)
        E_reag = sp_energies.get(molAnoext)
       

        if E_solvent is None or E_reag is None:
            err_msg = "❌ Error: Could not find reference energies.\n"
            print(Fore.RED + err_msg + Style.RESET_ALL)
            with open("complax_summary.txt", "a", encoding="utf-8") as f_out:
                f_out.write(err_msg)
        else:
            headers = [str(i) for i in range(1, n_copies + 1)]
            somma_reag_solv = [E_reag + i * E_solvent for i in range(1, n_copies + 1)]
            
            diff = []
            for theo, calc in zip(somma_reag_solv, best_energies_per_n):
                if isinstance(calc, float):
                    diff.append((calc - theo) * 627.51)
                else:
                    diff.append("N/A")

            table = [
                ["Sum of reag + solv (Hartree)"] + [f"{x:.6f}" for x in somma_reag_solv],
                ["Tot Energy (Complex) (Hartree)"] + [f"{x:.6f}" if isinstance(x, float) else x for x in best_energies_per_n],
                ["ΔE (kcal/mol)"] + [f"{x:.2f}" if isinstance(x, float) else x for x in diff]
            ]

            table_str = tabulate(table, headers=[""] + headers, tablefmt="grid", disable_numparse=True)
            solvfx_output += table_str + "\n"
            print(solvfx_output)

            with open("complax_summary.txt", "a", encoding="utf-8") as f_out:
                f_out.write(solvfx_output)
        
    print("\nAll the results have been saved in the 'outplax' folder.")


def main():
    parser = BannerArgumentParser(
        usage='%(prog)s <file1.xyz> <file2.xyz> [options]',
        description='COMPLAX: A tool for automated microsolvation and geometry optimization using xTB.',
        epilog="For further information, please contact the programme author.",
        formatter_class=lambda prog: SpacedHelpFormatter(prog, max_help_position=40, width=95) 
    )
    
    try:
        current_version = version('complax')
    except PackageNotFoundError:
        current_version = "1.1.1"
    
    parser.add_argument('-v', '--version', action='version', version=f'complax {current_version}')
    # --- Positional Arguments ---
    parser.add_argument('file1', type=str, nargs='?', help='XYZ file of the solute molecule.')
    parser.add_argument('file2', type=str, nargs='?', help='XYZ file of the solvent molecule.')

    # --- Solvation Options ---
    parser.add_argument('-c', type=int, default=1, metavar="INT", help='Number of solvent molecules to add (default: 1).')
    parser.add_argument('-t', type=float, default=2.0, metavar="FLOAT", help='Target distance between atoms for placement (default: 2.0 Å).')
    parser.add_argument('-a', type=int, nargs=2, metavar=('MOLECULE_ATOM', 'SOLVENT_ATOM'), help='Indices of the two atoms of molecule and solvent files to keep at distance -t. The number has to be specified using 1-based indexing.')
    parser.add_argument('--res', type=int, default=10000, metavar="INT", help='Number of points in the Fibonacci spherical grid (default: 10000).')
    
    # --- xTB Options ---
    parser.add_argument('--alpb', type=str, metavar='[SOLVENT]', choices=['acetone', 'acetonitrile', 'aniline', 'benzaldehyde', 'benzene','ch2cl2', 'chcl3', 'cs2', 'dioxane', 'dmf', 'dmso', 'ether','ethylacetate', 'furane', 'hexandecane', 'hexane', 'methanol','nitromethane', 'octanol', 'woctanol', 'phenol', 'toluene', 'thf', 'water'], help='Analytical linearized Poisson-Boltzmann (ALPB) model, available solvents are acetone, acetonitrile, aniline, benzaldehyde, benzene, ch2cl2, chcl3, cs2, dioxane, dmf, dmso, ether, ethylacetate, furane, hexandecane, hexane, methanol, nitromethane, octanol, woctanol, phenol, toluene, thf, water.')
    parser.add_argument('--gbsa', type=str, metavar='[SOLVENT]', choices=['acetone', 'acetonitrile', 'benzene', 'CH2Cl2', 'CHCl3', 'CS2', 'DMF', 'DMSO', 'ether', 'H2O', 'methanol', 'n-hexan', 'THF', 'toluene','ch2cl2', 'chcl3', 'cs2', 'dmf', 'dmso', 'h2o', 'thf'], help='Generalized Born model with a simple switching function (GBSA), available solvents are acetone, acetonitrile, benzene (only GFN1-xTB), ch2cl2, chcl3, cs2, dmf (only GFN2-xTB), dmso, ether, h2o, methanol, n-hexane (only GFN2-xTB), thf and toluene.')
    parser.add_argument('--gbe', type=str, metavar='[SOLVENT]', help='Generalized Born model with finite epsilon for solvation. Requires --lev gxtb.')
    parser.add_argument('--lev', type=str, metavar='[METHOD]', default="gfn2", choices=['gfn0', 'gfn1', 'gfn2', 'gfnff', 'gxtb'], help='xTB optimization level (default: gfn2). Other options include gfn0, gfn1, gfn2, gfnff, gxtb. E.g.: --lev gfn1')
    parser.add_argument('--chrg', type=int, default=0, help='Total system charge (default: 0).')
    parser.add_argument('-u', '--uhf', type=int, default=0, help='Number of unpaired electrons (default: 0).')
    parser.add_argument('-p', type=int, default=1, help='Number of processors for parallel execution (default: 1).')

    # --- Analysis and Ranking ---
    parser.add_argument('--nstruct', type=int, default=1, help='Number of different stochastic structures to generate.')
    parser.add_argument('--keep-best', type=int, metavar='N', help='Keep only the N lowest-energy optimized structures.')
    parser.add_argument('--solvfx', action='store_true', help='Calculate and display the solvation energy effect table.')
    parser.add_argument('--cutoff', type=float, default=1.6, metavar="FLOAT", help='Steric overlap distance threshold in Å (default: 1.6).')
    parser.add_argument('--seed', type=int, default=42, help='Seed for the pseudo-random number generator to ensure reproducibility. Set to -1 for pure stochasticity (default: 42).')
    
    args = parser.parse_args()
    
    if args.gbe and args.lev != "gxtb":
        print(Fore.RED + " Error: The --gbe model is only compatible with --lev gxtb." + Style.RESET_ALL)
        sys.exit(1)
        
    if args.lev == "gxtb":
        if args.alpb or args.gbsa:
            print(Fore.RED + " Error: When using --lev gxtb, only --gbe is supported for solvation." + Style.RESET_ALL)
            print(Fore.YELLOW + " Please remove --alpb or --gbsa and use --gbe if needed." + Style.RESET_ALL)
            sys.exit(1)

    if not args.file1 or not args.file2:
        parser.print_help()
        sys.exit(1)
        
    run_complax_workflow(args)

if __name__ == "__main__":
    main()