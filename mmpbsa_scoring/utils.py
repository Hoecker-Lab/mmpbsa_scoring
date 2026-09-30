import subprocess
import numpy as np
from Bio.Data import IUPACData
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
from Bio.PDB.vectors import rotmat, Vector
from pdbfixer import PDBFixer
from openmm.app import *
from openmm import *
from openmm.unit import *
from copy import deepcopy
from typing import List, Tuple, NoReturn
import logging

logger = logging.getLogger(__name__)

foldx_bin = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bin', 'foldx')

if not os.path.isfile(foldx_bin):
    logger.error(f'Cannot locate FoldX binary. Please make sure it is located at the following path: {foldx_bin}.')
    raise FileNotFoundError(f'Cannot locate FoldX binary. Please make sure it is located at the following path: {foldx_bin}.')


def cif2pdb(cif_file: str, pdb_file: str) -> NoReturn:
    """
    Converts a CIF (mmCIF) file into a PDB file

    Parameters
    ----------
    cif_file: str
        Path to the input mmCIF file
    pdb_file: str
        Path to the output PDB file
    """
    cif_parser = MMCIFParser(QUIET=True)

    structure = cif_parser.get_structure('structure', cif_file)

    pdb_writer = PDBIO()
    pdb_writer.set_structure(structure)

    pdb_writer.save(pdb_file)


def clean(input_pdb_file: str, output_pdb_file: str, keep_chains: List[str]) -> NoReturn:
    """
    Cleans input PDB file and removes unwanted chains.

    Parameters
    ----------
    input_pdb_file: str
        Path to the input PDB file
    output_pdb_file : str
        Path to the output PDB file
    keep_chains : list of str
        List of chain IDs to keep in the output file
    """
    class CleanSelect(Select):
        def accept_residue(self, residue):
            return residue.id[0] == ' ' and residue.get_parent().id in keep_chains

        def accept_atom(self, atom):
            return atom.element != 'H'

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_pdb_file)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_file, select=CleanSelect())


def repack(input_pdb_file: str, output_pdb_file: str) -> NoReturn:
    """
    Repacks an input pdb structure using FoldX

    Parameters
    ----------
    input_pdb_file: str
        Path to the input pdb file
    output_pdb_file: str
        Path to the output pdb file
    """
    command = [foldx_bin, '--command=RepairPDB', f'--pdb={input_pdb_file}']
    subprocess.run(command, check=True, capture_output=True)
    os.rename(os.path.splitext(input_pdb_file)[0] + '_Repair.pdb', output_pdb_file)


def write_mutation_file(input_pdb_file: str, mutation: Tuple[str, int, str]) -> NoReturn:
    """
    Writes a mutation file for a single mutation.

    Parameters
    ----------
    input_pdb_file : str
        Path to the input PDB file
    mutation : tuple
        Mutation in the form (chain_id, residue_id, mutated_resname) where mutated_resname is 3-letter code.
    """
    three_to_one = IUPACData.protein_letters_3to1

    chain_id, res_id, mut_3letter = mutation

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_pdb_file)
    model = structure[0]

    try:
        chain = model[chain_id]
        residue = chain[(' ', res_id, ' ')]
        orig_3letter = residue.get_resname()
    except KeyError:
        logger.error(f'Residue {res_id} not found in chain {chain_id}.')
        raise ValueError(f'Residue {res_id} not found in chain {chain_id}.')

    try:
        orig_1letter = three_to_one[orig_3letter.capitalize()]
        mut_1letter = three_to_one[mut_3letter.capitalize()]
    except KeyError:
        logger.error(f'Invalid residue name: {orig_3letter} or {mut_3letter}')
        raise ValueError(f'Invalid residue name: {orig_3letter} or {mut_3letter}')

    mutation_str = f'{orig_1letter}{chain_id}{res_id}{mut_1letter}'

    with open('individual_list.txt', 'w') as file:
        file.write(mutation_str + ';')


def mutate(input_pdb_file: str, output_mut_pdb_file: str, output_wt_pdb_file: str) -> NoReturn:
    """
    Builds mutated and corresponding wildtype model of input PDB structure using FoldX.

    Parameters
    ----------
    input_pdb_file: str
        Path to the input PDB file
    output_mut_pdb_file: str
        Path to the mutated output PDB file
    output_wt_pdb_file: str
        Path to the wildtype output PDB file
    """
    command = [foldx_bin, '--command=BuildModel', f'--pdb={input_pdb_file}', '--mutant-file=individual_list.txt', '--order=_USERDEFINED']
    subprocess.run(command, check=True, capture_output=True)

    os.rename(os.path.splitext(input_pdb_file)[0] + '_1.pdb', output_mut_pdb_file)
    os.rename('WT_' + os.path.splitext(input_pdb_file)[0] + '_1.pdb', output_wt_pdb_file)


def calculate_stability(input_pdb_file: str) -> float:
    """
    Calculate protein stability using FoldX.

    Parameters
    ----------
    input_pdb_file: str
        Path to the input PDB file

    Returns
    -------
    Protein stability in kcal/mol.
    """
    command = [foldx_bin, '--command=Stability', f'--pdb={input_pdb_file}']
    output = subprocess.run(command, check=True, capture_output=True)

    stability = 0.0
    for line in output.stdout.decode('ASCII').splitlines():
        if line.startswith('Total'):
            stability = float(line.split()[2])
            break

    return stability


def minimize(input_pdb_file: str, output_pdb_file: str) -> NoReturn:
    """
    Minimize structure with backbone restraints.

    Parameters
    ----------
    input_pdb_file : str
        Path to the input PDB file
    output_pdb_file : str
        Path to the output PDB file
    """
    fixer = PDBFixer(filename=input_pdb_file)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=0)
    PDBFile.writeFile(fixer.topology, fixer.positions, open('tmp.pdb', 'w'), keepIds=True)

    command = ['reduce', '-BUILD', '-NUClear', '-NOCon', 'tmp.pdb']
    output = subprocess.run(command, capture_output=True)

    with open('tmp_protonated.pdb', 'w') as file:
        for line in output.stdout.decode('ASCII').splitlines():
            file.write(line + '\n')
    
    fixer = PDBFixer(filename='tmp_protonated.pdb')
    fixer.addMissingHydrogens(pH=7.0)
    
    forcefield = ForceField('amber19-all.xml', 'implicit/obc2.xml')
    system = forcefield.createSystem(fixer.topology, nonbondedMethod=NoCutoff, constraints=HBonds)

    # Restrain backbone atoms
    restraint = CustomExternalForce('0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)')
    restraint.addGlobalParameter('k', 10.0 * kilocalories_per_mole / angstroms ** 2)
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')

    for i, atom in enumerate(fixer.topology.atoms()):
        if atom.name in ('N', 'CA', 'C', 'O'):
            pos = fixer.positions[i]
            restraint.addParticle(i, [pos.x, pos.y, pos.z])

    system.addForce(restraint)

    platform = Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed'}
    integrator = LangevinIntegrator(300.0 * kelvin, 1.0 / picosecond, 2.0 * femtoseconds)

    simulation = Simulation(fixer.topology, system, integrator, platform, properties)
    simulation.context.setPositions(fixer.positions)

    simulation.minimizeEnergy()

    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(output_pdb_file, 'w'), keepIds=True)


def calc_nonbonded_nrgs(input_pdb_file: str, epsilon_protein: float = 4.0) -> Tuple[float, float]:
    """
    Read in a PDB file using OpenMM and calculate coulomb and vdw energy.

    Parameters
    ----------
    input_pdb_file: str
        Path to the input PDB file
    epsilon_protein: float
        Protein dielectric constant [default: 4.0]

    Returns
    -------
    Coulomb and vdw energy in kcal/mol.
    """
    pdb = PDBFile(input_pdb_file)

    forcefield = ForceField('amber19-all.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=None)

    coulomb = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
    vdw = deepcopy(coulomb)
    for i in range(coulomb.getNumParticles()):
        charge, sigma, epsilon = coulomb.getParticleParameters(i)
        coulomb.setParticleParameters(i, charge, 1.0, 0.0)
        vdw.setParticleParameters(i, 0.0, sigma, epsilon)
    for i in range(coulomb.getNumExceptions()):
        p1, p2, chargeProd, sigma, epsilon = coulomb.getExceptionParameters(i)
        coulomb.setExceptionParameters(i, p1, p2, chargeProd, 1.0, 0.0)
        vdw.setExceptionParameters(i, p1, p2, 0.0, sigma, epsilon)
    coulomb.setForceGroup(1)
    vdw.setForceGroup(2)
    system.addForce(vdw)

    integrator = LangevinIntegrator(300.0 * kelvin, 1.0 / picosecond, 2.0 * femtoseconds)
    context = Context(system, integrator)
    context.setPositions(pdb.positions)

    ecoulomb = context.getState(getEnergy=True, groups={1}).getPotentialEnergy().value_in_unit(kilocalorie_per_mole) / epsilon_protein
    evdw = context.getState(getEnergy=True, groups={2}).getPotentialEnergy().value_in_unit(kilocalorie_per_mole)

    return ecoulomb, evdw


def pdb2pqr(input_pdb_file: str, output_pqr_file: str) -> NoReturn:
    """
    Convert a PDB file into a PQR file:

    * Assign AMBER partial atomic charges and modified radii

    Parameters
    ----------
    input_pdb_file: str
        Path to the input PDB file
    output_pqr_file: str
        Path to the output PQR file
    """
    # Rename H3 of N-terminal PRO to enable correct recognition by pdb2pqr
    lines = []
    with open(input_pdb_file, 'r') as file:
        for line in file:
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                res_name = line[17:20]
                if res_name == 'PRO' and atom_name == 'H3':
                    line = line.replace('H3', 'H ')
            lines.append(line)

    with open('tmp.pdb', 'w') as file:
        file.writelines(lines)

    command = ['pdb2pqr', 'tmp.pdb', output_pqr_file, '--keep-chain', '--assign-only', '--whitespace',
               '--userff', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'parameters', 'AMBER_modified_radii.DAT'),
               '--usernames', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'parameters', 'AMBER_modified_radii.names')]
    subprocess.run(command, check=True, capture_output=True)


def extract_chains(input_pdb_file: str, output_pdb_file: str, extract_chain_ids: List[str]) -> NoReturn:
    """
    Extract specific chains from a PDB file and save them to a new PDB file.

    Parameters
    ----------
    input_pdb_file : str
        Path to the input PDB file
    output_pdb_file : str
        Path to the output PDB file
    extract_chain_ids : list
        Chain IDs to extract (e.g., ['A', 'B']).
    """
    class ChainSelect(Select):
        """Select only specified chains."""

        def __init__(self, chain_ids: List[str]):
            self.chain_ids = set(chain_ids)

        def accept_chain(self, chain):
            return chain.id in self.chain_ids

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_pdb_file)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_file, ChainSelect(extract_chain_ids))


def align(input_pdb_file: str, output_pdb_file: str) -> NoReturn:
    """
    Align entire PDB structure on x and y-axis and move to origin based on center of geometry.

    Parameters
    ----------
    input_pdb_file : str
        Path to the input PDB file
    output_pdb_file : str
        Path to the output PDB file
    """
    def pca(coords: np.ndarray) -> np.ndarray:
        """
        Calculate PCA of coordinates.

        Parameters
        ----------
        coords : np.ndarray
            Coordinates to calculate PCA

        Returns
        -------
        PCA
        """
        mean = coords.mean(axis=0)
        centered = coords - mean
        cov = np.cov(centered.T)
        eigvals, eigvecs = np.linalg.eigh(cov)
        order = np.argsort(eigvals)[::-1]
        return eigvecs[:, order]
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_pdb_file)

    center = structure.center_of_mass(geometric=True)
    structure.transform(np.identity(3), -center)

    coords = np.array([atom.get_coord() for atom in structure.get_atoms()])
    components = pca(coords)

    R1 = np.array(rotmat(Vector(components[:, 0]), Vector(1, 0, 0))).T
    structure.transform(R1, [0, 0, 0])

    coords = np.array([atom.get_coord() for atom in structure.get_atoms()])
    components = pca(coords)

    R2 = np.array(rotmat(Vector(components[:, 1]), Vector(0, 1, 0))).T
    structure.transform(R2, [0, 0, 0])

    center = structure.center_of_mass(geometric=True)
    structure.transform(np.identity(3), -center)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_file)


def get_grid_lengths(input_pdb_file: str, grid_spacings: np.ndarray, padding_distances: np.ndarray = np.array([10.0, 10.0, 10.0])) -> np.ndarray:
    """
    Calculate lengths in x, y and z of grid based on maximum dimensions of PDB file, grid spacings and padding distances.

    Parameters
    ----------
    input_pdb_file: str
        Path to the input PDB file
    grid_spacings: np.ndarray
        Grid spacing in each dimension
    padding_distances: np.ndarray
        Padding distance in each dimension [default: np.array([10.0, 10.0, 10.0])]

    Returns
    -------
    Number of grid points in each dimension
    """
    def next_multiple_64_plus_1(value: float) -> int:
        """
        Calculates the next multiple 64 plus 1

        Parameters
        ----------
        value: float
            Input value

        Returns
        -------
        Next multiple 64 plus 1
        """
        n = (value - 1) // 64
        candidate = n * 64 + 1

        if candidate < value:
            candidate = (n + 1) * 64 + 1

        return candidate

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_pdb_file)
    model = structure[0]

    coordinates = []
    for chain in model:
        for residue in chain:
            for atom in residue:
                coordinates.append(atom.get_coord())

    min_values = np.min(coordinates, axis=0) - padding_distances
    max_values = np.max(coordinates, axis=0) + padding_distances

    next_multiple_64_plus_1_vec = np.vectorize(next_multiple_64_plus_1)
    glen = next_multiple_64_plus_1_vec(np.abs(max_values - min_values) / grid_spacings)

    return glen