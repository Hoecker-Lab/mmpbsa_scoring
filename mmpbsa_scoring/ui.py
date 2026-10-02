import os
import numpy as np
import shutil
import logging.config
from typing import List, Tuple, NoReturn

logging.config.fileConfig(os.path.join(os.path.dirname(__file__), 'logging.ini'), disable_existing_loggers=False)
logger = logging.getLogger('mmpbsa_scoring')

weights = np.array([0.206, 0.254, 0.022, 0.277])
intercept = 0.495


def run_preparation(structure: str, chain_group_a: List[str], chain_group_b: List[str]) -> NoReturn:
    """
    Prepares input structure by cleaning and repacking it.
    
    Parameters
    ----------
    structure : str
        Input structure file in mmCIF or PDB format
    chain_group_a: list
        List of chain IDs for chain group A
    chain_group_b: list
        List of chain IDs for chain group B
    """
    from mmpbsa_scoring.utils import cif2pdb, clean, repack

    if structure.endswith('.cif'):
        shutil.copy(structure, 'structure.cif')
        cif2pdb('structure.cif', 'structure.pdb')
    elif structure.endswith('.pdb'):
        shutil.copy(structure, 'structure.pdb')
    else:
        logger.error(f'Structural datatype not supported.')
        raise RuntimeError(f'Structural datatype not supported.')

    logger.info('Cleaning structure.')
    clean('structure.pdb', 'structure_cleaned.pdb', keep_chains=chain_group_a + chain_group_b)

    logger.info('Repacking structure.')
    repack('structure_cleaned.pdb', 'structure_cleaned_repacked.pdb')


def run_scoring(structure: str, chain_group_a: List[str], chain_group_b: List[str], mutation: Tuple[str, int, str]) -> NoReturn:
    """
    Calculates binding free energy difference for mutation between chain groups in input structure.

    Parameters
    ----------
    structure : str
        Prepared input structure file in PDB format
    chain_group_a: list
        List of chain IDs for chain group A
    chain_group_b: list
        List of chain IDs for chain group B
    mutation: tuple
        Specific mutation to evaluate
    """
    from mmpbsa_scoring.utils import write_mutation_file, mutate, minimize, calc_nonbonded_nrgs, align, pdb2pqr, extract_chains, get_grid_lengths, calculate_stability
    from mmpbsa_scoring.apbs import run_apbs
    
    if structure.endswith('.pdb'):
        shutil.copy(structure, 'structure.pdb')
    else:
        logger.error(f'Structural datatype not supported.')
        raise RuntimeError(f'Structural datatype not supported.')

    logger.info('Mutating structure.')
    write_mutation_file('structure.pdb', mutation)
    mutate('structure.pdb', 'structure_mutated.pdb', 'structure_wildtype.pdb')

    for structure_type in ['mutated', 'wildtype']:
        logger.info(f'Minimizing {structure_type} structure.')
        minimize(f'structure_{structure_type}.pdb', f'structure_{structure_type}_minimized.pdb')

        logger.info(f'Aligning {structure_type} structure.')
        align(f'structure_{structure_type}_minimized.pdb', f'structure_{structure_type}_minimized_aligned.pdb')

        logger.info(f'Calculating binding free energy of {structure_type} structure.')
        extract_chains(f'structure_{structure_type}_minimized_aligned.pdb',f'chainA_{structure_type}_minimized_aligned.pdb', extract_chain_ids=chain_group_a)
        extract_chains(f'structure_{structure_type}_minimized_aligned.pdb',f'chainB_{structure_type}_minimized_aligned.pdb', extract_chain_ids=chain_group_b)

        coulomb_complex, vdw_complex = calc_nonbonded_nrgs(f'structure_{structure_type}_minimized_aligned.pdb')
        coulomb_chainA, vdw_chainA = calc_nonbonded_nrgs(f'chainA_{structure_type}_minimized_aligned.pdb')
        coulomb_chainB, vdw_chainB = calc_nonbonded_nrgs(f'chainB_{structure_type}_minimized_aligned.pdb')

        pdb2pqr(f'structure_{structure_type}_minimized_aligned.pdb',f'structure_{structure_type}_minimized_aligned.pqr')
        pdb2pqr(f'chainA_{structure_type}_minimized_aligned.pdb',f'chainA_{structure_type}_minimized_aligned.pqr')
        pdb2pqr(f'chainB_{structure_type}_minimized_aligned.pdb',f'chainB_{structure_type}_minimized_aligned.pqr')

        grid_spacings = np.array([0.5, 0.5, 0.5])
        grid_lengths = get_grid_lengths(f'structure_{structure_type}_minimized_aligned.pdb', grid_spacings)

        solvation_complex, apolar_complex = run_apbs(f'structure_{structure_type}_minimized_aligned.pqr', grid_lengths, grid_spacings)
        solvation_chainA, apolar_chainA = run_apbs(f'chainA_{structure_type}_minimized_aligned.pqr', grid_lengths, grid_spacings)
        solvation_chainB, apolar_chainB = run_apbs(f'chainB_{structure_type}_minimized_aligned.pqr', grid_lengths, grid_spacings)

        vdw = vdw_complex - vdw_chainA - vdw_chainB
        solvation = solvation_complex - solvation_chainA - solvation_chainB
        coulomb = coulomb_complex - coulomb_chainA - coulomb_chainB
        apolar = apolar_chainA + apolar_chainB - apolar_complex
        stability = calculate_stability(f'structure_{structure_type}.pdb')

        if structure_type == 'mutated':
            scores = np.array([vdw, coulomb + solvation, apolar, stability])
        else:
            scores -= np.array([vdw, coulomb + solvation, apolar, stability])

    scores *= weights
    scores_rounded = np.round(scores, 3)
    ddG_bind = np.round(np.sum(scores) + intercept, 3)

    logger.info(f'Estimated binding free energy difference: {ddG_bind} kcal/mol.')
    logger.info(f'Estimated change in VdW energy: {scores_rounded[0]} kcal/mol.')
    logger.info(f'Estimated change in electrostatic energy: {scores_rounded[1]} kcal/mol.')
    logger.info(f'Estimated change in apolar solvation energy: {scores_rounded[2]} kcal/mol.')
    logger.info(f'Estimated change in folding stability: {scores_rounded[3]} kcal/mol.')


def main():
    import sys
    import tempfile
    import argparse

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

    parser = argparse.ArgumentParser(description='Pipeline for estimating single residue mutation effects on protein-protein interactions.')
    parser.add_argument('--structure', type=str, help='Input structure file in mmCIF or PDB format', required=True)
    parser.add_argument('--prepare', action='store_true', help='Prepare input structure instead of evaluating mutation')
    parser.add_argument('--chain_group_a', type=str, help='Chain IDs for chain group A or first interaction partner (For example: A)', required=True)
    parser.add_argument('--chain_group_b', type=str, help='Chain IDs for chain group B or second interaction partner (For example: B,C)', required=True)
    parser.add_argument('--mutation', type=str, help='Specific mutation to evaluate, in the format: Chain_ID:Residue_ID:Residue_Name (For example: A:15:ALA)', required=False)
    parser.add_argument('--tmp_dir', type=str, help='Path to temporary directory (Will be deleted afterwards!)', required=False)
    args = parser.parse_args()

    structure = os.path.abspath(args.structure)
    if not os.path.isfile(structure):
        logger.error(f'No structure file found: {structure}.')
        raise FileNotFoundError(f'No structure file found: {structure}.')

    if args.tmp_dir:
        tmp_dir = args.tmp_dir
    else:
        tmp_dir = os.path.join(tempfile.gettempdir(), 'mmpbsa_scoring')

    os.makedirs(tmp_dir, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(tmp_dir)
    os.environ['TMPDIR'] = tmp_dir

    if not args.chain_group_a or not args.chain_group_b:
        logger.error('Chain group A or B are not defined.')
        raise RuntimeError('Chain group A or B are not defined.')
    else:
        chain_group_a = args.chain_group_a.split(',')
        chain_group_b = args.chain_group_b.split(',')
    
    if args.prepare:
        run_preparation(structure=structure,
                        chain_group_a=chain_group_a,
                        chain_group_b=chain_group_b)

        if os.path.isfile('structure_cleaned_repacked.pdb'):
            shutil.copy('structure_cleaned_repacked.pdb', os.path.join(cwd, 'structure_prepared.pdb'))
        else:
            logger.error('Structure preparation failed.')
            raise RuntimeError('Structure preparation failed.')
        
    else:
        if args.mutation:
            allowed_res_names = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                                 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                                 'THR', 'TRP', 'TYR', 'VAL']
            try:
                chain_id, res_id, res_name = args.mutation.split(':')
                res_name = res_name.upper()
                if res_name not in allowed_res_names:
                    logger.error(f'Amino acid {res_name} not in allowed amino acids: {allowed_res_names}.')
                    raise RuntimeError(f'Amino acid {res_name} not in allowed amino acids: {allowed_res_names}.')
                try:
                    res_id = int(res_id)
                except ValueError:
                    logger.error('Residue ID should be a number.')
                    raise ValueError('Residue ID should be a number.')
            except ValueError:
                logger.error('Define mutation in the following format: Chain_ID:Residue_ID:Residue_Name.')
                raise ValueError('Define mutation in the following format: Chain_ID:Residue_ID:Residue_Name.')
    
            mutation = (chain_id, res_id, res_name)
        else:
            logger.error('No mutation to evaluate specified.')
            raise RuntimeError('No mutation to evaluate specified.')
    
        run_scoring(structure=structure,
                   chain_group_a=chain_group_a,
                   chain_group_b=chain_group_b,
                   mutation=mutation)

        shutil.copy('structure_mutated_minimized.pdb', os.path.join(cwd, 'structure_mutated.pdb'))
        shutil.copy('structure_wildtype_minimized.pdb', os.path.join(cwd, 'structure_wildtype.pdb'))

    os.chdir(cwd)
    shutil.rmtree(tmp_dir)


if __name__ == '__main__':
    main()