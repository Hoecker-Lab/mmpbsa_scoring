import os
import re
import subprocess
import numpy as np
from typing import Tuple, NoReturn
import logging

logger = logging.getLogger(__name__)


def parse_energies(output: str) -> Tuple[float]:
    """
    Parses solvation and apolar energy from a string.

    Parameters
    ----------
    output: str
        The output string from the APBS calculation

    Returns
    -------
    Energy values
    """
    energy_pattern_solvation = r'Global net ELEC energy = (-?\d+\.\d+)E.(\d+) kJ/mol'
    energy_pattern_apolar = r'Global net APOL energy = (-?\d+\.\d+)E.(\d+) kJ/mol'

    match_solvation = re.search(energy_pattern_solvation, output)
    match_apolar = re.search(energy_pattern_apolar, output)

    if match_solvation and match_apolar:
        match, exponent = match_solvation.groups()
        nrg_solvation = float(match) * 10 ** int(exponent)

        match, exponent = match_apolar.groups()
        nrg_apolar = float(match) * 10 ** int(exponent)

        return nrg_solvation, nrg_apolar
    else:
        logger.error(f'Could not parse energies, APBS calculation failed.')
        raise RuntimeError(f'Could not parse energies, APBS calculation failed.')


def create_apbs_input_file(input_pqr_file: str, grid_lengths: np.ndarray, grid_spacings: np.ndarray) -> NoReturn:
    """
    Create an input file for APBS solvation energy calculation.

    Parameters
    ----------
    input_pqr_file: str
        Name of the input PQR file
    grid_lengths: np.ndarray
        Number of grid points in x, y and z
    grid_lengths: np.ndarray
        Grid spacings for x, y and z
    """
    input_string = [f"""read
    mol pqr {input_pqr_file}
end

elec name solv
    mg-manual
    mol 1
    gcent 0 0 0
    grid {grid_spacings[0]} {grid_spacings[1]} {grid_spacings[2]}
    dime {grid_lengths[0]} {grid_lengths[1]} {grid_lengths[2]}
    lpbe
    bcfl mdh
    pdie 4.0
    sdie 80.0
    chgm spl2
    srfm smol
    swin 0.3
    sdens 10.0
    srad 1.4
    temp 300.0
    ion charge 1 conc 0.15 radius 2.0
    ion charge -1 conc 0.15 radius 2.0
    calcforce no
    calcenergy total
end

elec name vac
    mg-manual
    mol 1
    gcent 0 0 0
    grid {grid_spacings[0]} {grid_spacings[1]} {grid_spacings[2]}
    dime {grid_lengths[0]} {grid_lengths[1]} {grid_lengths[2]}
    lpbe
    bcfl mdh
    pdie 4.0
    sdie 4.0
    chgm spl2
    srfm smol
    swin 0.3
    sdens 10.0
    srad 1.4
    temp 300.0
    calcforce no
    calcenergy total
end

apolar name area
    grid 0.1 0.1 0.1
    mol 1
    srfm sacc
    swin 0.3
    srad 1.4
    press 0.0
    gamma 1.0
    bconc 0.0
    sdens 10.0
    dpos 0.2
    temp 300.0
    calcforce no
    calcenergy total
end

print elecEnergy solv - vac end 
print apolEnergy area end

quit"""]

    with open(os.path.splitext(os.path.basename(input_pqr_file))[0] + '.in', 'w') as file:
        for line in input_string:
            file.write(line + '\n')


def run_apbs(input_pqr_file: str, grid_lengths: np.ndarray, grid_spacings: np.ndarray) -> float:
    """
    Calculates solvation and apolar energy of a structure in PQR format.

    Parameters
    ----------
    input_pqr_file: str
        Name of the input PQR file
    grid_lengths: np.ndarray
        Number of grid points in x, y and z
    grid_lengths: np.ndarray
        Grid spacings for x, y and z

    Returns
    -------
    Solvation and apolar energy of the structure in kcal/mol.
    """

    create_apbs_input_file(input_pqr_file, grid_lengths, grid_spacings)
    output = subprocess.run(['apbs', os.path.splitext(os.path.basename(input_pqr_file))[0] + '.in'], check=True, capture_output=True)
    output = output.stdout.decode('utf-8')
    solvation_energy, apolar_energy = parse_energies(output)

    return solvation_energy / 4.184, apolar_energy / 4.184