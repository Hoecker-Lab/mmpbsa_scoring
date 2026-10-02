MMPBSA_Scoring - A python pipeline for estimating single residue mutation effects on protein-protein interactions
=================================================================================================================

Install the conda environment from the file ``environment.yml`` and execute the python script ``mmpbsa_scoring/ui.py`` to run the pipeline.

**Important**: Edit the line ``cuda-version=x.x`` in the file ``environment.yml`` to make sure the correct version of CUDA is installed for your system!

.. code-block:: bash

  conda env create -f environment.yml


The pipeline consists of:

1. Structural modelling using FoldX

2. Energy minimization using OpenMM

3. Binding free energy calculation using OpenMM, APBS and FoldX

**Important**: First prepare the starting structure in mmCIF or PDB format by using the flag ``--prepare`` and use the prepared output model (``structure_prepared.pdb``) as input for all mutational modelling.

.. code-block:: bash

    Pipeline for estimating single residue mutation effects on protein-protein interactions.

    options:
      -h, --help            show this help message and exit
      --structure STRUCTURE
                            Input structure file in mmCIF or PDB format
      --prepare             Prepare input structure instead of evaluating mutation
      --chain_group_a CHAIN_GROUP_A
                            Chain IDs for chain group A or first interaction partner (For example: A)
      --chain_group_b CHAIN_GROUP_B
                            Chain IDs for chain group B or second interaction partner (For example: B,C)
      --mutation MUTATION   Specific mutation to evaluate, in the format: Chain_ID:Residue_ID:Residue_Name (For example: A:15:ALA)
      --tmp_dir TMP_DIR     Path to temporary directory (Will be deleted afterwards!)

Example
-------

See also the provided example in the folder ``example``, you can recreate the output by executing:

.. code-block:: bash

    python3 mmpbsa_scoring/ui.py --prepare --structure 2FTL.pdb --chain_group_a E --chain_group_b I

    2026-10-02 11:00:48,375 - mmpbsa_scoring - INFO - Cleaning structure.
    2026-10-02 11:00:48,397 - mmpbsa_scoring - INFO - Repacking structure.

.. code-block:: bash

    python3 mmpbsa_scoring/ui.py --structure structure_prepared.pdb --chain_group_a E --chain_group_b I --mutation I:15:ASP

    2026-10-02 12:03:08,875 - mmpbsa_scoring - INFO - Mutating structure.
    2026-10-02 12:03:35,064 - mmpbsa_scoring - INFO - Minimizing mutated structure.
    2026-10-02 12:03:38,998 - mmpbsa_scoring - INFO - Aligning mutated structure.
    2026-10-02 12:03:39,072 - mmpbsa_scoring - INFO - Calculating binding free energy of mutated structure.
    2026-10-02 12:04:10,341 - mmpbsa_scoring - INFO - Minimizing wildtype structure.
    2026-10-02 12:04:13,817 - mmpbsa_scoring - INFO - Aligning wildtype structure.
    2026-10-02 12:04:13,892 - mmpbsa_scoring - INFO - Calculating binding free energy of wildtype structure.
    2026-10-02 12:04:44,580 - mmpbsa_scoring - INFO - Estimated binding free energy difference: 9.458 kcal/mol.
    2026-10-02 12:04:44,580 - mmpbsa_scoring - INFO - Estimated change in VdW energy: 0.966 kcal/mol.
    2026-10-02 12:04:44,580 - mmpbsa_scoring - INFO - Estimated change in electrostatic energy: 7.225 kcal/mol.
    2026-10-02 12:04:44,580 - mmpbsa_scoring - INFO - Estimated change in apolar solvation energy: -0.422 kcal/mol.
    2026-10-02 12:04:44,580 - mmpbsa_scoring - INFO - Estimated change in folding stability: 1.194 kcal/mol.

**Important:** Be aware that there might be small differences in the results, which are due to the minimization!

References / Licenses
---------------------

See also our preprint about the method: `bioRXiv <https://www.biorxiv.org/content/10.64898/2026.09.22.753561v1>`_

This project is licensed under the MIT license. Keep in mind that third-party software components required to run the tool are subject to their own license agreements.

**Important:** Users have to obtain a copy of FoldX and place the binary in the folder ``mmpbsa_scoring/bin`` under the name ``foldx``, be aware that FoldX is subject to its own `license agreements <https://foldxsuite.crg.eu/licensing-and-services>`_.