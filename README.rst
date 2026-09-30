MMPBSA_Scoring - A python pipeline for estimating single residue mutation effects on protein-protein interactions
=================================================================================================================

Install the conda environment from the file ``environment.yml`` and execute the python script ``mmpbsa_scoring/ui.py`` to run the pipeline.

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


References / Licenses
---------------------

See also our preprint about the method: `bioRXiv <https://www.biorxiv.org/content/10.64898/2026.09.22.753561v1>`_

This project is licensed under the MIT license. Keep in mind that third-party software components required to run the tool are subject to their own license agreements.

**Important:** Users have to obtain a copy of FoldX and place it in the folder ``mmpbsa_scoring/bin`` as ``foldx``, be aware that FoldX is subject to its own `license agreements <https://foldxsuite.crg.eu/licensing-and-services>`_.