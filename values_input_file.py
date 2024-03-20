CONTROL_VALUES = {
    'calculation': {
        'values': ['scf', 'nscf', 'bands', 'relax', 'md', 'vc-relax', 'vc-md'],
        'type': 'Character',
        'mandatory': True,
        'input_type' : 'select_multiple'
    },
    'title': {
        'values': [''],
        'type': 'Character',
        'mandatory': False,
        'input_type' : 'text'
    },
    'verbosity': {
        'values': ['low', 'high'],
        'type': 'Character',
        'mandatory': False,
        'input_type' : 'select_multiple'
    },
    'restart_mode': {
        'values': ['from_scratch', 'restart'],
        'type': 'Character',
        'mandatory': False,
        'input_type' : 'select_multiple'
    },
    # 'wf_collect': {
    #     'values': [True],
    #     'type': None,  # No longer implemented
    #     'mandatory': False
    # },
    'nstep': {
        'values': ['1','50'], #1 if calculation == 'scf', 'nscf', 'bands'; 50 for the other cases 
        'type': 'Integer',
        'mandatory': True,
        'input_type' : 'automatic'
    },
    'iprint': {
        'values': [''],
        'type': 'Integer',
        'mandatory': False,
        'input_type' : 'automatic'
    },
    'tstress': {
        'values': ['.false.', '.true.'],
        'type': 'Logical',
        'mandatory': False,
        'input_type' : 'select_multiple'
    },
    'tprnfor': {
        'values': ['.true.', '.false'], #if calculation == 'relax', 'md', 'vc-md' is .TRUE.
        'type': 'Logical',
        'mandatory': False,
        'input_type' : 'automatic'
    },
    'dt': {
        'values': ['20.D0'],
        'type': 'Real',
        'mandatory': False,
        'input_type' : 'text'
    },
    'outdir' : {
        'values' : ['./'], #value of the ESPRESSO_TMPDIR environment variable if set
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'wfcdir' : {
        'values' : ['./'], #same as outdir
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'automatic'
    },
    'prefix' : {
        'values' : ['pwscf'],
        'type' : 'Character',
        'mandatory' : True,
        'input_type' : 'text'
    },
    'max_seconds' : {
        'values' : ['1.D+7'],
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'etot_conv_thr' : {
        'values' : ['1.0D-4'],
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'forc_conv_thr' : {
        'values' : ['1.0D-3'],
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'disk_io' : {
        'values' : ['see below','high','medium','low','nowf','none'],
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'pseudo_dir' : {
        'values' : ['$ESPRESSO_PSEUDO', '$HOME/espresso/pseudo/'],
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'tefield' : {
        'values' : ['.FALSE.', '.TRUE.'],
        'type' : 'Logical',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'dipfield' : {
        'values' : ['.FALSE.', '.TRUE.'],
        'type' : 'Logical',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'lelfield' : {
        'values' : ['.FALSE.', '.TRUE.'],
        'type' : 'Logical',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'nberrycyc' : {
        'values' : ['1'], # if lelfield == .TRUE.
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'lorbm' : {
        'values' : ['.FALSE.','.TRUE.'], # if lelfield == .TRUE.
        'type' : 'Logical',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'lberry' : {
        'values' : ['.FALSE.','.TRUE.'],
        'type' : 'Logical',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'gdir' : {
        'values' : ['1','2','3'], # if lelfield == .TRUE.
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'nppstr' : {
        'values' : [''], # if lelfield == .TRUE.
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'gate' : {
        'values' : ['.FALSE.','.TRUE.'], # if dipfield=.true.
        'type' : 'Logical',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'twochem' : {
        'values' : ['.FALSE.','.TRUE.'],
        'type' : 'Logical',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'lfcp' : {
        # if calculation = relax or md
        # assume_isolated='esm' and 'esm_bc'='bc2' or 'bc3'
        # ignore_wolfe is always .TRUE. for BFGS
        'values' : ['.FALSE.','.TRUE.'],
        'type' : 'Logical',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'trism' : {
        'values' : ['.FALSE.','.TRUE.'],
        'type' : 'Logical',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    }
}



SYSTEM_VALUES = {
    'ibrav' : {
        'values' : ['0','1','2','3','-3','4','5',
                    '-5','6','7','8','9','-9','91',
                    '10','11','12','-12','13','-13','14'],
        'type' : 'Integer',
        'mandatory' : True,
        'input_type' : 'select_multiple'
    },
    # ========= AGREGAR celldmi(i) o A, B, C, cosAB, cosAC, cosBC =========
    'nat' : {
        'values' : [''],
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'ntyp' : {
        'values' : [''],
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'nbnd' : {
        'values' : [''], # for an insulator, nbnd = number of valence bands (nbnd = # of electrons /2); for a metal, 20% more (minimum 4 more) 
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'nbnd_cond' : {
        'values' : [''], #nbnd_cond = nbnd - # of electrons / 2 in the collinear case; nbnd_cond = nbnd - # of electrons in the noncollinear case. 
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'tot_charge' : {
        'values' : ['0.0'],
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'starting_charge' : {
        'values' : ['0.0'], # Different (i), starting_charge(i), i=1,ntyp
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'tot_magnetization' : {
        'values' : ['-10000'],
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'starting_magnetization' : {
        'values' : ['0'], # Differente (i), starting_magnetization(i), i=1,ntyp
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'ecutwfc' : {
        'values' : ['0'],
        'type' : 'Real',
        'mandatory' : True,
        'input_type' : 'text'
    },
    'ecutrho' : {
        'values' : ['4'], # IMPORTANT 4*ecutwfc
        'type' : 'Real',
        'mandatory' : True,
        'input_type' : 'text'
    },
    'ecutfock' : {
        'values' : [''], # Default: Same as ecutrho
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'nr1' : {
        'values' : [''],
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'nr2' : {
        'values' : [''],
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'nr3' : {
        'values' : [''],
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'nr1s' : {
        'values' : [''], # Coincides with nr1, nr2, nr3 if ecutrho = 4 * ecutwfc ( default )
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'nr2s' : {
        'values' : [''], # Coincides with nr1, nr2, nr3 if ecutrho = 4 * ecutwfc ( default )
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'nr3s' : {
        'values' : [''], # Coincides with nr1, nr2, nr3 if ecutrho = 4 * ecutwfc ( default )
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'text'
    }

    # ========== TERMINAR DE ESCRIBIR ==========
}


ELECTRONS_VALUES = {
    'electron_maxstep' : {
        'values' : ['100'],
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'exx_maxstep' : {
        'values' : ['100'],
        'type' : 'Integer',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'scf_must_converge' : {
        'values' : ['.TRUE.', '.FALSE.'],
        'type' : 'Logical',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'conv_thr' : {
        'values' : ['1.D-6'], # estimated energy error < conv_thr
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'adaptive_thr' : {
        'values' : ['.FALSE.', '.TRUE.'],
        'type' : 'Logical',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'conv_thr_init' : {
        'values' : ['1.D-3'],
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'conv_thr_multi' : {
        'values' : ['1.D-1'], # When adaptive_thr = .TRUE. the convergence threshold for each scf cycle is given by: max( conv_thr, conv_thr_multi * dexx )
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'automatic'
    },
    'mixing_mode' : {
        'values' : ['plain', 'TF', 'local-TF'],
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    }
    # =============== TERMINAR ESTA SECCION ===========
}

# REQUIRED if calculation == 'relax', 'md', 'vc-relax', or 'vc-md' OPTIONAL for calculation == 'scf' (only ion_positions is used) 
IONS_VALUES = {
    'ion_positions' : {
        'values': ['default', 'from_input'],
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'ion_velocities' : {
        'values' : ['default', 'from_input'],
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'ion_dynamic' : {
        'values' : [['bfgs', 'damp', 'fire'], # if calculation == 'relax'
                    ['verlet', 'langevin', 'langevin-smc'], # if calculation == 'md'
                    ['bfgs', 'damp'], # if calculation == 'vc-relax'
                    ['beeman'] # if calculation == 'vc-md'
                    ],
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'pot_extrapolation' : {
        'values' : ['atomic', 'none', 'first_order', 'second_order'],
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'wfc_extrapolation' : {
        'values' : ['none', 'first_order', 'second_order'],
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'remove_grid_rot' : {
        'values' : ['.FALSE.', '.TRUE.'],
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    # ======== variables used for molecular dynamics =======
    'ion_temperature' : {
        'values' : ['not_controlled', 'rescaling', 'rescale-v', 'rescale-T',
                    'berendsen', 'andersen', 'svr', 'initial'],
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    }
    # =============== TERMINAR ESTA SECCION ===========
}

# input this namelist only if calculation == 'vc-relax' or 'vc-md' 
CELL_VALUES = {
    'cell_dynamics' : {
        'values' : [['none', 'sd', 'damp-pr', 'damp-w', 'bfgs'], # if calculation == 'vc-relax
                    ['none', 'pr', 'w'] # if calculation = 'vc-md'
                    ],
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'press' : {
        'values' : ['0.D0'],
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'wmass' : {
        'values' : [''], # 0.75*Tot_Mass/pi**2 for Parrinello-Rahman MD; 0.75*Tot_Mass/pi**2/Omega**(2/3) for Wentzcovitch MD
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'cell_factor' : {
        'values' : ['2.0', '1.0'], # 2.0 for variable-cell calculations, 1.0 otherwise
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'automatic'
    },
    'press_conv_thr' : {
        'values' : ['0.5D0 Kbar'],
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'all' : {
        'values' : ['all', 'ibrav', 'a', 'b', 'c', 'fixa', 'fixb', 'fixc',
                    'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz', 'shape', 'volume',
                    '2Dxy', '2Dshape', 'epitaxial_ab', 'epitaxial_ac', 'epitaxial_bc'],
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    }
}

# Input this namelist only if lfcp = .TRUE. 
FSP_VALUES = {
    'fcp_mu' : {
        'values' : [''],
        'type' : 'Character',
        'mandatory' : True,
        'input_type' : 'text'
    },
    'fcp_dynamics' : {
        'values' : [['bfgs', 'newton', 'damp', 'lm'], # if calculation == 'relax'
                    ['velocity-verlet', 'verlet'] # if calculation == 'md'
                    ],
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'select_multiple'
    },
    'fcp_conv_thr' : {
        'values' : ['1.D-2'],
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'fcp_ndiis' : { # used only if dcp_dynamics == 'newton'
        'values' : ['4'], 
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    # ============ Variables used for FCP dynamics. ===========
    'fcp_mass' : {
        'values' : [''], # 5.D+6 / (xy area) for ESM only; 5.D+4 / (xy area) for ESM-RISM
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'fcp_velocity' : {
        'values' : [''], # Determined by fcp_temperature
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'fcp_temperature' : {
        'values' : ['rescaling', 'rescale-v', 'rescale-T', 'reduce-T', 'berendsen',
                    'andersen', 'initial', 'not_controlled'], # Default : ion_temperature
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'fcp_tempw' : {
        'values' : [''], # Default : tempw
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    },
    'fcp_tolp' : {
        'values' : [''], # Default : tolp
        'type' : 'Real',
        'mandatory' : False,
        'input_type' : 'text'
    }

    # ================ TERMINAR =============
}

# Input this namelist only if trism = .TRUE. 
RISM_VALUES = {
    'nsolv' : {
        'values' : [''],
        'type' : 'Integer',
        'mandatory' : True,
        'input_type' : 'text'
    },
    'closure' : {
        'values' : ['kh', 'hnc'],
        'type' : 'Character',
        'mandatory' : False,
        'input_type' : 'text'
    }
    # ================ TERMINAR =============
}