ELECTRONS_DICT = {'adaptive_thr': {'default': '.FALSE',
                  'description': 'Adaptive threshold for electronic convergence',
                  'info': ' If .TRUE. this turns on the use of an adaptive '
                          '@ref conv_thr for the inner scf loops when using '
                          'EXX. ',
                  'input_type': 'select_multiple',
                  'options': ['.FALSE.', '.TRUE.'],
                  'type': 'LOGICAL'},
 'conv_thr': {'default': '1.D-6',
              'description': 'Convergence threshold',
              'info': ' Convergence threshold for selfconsistency: estimated '
                      'energy error < conv_thr (note that conv_thr is '
                      'extensive, like the total energy).  For '
                      'non-self-consistent calculations, conv_thr is used to '
                      'set the default value of the threshold (ethr) for '
                      'iterative diagonalization: see @ref diago_thr_init ',
              'type': 'REAL'},
 'conv_thr_init': {'default': '1.D-3',
                   'description': 'Initial convergence threshold',
                   'info': ' When @ref adaptive_thr = .TRUE. this is the '
                           'convergence threshold used for the first scf '
                           'cycle. ',
                   'type': 'REAL'},
 'conv_thr_multi': {'default': '1.D-1',
                    'description': 'Convergence threshold multiplier',
                    'info': ' When @ref adaptive_thr = .TRUE. the convergence '
                            'threshold for each scf cycle is given by: max( '
                            '@ref conv_thr, @ref conv_thr_multi * dexx ) ',
                    'type': 'REAL'},
 'diago_cg_maxiter': {'description': 'Maximum number of CG iterations for '
                                     'diagonalization',
                      'info': ' For conjugate gradient diagonalization:  max '
                              'number of iterations ',
                      'type': 'INTEGER'},
 'diago_david_ndim': {'default': '2',
                      'description': 'Number of Davidson dimensions for '
                                     'diagonalization',
                      'info': ' For Davidson diagonalization: dimension of '
                              'workspace (number of wavefunction packets, at '
                              'least 2 needed). A larger value may yield a '
                              'smaller number of iterations in the algorithm '
                              'but uses more memory and more CPU time in '
                              'subspace diagonalization (cdiaghg/rdiaghg). You '
                              'may try @ref diago_david_ndim=4 if you are not '
                              'tight on memory and if the time spent in '
                              'subspace diagonalization is small compared to '
                              'the time spent in h_psi ',
                      'type': 'INTEGER'},
 'diago_full_acc': {'default': '.FALSE.',
                    'description': 'Full accuracy flag for diagonalization',
                    'info': ' If .TRUE. all the empty states are diagonalized '
                            'at the same level of accuracy of the occupied '
                            'ones. Otherwise the empty states are diagonalized '
                            'using a larger threshold (this should not affect '
                            'total energy, forces, and other ground-state '
                            'properties). ',
                    'input_type': 'select_multiple',
                    'options': ['.FALSE.', '.TRUE.'],
                    'type': 'LOGICAL'},
 'diago_gs_nblock': {'default': '16',
                     'description': 'Number of block iterations for GS '
                                    'diagonalization',
                     'info': ' For RMM-DIIS diagonalization: blocking size of '
                             'Gram-Schmidt orthogonalization ',
                     'type': 'INTEGER'},
 'diago_ppcg_maxiter': {'description': 'Maximum number of PPCG iterations for '
                                       'diagonalization',
                        'info': ' For @b ppcg diagonalization:  max number of '
                                'iterations ',
                        'type': 'INTEGER'},
 'diago_rmm_conv': {'default': '.FALSE.',
                    'description': 'Convergence threshold for residual '
                                   'minimization method diagonalization',
                    'info': ' If .TRUE., RMM-DIIS is performed up to converge. '
                            'If .FALSE., RMM-DIIS is performed only once. ',
                    'input_type': 'select_multiple',
                    'options': ['.FALSE.', '.TRUE.'],
                    'type': 'LOGICAL'},
 'diago_rmm_ndim': {'default': '4',
                    'description': 'Number of RMM iterations for '
                                   'diagonalization',
                    'info': ' For RMM-DIIS diagonalization: dimension of '
                            'workspace (number of wavefunction packets, at '
                            'least 2 needed). ',
                    'type': 'INTEGER'},
 'diago_thr_init': {'description': 'Initial diagonalization threshold',
                    'info': ' Convergence threshold (ethr) for iterative '
                            'diagonalization (the check is on eigenvalue '
                            'convergence).  For scf calculations: default is '
                            '1.D-2 if starting from a superposition of atomic '
                            'orbitals; 1.D-5 if starting from a charge '
                            'density. During self consistency the threshold is '
                            'automatically reduced (but never below 1.D-13) '
                            'when approaching convergence.  For non-scf '
                            'calculations: default is (@ref conv_thr/N '
                            'elec)/10. ',
                    'type': 'REAL'},
 'diagonalization': {'default': 'david',
                     'description': 'Diagonalization algorithm',
                     'info': 'Available options are:'
                        '\'david\' : Davidson iterative diagonalization with overlap matrix'
                        '(default). Fast, may in some rare cases fail.'
                        '\'cg\' : Conjugate-gradient-like band-by-band diagonalization.' 
                        'MUCH slower than \'david\' but uses less memory and is (a little bit) more robust.'
                        '\'ppcg\' : PPCG iterative diagonalization'
                        '\'paro\', \'ParO\' : ParO iterative diagonalization'
                        '\'rmm-davidson\', \'rmm-paro\' : RMM-DIIS iterative diagonalization.'
                        'To stabilize the SCF loop'
                        'RMM-DIIS is alternated with calls to Davidson or'
                        'ParO  solvers depending on the string used.'
                        'Other variables that can be used to tune the behavior of'
                        'RMM-DIIS are:  diago_rmm_ndim and diago_rmm_conv',
                     'input_type': 'select_multiple',
                     'options': ['david', 'cg','ppcg','paro','ParO','rmm-davidson','rmm-paro'],
                     'type': 'CHARACTER'},
 'efield': {'default': '0.D0',
            'description': 'Electric field amplitude',
            'info': ' Amplitude of the finite electric field (in Ry a.u.; 1 '
                    'a.u. = 36.3609*10^10 V/m). Used only if @ref '
                    'lelfield==.TRUE. and if k-points (@ref K_POINTS card) are '
                    'not automatic. ',
            'type': 'REAL'},
 'efield_phase': {'default': 'none',
                  'description': 'Phase of the electric field',
                  'info': '',
                  'input_type': 'select_multiple',
                  'options': ['read', 'write', 'none'],
                  'type': 'CHARACTER'},
 'electron_maxstep': {'default': '100',
                      'description': 'Maximum electron steps',
                      'info': ' maximum number of iterations in a scf step. If '
                              'exact exchange is active, this will affect the '
                              'inner loops. ',
                      'type': 'INTEGER'},
 'exx_maxstep': {'default': '100',
                 'description': 'Maximum number of steps for EXX '
                                'diagonalization',
                 'info': ' maximum number of outer iterations in a scf '
                         'calculation with exact exchange. ',
                 'type': 'INTEGER'},
 'mixing_beta': {'default': '0.7D0',
                 'description': 'Mixing parameter',
                 'info': ' mixing factor for self-consistency ',
                 'type': 'REAL'},
 'mixing_fixed_ns': {'default': '0',
                     'description': 'Fixed number of states for mixing',
                     'info': ' For DFT+U : number of iterations with fixed ns '
                             '( ns is the atomic density appearing in the '
                             'Hubbard term ). ',
                     'type': 'INTEGER'},
 'mixing_mode': {'default': 'plain',
                 'description': 'Mixing mode',
                 'info': '',
                 'input_type': 'select_multiple',
                 'options': ['plain', 'TF', 'local-TF'],
                 'type': 'CHARACTER'},
 'mixing_ndim': {'default': '8',
                 'description': 'Number of mixing iterations',
                 'info': ' number of iterations used in mixing scheme. If you '
                         'are tight with memory, you may reduce it to 4 or so. '
                         '',
                 'type': 'INTEGER'},
 'real_space': {'default': '.FALSE.',
                'description': 'Real space diagonalization',
                'info': ' If .true., exploit real-space localization to '
                        'compute matrix elements for nonlocal projectors. '
                        'Faster and in principle better scaling than the '
                        'default G-space algorithm, but numerically less '
                        'accurate, may lead to some loss of translational '
                        'invariance. Use with care and after testing! ',
                'input_type': 'select_multiple',
                'options': ['.FALSE.', '.TRUE.'],
                'type': 'LOGICAL'},
 'scf_must_converge': {'default': '.TRUE.',
                       'description': 'SCF must converge flag',
                       'info': ' If .false. do not stop molecular dynamics or '
                               'ionic relaxation when electron_maxstep is '
                               'reached. Use with care. ',
                       'input_type': 'select_multiple',
                       'options': ['.FALSE.', '.TRUE.'],
                       'type': 'LOGICAL'},
 'startingpot': {'description': 'Starting potentials',
                 'info': '',
                 'input_type': 'select_multiple',
                 'options': ['atomic', 'file'],
                 'type': 'CHARACTER'},
 'startingwfc': {'default': 'atomic+random',
                 'description': 'Starting wavefunctions',
                 'info': '',
                 'input_type': 'select_multiple',
                 'options': ['atomic', 'atomic+random', 'random', 'file'],
                 'type': 'CHARACTER'},
 'tqr': {'default': '.FALSE.',
         'description': 'TQR algorithm',
         'info': ' If .true., use a real-space algorithm for augmentation '
                 'charges of ultrasoft pseudopotentials and PAWsets. Faster '
                 'but numerically less accurate than the default G-space '
                 'algorithm. Use with care and after testing! ',
         'input_type': 'select_multiple',
         'options': ['.FALSE.', '.TRUE.'],
         'type': 'LOGICAL'}}
