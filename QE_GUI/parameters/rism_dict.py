RISM_DICT = {'closure': {'default': 'kh',
             'description': 'Closure relation for the Ornstein-Zernike '
                            'equation',
             'info': ' Specify the type of closure equation: }',
             'input_type': 'select_multiple',
             'options': ['kh', 'hnc'],
             'type': 'CHARACTER'},
 'ecutsolv': {'default': '4 * @ref ecutwfc',
              'description': 'Cutoff energy for the solvent',
              'info': " Kinetic energy cutoff (Ry) for solvent's correlation "
                      'functions. If a solute is an isolated system or slab, '
                      'you may allowed to use default value. For a frameworked '
                      'or porous solute (e.g. Zeolite, MOF), it is desirable '
                      'to apply a larger value. Solvents confined in a '
                      'framework often have a high frequency. }',
              'type': 'REAL'},
 'laue_both_hands': {'default': '.FALSE.',
                     'description': 'Use both hands in the Laue method',
                     'info': ' If .TRUE., you can set different densities to '
                             'the solvent regions of right-hand side and '
                             'left-hand side. See @ref SOLVENTS card. }',
                     'input_type': 'select_multiple',
                     'options': ['.FALSE.', '.TRUE.'],
                     'type': 'LOGICAL'},
 'laue_buffer_left': {'description': 'Left buffer for the Laue method',
                      'info': ' If positive value, set the buffering length '
                              '[in a.u.] of the solvent region on left-hand '
                              'side of the unit cell. Then correlation '
                              'functions are defined inside of [ -L_z/2 - @ref '
                              'laue_expand_left , @ref laue_starting_left + '
                              '@ref laue_buffer_left ]. This is only for '
                              'Laue-RISM. }',
                      'type': 'REAL'},
 'laue_buffer_right': {'description': 'Right buffer for the Laue method',
                       'info': ' If positive value, set the buffering length '
                               '[in a.u.] of the solvent region on right-hand '
                               'side of the unit cell. Then correlation '
                               'functions are defined inside of [ @ref '
                               'laue_starting_right - @ref laue_buffer_right , '
                               'L_z/2 + @ref laue_expand_right ]. This is only '
                               'for Laue-RISM. }',
                       'type': 'REAL'},
 'laue_expand_left': {'default': '-1.0',
                      'description': 'Left expansion for the Laue method',
                      'info': ' If positive value, set the ending position '
                              'offset [in a.u.] of the solvent region on '
                              'left-hand side of the unit cell, measured '
                              'relative to the unit cell edge. (the solvent '
                              'region ends at z = - [L_z/2 + @ref '
                              'laue_expand_left].) This is only for Laue-RISM. '
                              '}',
                      'type': 'REAL'},
 'laue_expand_right': {'default': '-1.0',
                       'description': 'Right expansion for the Laue method',
                       'info': ' If positive value, set the ending position '
                               'offset [in a.u.] of the solvent region on '
                               'right-hand side of the unit cell, measured '
                               'relative to the unit cell edge. (the solvent '
                               'region ends at z = + [L_z/2 + @ref '
                               'laue_expand_right].) This is only for '
                               'Laue-RISM. }',
                       'type': 'REAL'},
 'laue_nfit': {'default': '4',
               'description': 'Number of fitting coefficients for the Laue '
                              'method',
               'info': ' The number of z-grid points for the polynomial fit '
                       'along the cell edge. This is only for Laue-RISM. }',
               'type': 'INTEGER'},
 'laue_starting_left': {'default': '0.0',
                        'description': 'Starting left value for the Laue '
                                       'method',
                        'info': ' Set the starting position [in a.u.] of the '
                                'solvent region on left-hand side of the unit '
                                'cell. Then the solvent region is defined as [ '
                                '-L_z/2 - @ref laue_expand_left , @ref '
                                'laue_starting_left ], where distribution '
                                'functions are finite. This is only for '
                                'Laue-RISM. }',
                        'type': 'REAL'},
 'laue_starting_right': {'default': '0.0',
                         'description': 'Starting right value for the Laue '
                                        'method',
                         'info': ' Set the starting position [in a.u.] of the '
                                 'solvent region on right-hand side of the '
                                 'unit cell. Then the solvent region is '
                                 'defined as [ @ref laue_starting_right , '
                                 'L_z/2 + @ref laue_expand_right ], where '
                                 'distribution functions are finite. This is '
                                 'only for Laue-RISM. }',
                         'type': 'REAL'},
 'laue_wall': {'default': 'auto',
               'description': 'Type of wall for the Laue method',
               'info': ' Set the repulsive wall with (1/r)^12 term of '
                       'Lennard-Jones potential. This is only for Laue-RISM. }',
               'input_type': 'select_multiple',
               'options': ['none', 'auto', 'manual'],
               'type': 'CHARACTER'},
 'laue_wall_epsilon': {'default': '0.1',
                       'description': 'Epsilon parameter for the Laue wall',
                       'info': ' The Lennard-Jones potential of the repulsive '
                               'wall. Here, you can set the parameter '
                               "'epsilon' (kcal/mol). This is only for "
                               "Laue-RISM and @ref laue_wall /= 'none' . }",
                       'type': 'REAL'},
 'laue_wall_lj6': {'default': '.FALSE.',
                   'description': 'Lennard-Jones 6-12 parameter for the Laue '
                                  'wall',
                   'info': ' If .TRUE., the attractive term -(1/r)^6 of '
                           'Lennard-Jones potential is added. This is only for '
                           "Laue-RISM and @ref laue_wall /= 'none' . }",
                   'input_type': 'select_multiple',
                   'options': ['.FALSE.', '.TRUE.'],
                   'type': 'LOGICAL'},
 'laue_wall_rho': {'default': '0.01',
                   'description': 'Density parameter for the Laue wall',
                   'info': ' The density (1/bohr^3) of the repulsive wall. '
                           'This is only for Laue-RISM and @ref laue_wall /= '
                           "'none' . }",
                   'type': 'REAL'},
 'laue_wall_sigma': {'default': '4.0',
                     'description': 'Sigma parameter for the Laue wall',
                     'info': ' The Lennard-Jones potential of the repulsive '
                             "wall. Here, you can set the parameter 'sigma' "
                             '(Angstrom). This is only for Laue-RISM and @ref '
                             "laue_wall /= 'none' . }",
                     'type': 'REAL'},
 'laue_wall_z': {'default': '0.0',
                 'description': 'Z parameter for the Laue wall',
                 'info': ' Set the edge position [in a.u.] of the repulsive '
                         'wall. If @ref laue_expand_right > 0.0, the repulsive '
                         'wall is defined on [ -inf , @ref laue_wall_z ]. If '
                         '@ref laue_expand_left > 0.0, the repulsive wall is '
                         'defined on [ @ref laue_wall_z , inf ]. This is only '
                         "for Laue-RISM and @ref laue_wall == 'manual' . }",
                 'type': 'REAL'},
 'mdiis1d_size': {'default': '20',
                  'description': 'Size of the 1D-MDIIS buffer',
                  'info': ' Size of Modified DIIS (MDIIS) for 1D-RISM. }',
                  'type': 'INTEGER'},
 'mdiis1d_step': {'default': '0.5D0',
                  'description': 'Step of the 1D-MDIIS buffer',
                  'info': ' Step of Modified DIIS (MDIIS) for 1D-RISM. }',
                  'type': 'REAL'},
 'mdiis3d_size': {'default': '10',
                  'description': 'Size of the 3D-MDIIS buffer',
                  'info': ' Size of Modified DIIS (MDIIS) for 3D-RISM. }',
                  'type': 'INTEGER'},
 'mdiis3d_step': {'default': '0.8D0',
                  'description': 'Step of the 3D-MDIIS buffer',
                  'info': ' Step of Modified DIIS (MDIIS) for 3D-RISM. }',
                  'type': 'REAL'},
 'nsolv': {'description': 'Number of solvent atoms',
           'info': ' The number of solvents (i.e. molecular species) in the '
                   'unit cell }',
           'status': '',
           'type': 'INTEGER'},
 'rism1d_bond_width': {'description': 'Bond width for 1D-RISM calculation',
                       'info': ' Gaussian width of bonds to smear '
                               'intra-molecular correlation for 1D-RISM. If '
                               '3D-RISM calculation, default is 0. If '
                               'Laue-RISM calculation, default is 2 / '
                               'SQRT(@ref ecutwfc). }',
                       'type': 'REAL'},
 'rism1d_conv_thr': {'default': '1.D-8',
                     'description': 'Convergence threshold for 1D-RISM '
                                    'calculation',
                     'info': ' Convergence threshold for 1D-RISM. }',
                     'type': 'REAL'},
 'rism1d_dielectric': {'default': '-1.0D0',
                       'description': 'Dielectric constant for 1D-RISM '
                                      'calculation',
                       'info': ' Dielectric constant for 1D-RISM. If @ref '
                               'rism1d_dielectric > 0, dielectrically '
                               'consistent RISM (DRISM) is performed.  For '
                               'details of DRISM, see: J.S.Perkyns and '
                               'B.M.Pettitt, CPL 1992, 190, 626, '
                               'doi:10.1016/0009-2614(92)85201-K }',
                       'type': 'REAL'},
 'rism1d_maxstep': {'default': '50000',
                    'description': 'Maximum steps for 1D-RISM calculation',
                    'info': ' Maximum number of iterations in a 1D-RISM step. '
                            '}',
                    'type': 'INTEGER'},
 'rism1d_molesize': {'default': '2.0D0',
                     'description': 'Molecule size for 1D-RISM calculation',
                     'info': ' Size of solvent molecules (a.u.) for 1D-RISM. '
                             'This is used only if @ref rism1d_dielectric > 0. '
                             'If you have large molecules, you have to set ~ '
                             '20 a.u. . }',
                     'type': 'REAL'},
 'rism1d_nproc': {'default': '128',
                  'description': 'Number of processors for 1D-RISM calculation',
                  'info': ' Number of processes to calculate 1D-RISM. }',
                  'type': 'INTEGER'},
 'rism3d_conv_level': {'description': 'Convergence level for 3D-RISM '
                                      'calculation',
                       'info': ' Convergence level of 3D-RISM. }',
                       'type': 'REAL'},
 'rism3d_conv_thr': {'description': 'Convergence threshold for 3D-RISM '
                                    'calculation',
                     'info': ' Convergence threshold for 3D-RISM. }',
                     'type': 'REAL'},
 'rism3d_maxstep': {'default': '5000',
                    'description': 'Maximum steps for 3D-RISM calculation',
                    'info': ' Maximum number of iterations in a 3D-RISM step. '
                            '}',
                    'type': 'INTEGER'},
 'rism3d_planar_average': {'description': 'Planar average for 3D-RISM '
                                          'calculation',
                           'info': ' If .TRUE., planar averages of solvent '
                                   'densities and potentials are calculated '
                                   "and written to 'prefix.rism1'. For "
                                   '3D-RISM, default is .FALSE. For Laue-RISM, '
                                   'default is .TRUE. }',
                           'input_type': 'select_multiple',
                           'options': ['.FALSE.', '.TRUE.'],
                           'type': 'LOGICAL'},
 'smear1d': {'default': '2.D0',
             'description': 'Smearing parameter for 1D-RISM calculation',
             'info': ' Coulomb smearing radius (a.u.) for 1D-RISM. }',
             'type': 'REAL'},
 'smear3d': {'default': '2.D0',
             'description': 'Smearing parameter for 3D-RISM calculation',
             'info': ' Coulomb smearing radius (a.u.) for 3D-RISM. }',
             'type': 'REAL'},
 'starting1d': {'description': 'Starting point for 1D-RISM calculation',
                'input_type': 'select_multiple',
                'options': ['zero', 'file', 'fix'],
                'type': 'CHARACTER'},
 'starting3d': {'description': 'Starting point for 3D-RISM calculation',
                'input_type': 'select_multiple',
                'options': ['zero', 'file'],
                'type': 'CHARACTER'},
 'tempv': {'default': '300.D0',
           'description': 'Temperature for the volume',
           'info': ' Temperature (Kelvin) of solvents. }',
           'type': 'REAL'}}
