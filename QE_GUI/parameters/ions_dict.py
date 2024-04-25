IONS_DICT = {'ion_dynamics': {'description': 'Ion dynamics type',
                  'info': ' Specify the type of ionic dynamics.  For different '
                          'type of calculation different possibilities are '
                          'allowed and different default values apply:  @b '
                          "CASE ( @ref calculation == 'relax' ) } @b CASE ( "
                          "@ref calculation == 'md' ) } @b CASE ( @ref "
                          "calculation == 'vc-relax' ) } @b CASE ( @ref "
                          "calculation == 'vc-md' ) }",
                  'input_type': 'select_multiple',
                  'options': ['bfgs',
                              'damp',
                              'fire',
                              'verlet',
                              'langevin',
                              'langevin-smc',
                              'bfgs',
                              'damp',
                              'beeman'],
                  'type': 'CHARACTER'},
 'ion_positions': {'default': 'default',
                   'description': 'Ion positions',
                   'info': '',
                   'input_type': 'select_multiple',
                   'options': ['from_input'],
                   'type': 'CHARACTER'},
 'ion_velocities': {'default': 'default',
                    'description': 'Ion velocities',
                    'input_type': 'select_multiple',
                    'options': ['from_input'],
                    'type': 'CHARACTER'},
 'pot_extrapolation': {'default': 'atomic',
                       'description': 'Potential extrapolation',
                       'info': ' Used to extrapolate the potential from '
                               "preceding ionic steps. } Note: 'first_order' "
                               "and 'second-order' extrapolation make sense "
                               'only for molecular dynamics calculations ',
                       'input_type': 'select_multiple',
                       'options': ['none',
                                   'atomic',
                                   'first_order',
                                   'second_order'],
                       'type': 'CHARACTER'},
 'remove_rigid_rot': {'default': '.FALSE.',
                      'description': 'Remove rigid rotation',
                      'info': ' This keyword is useful when simulating the '
                              'dynamics and/or the thermodynamics of an '
                              'isolated system. If set to true the total '
                              'torque of the internal forces is set to zero by '
                              'adding new forces that compensate the spurious '
                              'interaction with the periodic images. This '
                              'allows for the use of smaller supercells.  '
                              'BEWARE: since the potential energy is no longer '
                              'consistent with the forces (it still contains '
                              'the spurious interaction with the repeated '
                              'images), the total energy is not conserved '
                              'anymore. However the dynamical and '
                              'thermodynamical properties should be in closer '
                              'agreement with those of an isolated system. '
                              'Also the final energy of a structural '
                              'relaxation will be higher, but the relaxation '
                              'itself should be faster. ',
                      'input_type': 'select_multiple',
                      'options': ['.FALSE.', '.TRUE.'],
                      'type': 'LOGICAL'},
 'wfc_extrapolation': {'default': 'none',
                       'description': 'Wavefunction extrapolation',
                       'info': ' Used to extrapolate the wavefunctions from '
                               'preceding ionic steps. } Note: @b '
                               "'first_order' and @b 'second-order' "
                               'extrapolation make sense only for molecular '
                               'dynamics calculations ',
                       'input_type': 'select_multiple',
                       'options': ['none', 'first_order', 'second_order'],
                       'type': 'CHARACTER'}}
