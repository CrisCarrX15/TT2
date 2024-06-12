IONS_DICT = { 'ion_positions': {'default': 'default',
                   'description': 'Ion positions',
                   'info': 'Available options are:\n'
                        '\'default\' :\n'
                        'if restarting, use atomic positions read from the\n'
                        'restart file; in all other cases, use atomic\n'
                        'positions from standard input.\n'
                        '\'from_input\' :\n'
                        'read atomic positions from standard input, even if restarting.',
                   'input_type': 'select_multiple',
                   'options': ['default','from_input'],
                   'type': 'CHARACTER'},
'ion_dynamics': {'description': 'Ion dynamics type',
                  'info': ' Specify the type of ionic dynamics.  For different '
                          'type of calculation different possibilities are '
                          'allowed and different default values apply:  @b '
                          "CASE ( @ref calculation == 'relax' ) } @b CASE ( "
                          "@ref calculation == 'md' ) @b CASE ( @ref "
                          "calculation == 'vc-relax' ) @b CASE ( @ref "
                          "calculation == 'vc-md' )",
                  'input_type': 'select_multiple',
                  'options': ['bfgs',
                              'damp',
                              'fire'],
                  'type': 'CHARACTER'},
'ion_temperature' : {
                        'default' : 'not_controlled',
                        'description' : 'Ion temperature',
                        'info' : '\'rescaling\' : control ionic temperature via velocity rescaling'
                                '(first method) see parameters tempw, tolp, and'
                                'nraise (for VC-MD only). This rescaling method'
                                'is the only one currently implemented in VC-MD'
                                '\'rescale-v\' : control ionic temperature via velocity rescaling'
                                '(second method) see parameters tempw and nraise'
                                '\'rescale-T\' : scale temperature of the thermostat every nraise steps'
                                'by delta_t, starting from tempw.'
                                'The temperature is controlled via velocitiy rescaling.'
                                '\'reduce-T\' : reduce temperature of the thermostat every nraise steps'
                                'by the (negative) value delta_t, starting from tempw.'
                                'If  delta_t is positive, the target temperature is augmented.'
                                'The temperature is controlled via velocitiy rescaling.'
                                '\'berendsen\' : control ionic temperature using "soft" velocity'
                                'rescaling - see parameters tempw and nraise'
                                '\'andersen\' : control ionic temperature using Andersen thermostat'
                                'see parameters tempw and nraise'
                                '\'svr\' : control ionic temperature using stochastic-velocity rescaling'
                                '(Donadio, Bussi, Parrinello, J. Chem. Phys. 126, 014101, 2007),'
                                'with parameters tempw and nraise.'
                                '\'initial\' : initialize ion velocities to temperature tempw'
                                'and leave uncontrolled further on'
                                '\'not_controlled\' : (default) ionic temperature is not controlled',
                        'input_type' : 'select_multiple',
                        'options' : ['rescaling', 'rescale-v', 'rescale-T', 'reduce-T', 'berendsen',
                                     'andersen', 'svr', 'initial', 'not_controlled'],
                        'type' : 'CHARACTER'},
 'ion_velocities': {'default': 'default',
                    'description': 'Ion velocities',
                    'info' : 'Initial ionic velocities. Available options are:'
                                '\'default\' :\n'
                                'start a new simulation from random thermalized\n'
                                'distribution of velocities if tempw is set,\n'
                                'with zero velocities otherwise; restart from\n'
                                'atomic velocities read from the restart file\n'
                                '\'from_input\' :\n'
                                'start or continue the simulation with atomic'
                                'velocities read from standard input',
                    'input_type': 'select_multiple',
                    'options': ['default', 'from_input'],
                    'type': 'CHARACTER'},
 'pot_extrapolation': {'default': 'atomic',
                       'description': 'Potential extrapolation',
                       'info': 'Used to extrapolate the potential from '
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
                       'type': 'CHARACTER'},
 'tempw': {'default': '300.D0',
                    'description': 'Temperature for wavefunctions',
                    'info' : 'Starting temperature (Kelvin) in MD runs'
                        'target temperature for most thermostats.',
                    'type': 'REAL'},
 'tolp': {'default': '100.D0',
                    'description': 'Tolerance for forces',
                    'info' : 'Tolerance for velocity rescaling. Velocities are rescaled if'
                                'the run-averaged and target temperature differ more than tolp.',
                    'type': 'REAL'},
 'delta_t': {'default': '1.D0',
                    'description': 'Time step',
                    'info' : 'For delta_t < 0, the actual average rate of heating or cooling'
                        'should be roughly C*delta_t/(nraise*dt) (C=1 for an'
                        'ideal gas, C=0.5 for a harmonic solid, theorem of energy'
                        'equipartition between all quadratic degrees of freedom).',
                    'type': 'REAL'},
 'nraise': {'default': '1',
                    'description': 'Number of steps to raise delta_t',
                    'info' : 'if ion_temperature == \'reduce-T\' : every nraise steps the instantaneous temperature is'
                        'reduced by -delta_t (i.e. delta_t is added to the temperature)'
                        'if ion_temperature == \'rescale-v\' : every nraise steps the average temperature, computed from'
                        'the last nraise steps, is rescaled to tempw'
                        'if ion_temperature == \'rescaling\' and calculation == \'vc-md\' : every nraise steps the instantaneous temperature'
                        'is rescaled to tempw'
                        'if ion_temperature == \'berendsen\' : the \"rise time\" parameter is given in units of the time step:'
                        'tau = nraise*dt, so dt/tau = 1/nraise'
                        'if ion_temperature == \'andersen\' : the "collision frequency" parameter is given as nu=1/tau'
                        'defined above, so nu*dt = 1/nraise'
                        'if ion_temperature == \'svr\' : the \"characteristic time\" of the thermostat is set to'
                        'tau = nraise*dt',
                    'type': 'INTEGER'},
 'refold_pos': {'default': '.FALSE.',
                      'description': 'Fold positions into the first WS cell',
                      'info': 'This keyword applies only in the case of molecular dynamics or'
                                'damped dynamics. If true the ions are refolded at each step into'
                                'the supercell.',
                      'input_type': 'select_multiple',
                      'options': ['.FALSE.', '.TRUE.'],
                      'type': 'LOGICAL'},
 'upscale': {'default': '100.D0',
                    'description': 'Upscaling factor for CG',
                    'info' : 'Max reduction factor for conv_thr during structural optimization'
                                'conv_thr is automatically reduced when the relaxation'
                                'approaches convergence so that forces are still accurate,'
                                'but conv_thr will not be reduced to less that conv_thr / upscale.',
                    'type': 'REAL'},
 'bfgs_ndim': {'default': '1',
                    'description': 'Number of stored gradients in BFGS',
                    'info' : 'Number of old forces and displacements vectors used in the'
                                'PULAY mixing of the residual vectors obtained on the basis'
                                'of the inverse hessian matrix given by the BFGS algorithm.'
                                'When bfgs_ndim = 1, the standard quasi-Newton BFGS method is'
                                'used.(bfgs only)',
                    'type': 'INTEGER'},
 'trust_radius_max': {'default': '0.8D0',
                    'description': 'Maximum trust radius',
                    'info' : 'Maximum ionic displacement in the structural relaxation.'
                                '(bfgs only)',
                    'type': 'REAL'},
 'trust_radius_min': {'default': '1.D-3',
                    'description': 'Minimum trust radius',
                    'info' : 'Minimum ionic displacement in the structural relaxation'
                                'BFGS is reset when trust_radius < trust_radius_min.'
                                '(bfgs only)',
                    'type': 'REAL'},
 'trust_radius_ini': {'default': '0.5D0',
                    'description': 'Initial trust radius',
                    'info' : 'Initial ionic displacement in the structural relaxation.'
                                '(bfgs only)',
                    'type': 'REAL'},
 'w_1': {'default': '0.01D0',
                    'description': 'First weight for forces in optimization algorithms',
                    'info' : 'First weight for forces in optimization algorithms',
                    'type': 'REAL'},
 'w_2': {'default': '0.5D0',
                    'description': 'Second weight for forces in optimization algorithms',
                    'info' : 'Parameters used in line search based on the Wolfe conditions.'
                                '(bfgs only)',
                    'type': 'REAL'},
 'fire_alpha_init': {'default': '0.2D0',
                    'description': 'Initial value for alpha in FIRE algorithm',
                    'info' : 'Initial value of the alpha mixing factor in the FIRE minimization scheme;'
                                'recommended values are between 0.1 and 0.3',
                    'type': 'REAL'},
 'fire_falpha': {'default': '0.99D0',
                    'description': 'Decay factor for alpha in FIRE algorithm',
                    'info' : 'Scaling of the alpha mixing parameter for steps with P > 0;',
                    'type': 'REAL'},
 'fire_nmin': {'default': '5',
                    'description': 'Minimum number of steps in FIRE algorithm',
                    'info' : 'Minimum number of steps with P > 0 before increase of dt',
                    'type': 'INTEGER'},
 'fire_f_inc': {'default': '1.1D0',
                    'description': 'Factor for increasing delta_t in FIRE algorithm',
                    'info' : 'Factor for increasing dt',
                    'type': 'REAL'},
 'fire_f_dec': {'default': '0.5D0',
                    'description': 'Factor for decreasing delta_t in FIRE algorithm',
                    'info' : 'Factor for decreasing dt',
                    'type': 'REAL'},
 'fire_dtmax': {'default': '10.D0',
                    'description': 'Maximum time step in FIRE algorithm',
                    'info' : 'Determines the maximum value of dt in the FIRE minimization;'
                                'dtmax = fire_dtmax*dt',
                    'type': 'REAL'}
}
