CONTROL_DICT = {'calculation': {'default': 'scf',
                 'info': ' A string describing the task to be performed. '
                         ' (vc = variable-cell).',
                 'input_type': 'select_multiple',
                 'options': ['scf',
                             'nscf',
                             'bands',
                             'relax',
                             #'md',
                             #'vc-relax',
                             #'vc-md'
                             ],
                 'type': 'CHARACTER',
                 'status': 'REQUIRED',
                 'description' : 'Calculation type'},
 'prefix': {'default': 'pwscf',
            'info': ' prepended to input/output filenames: prefix.wfc, '
                    'prefix.rho, etc. ',
            'type': 'CHARACTER',
            #'status' : 'REQUIRED',
            'description' : 'File prefix.'},
 'pseudo_dir': {'info': ' directory containing pseudopotential files ',
                'type': 'CHARACTER',
                #'status' : 'REQUIRED',
                'description' : 'Pseudopotential directory.'},
 'outdir': {'info': ' input, temporary, output files are found in this '
                    'directory, see also @ref wfcdir ',
            'type': 'CHARACTER',
            #'status' : 'REQUIRED',
            'description' : 'Output directory.'},
 'dipfield': {'default': '.FALSE.',
              'info': ' If .TRUE. and @ref tefield==.TRUE. a dipole correction '
                      'is also added to the bare ionic potential - implements '
                      'the recipe of L. Bengtsson, PRB 59, 12301 (1999). See '
                      'variables @ref edir, @ref emaxpos, @ref eopreg for the '
                      'form of the correction. Must be used ONLY in a slab '
                      'geometry, for surface calculations, with the '
                      'discontinuity FALLING IN THE EMPTY SPACE.',
              'input_type': 'select_multiple',
              'options': ['.FALSE.', '.TRUE.'],
              'type': 'LOGICAL',
              'description' : 'Electric dipole field.'},
 'disk_io': {'default': 'see below',
             'info': ' Specifies the amount of disk I/O activity: (only for '
                     'binary files and xml data file in data directory; other '
                     'files printed at each molecular dynamics / structural '
                     'optimization step are not controlled by this option ) '
                     "@b Default is @b 'low' for the scf case, @b 'medium' "
                     'otherwise. Note that the needed RAM increases as disk '
                     'I/O decreases',
             'input_type': 'select_multiple',
             'options': ['high', 'medium', 'low', 'nowf', 'minimal', 'none'],
             'type': 'CHARACTER',
             'description' : 'Disk input/output.'},
 'dt': {'default': '20.D0',
        'info': ' time step for molecular dynamics, in Rydberg atomic units (1 '
                'a.u.=4.8378 * 10^-17 s : beware, the CP code uses Hartree '
                'atomic units, half that much!!!)',
        'type': 'REAL',
        'description' : 'Time step.'},
 'etot_conv_thr': {'default': '1.0D-4',
                   'info': ' Convergence threshold on total energy (a.u) for '
                           'ionic minimization: the convergence criterion is '
                           'satisfied when the total energy changes less than '
                           '@ref etot_conv_thr between two consecutive scf '
                           'steps. Note that @ref etot_conv_thr is extensive, '
                           'like the total energy. See also @ref forc_conv_thr '
                           '- both criteria must be satisfied',
                   'type': 'REAL',
                   'description' : 'Total energy convergence threshold.'},
 'forc_conv_thr': {'default': '1.0D-3',
                   'info': ' Convergence threshold on forces (a.u) for ionic '
                           'minimization: the convergence criterion is '
                           'satisfied when all components of all forces are '
                           'smaller than @ref forc_conv_thr. See also @ref '
                           'etot_conv_thr - both criteria must be satisfied',
                   'type': 'REAL',
                   'description' : 'Force convergence threshold.'},
 'gate': {'default': '.FALSE.',
          'info': ' In the case of charged cells (@ref tot_charge .ne. 0) '
                  'setting gate = .TRUE. represents the counter charge (i.e. '
                  '-tot_charge) not by a homogeneous background charge but '
                  'with a charged plate, which is placed at @ref zgate (see '
                  'below). Details of the gate potential can be found in T. '
                  'Brumme, M. Calandra, F. Mauri; PRB 89, 245406 (2014). Note, '
                  'that in systems which are not symmetric with respect to the '
                  'plate, one needs to enable the dipole correction! (@ref '
                  'dipfield=.true.). Currently, symmetry can be used with '
                  'gate=.true. but carefully check that no symmetry is '
                  'included which maps @i z to -@i z even if in principle one '
                  'could still use them for symmetric systems (i.e. no dipole '
                  'correction). For @ref nosym=.false. verbosity is set to '
                  '\'high\'. Note: this option was called "monopole" in v6.0 '
                  'and 6.1 of pw.x',
          'input_type': 'select_multiple',
          'options': ['.FALSE.', '.TRUE.'],
          'type': 'LOGICAL',
          'description' : 'Gate.'},
 'gdir': {'info': ' For Berry phase calculation: direction of the k-point '
                  'strings in reciprocal space. Allowed values: 1, 2, 3 '
                  '1=first, 2=second, 3=third reciprocal lattice vector For '
                  'calculations with finite electric fields (@ref '
                  'lelfield==.true.) "gdir" is the direction of the field.',
          'input_type': 'select_multiple',
          'options' : ['1','2','3'],
          'type': 'INTEGER',
          'description' : 'Propagation direction.'},
 'iprint': {'default': 'write only at convergence',
            'info': " When @ref calculation == 'md' (molecular dynamics) "
                    'trajectory is written every @i iprint md steps.',
            'type': 'INTEGER',
            'description' : 'Print level.'},
 'lberry': {'default': '.FALSE.',
            'info': ' If .TRUE. perform a Berry phase calculation. See the '
                    'header of PW/src/bp_c_phase.f90 for documentation.',
            'input_type': 'select_multiple',
            'options': ['.FALSE.', '.TRUE.'],
            'type': 'LOGICAL',
            'description' : 'Berry phase calculation.'},
 'lelfield': {'default': '.FALSE.',
              'info': ' If .TRUE. a homogeneous finite electric field '
                      'described through the modern theory of the polarization '
                      'is applied. This is different from @ref tefield == '
                      '.true. ! ',
              'input_type': 'select_multiple',
              'options': ['.FALSE.', '.TRUE.'],
              'type': 'LOGICAL',
              'description' : 'Electric field calculation.'},
 'lfcp': {'default': '.FALSE.',
          'info': ' If .TRUE. perform a constant bias potential (constant-mu) '
                  'calculation for a system with ESM method. See the header of '
                  'PW/src/fcp_module.f90 for documentation. To perform the '
                  'calculation, you must set a namelist FCP.  NB: - The total '
                  'energy displayed in output includes the potentiostat '
                  "contribution (-mu*N). - @ref calculation must be 'relax' or "
                  "'md'. - @ref assume_isolated = 'esm' and @ref esm_bc = "
                  "'bc2' or 'bc3' must be set in @ref SYSTEM namelist. - "
                  "ESM-RISM is also supported (@ref assume_isolated = 'esm' "
                  "and @ref esm_bc = 'bc1' and @ref trism = .TRUE.). - @ref "
                  'ignore_wolfe is always .TRUE., for BFGS. ',
          'input_type': 'select_multiple',
          'options': ['.FALSE.', '.TRUE.'],
          'type': 'LOGICAL',
          'description' : 'Linear response calculation.'},
 'lkpoint_dir': {'info': ' OBSOLETE - NO LONGER IMPLEMENTED ',
                 'input_type': 'select_multiple',
                 'options': ['.FALSE.', '.TRUE.'],
                 'type': 'LOGICAL',
                 'description' : 'K-point direction.'},
 'lorbm': {'default': '.FALSE.',
           'info': ' If @b .TRUE. perform orbital magnetization calculation. '
                   'If finite electric field is applied (@ref '
                   'lelfield==.true.) only Kubo terms are computed [for '
                   'details see New J. Phys. 12, 053032 (2010), '
                   'doi:10.1088/1367-2630/12/5/053032].  The type of '
                   "calculation is @b 'nscf' and should be performed on an "
                   'automatically generated uniform grid of k points.  Works '
                   'ONLY with norm-conserving pseudopotentials. ',
           'input_type': 'select_multiple',
           'options': ['.FALSE.', '.TRUE.'],
           'type': 'LOGICAL',
           'description' : 'Orbital magnetic moment calculation.'},
 'max_seconds': {'default': '1.D+7, or 150 days, i.e. no time limit',
                 'info': ' Jobs stops after @ref max_seconds CPU time. Use '
                         'this option in conjunction with option @ref '
                         'restart_mode if you need to split a job too long to '
                         'complete into shorter jobs that fit into your batch '
                         'queues. ',
                 'type': 'REAL',
                 'description' : 'Maximum seconds.'},
 'nberrycyc': {'default': '1',
               'info': ' In the case of a finite electric field  ( @ref '
                       'lelfield == .TRUE. ) it defines the number of '
                       'iterations for converging the wavefunctions in the '
                       'electric field Hamiltonian, for each external '
                       'iteration on the charge density ',
               'type': 'INTEGER',
               'description' : 'Berry phase cycles.'},
 'nppstr': {'info': ' For Berry phase calculation: number of k-points to be '
                    'calculated along each symmetry-reduced string. The same '
                    'for calculation with finite electric fields (@ref '
                    'lelfield==.true.). ',
            'type': 'INTEGER',
            'description' : 'Post-processing steps.'},
 'nstep': {'info': ' number of molecular-dynamics or structural optimization '
                   'steps performed in this run. If set to 0, the code '
                   'performs a quick "dry run", stopping just after '
                   'initialization. This is useful to check for input '
                   'correctness and to have the summary printed. NOTE: in MD '
                   'calculations, the code will perform "nstep" steps even if '
                   'restarting from a previously interrupted calculation. ',
           'type': 'INTEGER',
           'description' : 'Time steps.'},
 'restart_mode': {'default': 'from_scratch',
                  'info': '',
                  'input_type': 'select_multiple',
                  'options': ['from_scratch', 'restart'],
                  'type': 'CHARACTER',
                  'description' : 'Restart mode.'},
 'tefield': {'default': '.FALSE.',
             'info': ' If .TRUE. a saw-like potential simulating an electric '
                     'field is added to the bare ionic potential. See '
                     'variables @ref edir, @ref eamp, @ref emaxpos, @ref '
                     'eopreg for the form and size of the added potential.  ',
             'input_type': 'select_multiple',
             'options': ['.FALSE.', '.TRUE.'],
             'type': 'LOGICAL',
             'description' : 'Electric field.'},
 'title': {'info': ' reprinted on output. ', 'type': 'CHARACTER', 'description' : 'Simulation title.'},
 'tprnfor': {'info': ' calculate forces. It is set to .TRUE. automatically if '
                     "@ref calculation == 'relax','md','vc-md' }",
             'input_type': 'select_multiple',
             'options': ['.FALSE.', '.TRUE.'],
             'type': 'LOGICAL',
             'description' : 'Forces print.'},
 'trism': {'default': '.FALSE.',
           'info': ' If .TRUE. perform a 3D-RISM-SCF calculation [for details '
                   'see H.Sato et al., JCP 112, 9463 (2000), '
                   "doi:10.1063/1.481564]. The solvent's distributions are "
                   'calculated by 3D-RISM, though solute is treated as SCF. '
                   'The charge density and the atomic positions are optimized, '
                   'simultaneously with the solvents. To perform the '
                   'calculation, you must set a namelist @ref RISM and a card '
                   "@ref SOLVENTS.  If @ref assume_isolated = 'esm' and @ref "
                   "esm_bc = 'bc1', Laue-RISM is calculated instead of 3D-RISM "
                   'and coupled with ESM method (i.e. ESM-RISM). [for details '
                   'see S.Nishihara and M.Otani, PRB 96, 115429 (2017)].  The '
                   'default of @ref mixing_beta is 0.2 for both 3D-RISM and '
                   'Laue-RISM.  For structural relaxation with BFGS, @ref '
                   'ignore_wolfe is always .TRUE. . ',
           'input_type': 'select_multiple',
           'options': ['.FALSE.', '.TRUE.'],
           'type': 'LOGICAL',
           'description' : 'Tracing image.'},
 'tstress': {'default': '.false.',
             'info': ' calculate stress. It is set to .TRUE. automatically if '
                     "@ref calculation == 'vc-md' or 'vc-relax' }",
             'input_type': 'select_multiple',
             'options': ['.FALSE.', '.TRUE.'],
             'type': 'LOGICAL',
             'description' : 'Stress tensor.'},
 'twochem': {'default': '.FALSE.',
             'info': ' IF .TRUE. , a two chemical potential calculation for '
                     'the simulation of photoexcited systems is performed, '
                     'constraining a fraction of the electrons in the '
                     'conduction manifold. See G. Marini, M. Calandra; PRB '
                     '104, 144103 (2021). Note: requires @ref occupations to '
                     "be set to 'smearing'. }",
             'input_type': 'select_multiple',
             'options': ['.FALSE.', '.TRUE.'],
             'type': 'LOGICAL',
             'description' : 'Two-chemical-potential.'},
 'verbosity': {'default': 'low',
               'info': ' Currently two verbosity levels are implemented: } @b '
                       "'debug' and @b 'medium' have the same effect as @b "
                       "'high'; @b 'default' and @b 'minimal' as @b 'low' }",
               'input_type': 'select_multiple',
               'options': ['high', 'low'],
               'type': 'CHARACTER',
               'description' : 'Output verbosity.'},
 'wf_collect': {'info': '',
                'input_type': 'select_multiple',
                'options': ['.FALSE.', '.TRUE.'],
                'type': 'LOGICAL',
                'description' : 'Wavefunction collection.'},
 'wfcdir': {'default': 'same as @ref outdir',
            'info': ' This directory specifies where to store files generated '
                    'by each processor (*.wfc{N}, *.igk{N}, etc.). Useful for',
            'type': 'CHARACTER',
            'description' : 'Wavefunction directory.'}}
