SYSTEM_DICT = {'ibrav': {'description': 'Bravais lattice type',
           'info': ' Bravais-lattice index. Optional only if space_group is '
                   'set. If ibrav /= 0, specify EITHER [ @ref celldm(1)-@ref '
                   'celldm(6) ] OR [ @ref A, @ref B, @ref C, @ref cosAB, @ref '
                   'cosAC, @ref cosBC ] but NOT both. The lattice parameter '
                   '"alat" is set to alat = celldm(1) (in a.u.) or alat = A '
                   '(in Angstrom); see below for the other parameters. For '
                   'ibrav=0 specify the lattice vectors in @ref '
                   'CELL_PARAMETERS, optionally the lattice parameter alat = '
                   'celldm(1) (in a.u.) or = A (in Angstrom). If not '
                   'specified, the lattice parameter is taken from @ref '
                   'CELL_PARAMETERS IMPORTANT NOTICE 1: with ibrav=0 lattice '
                   'vectors must be given with a sufficiently large number of '
                   'digits and with the correct symmetry, or else symmetry '
                   'detection may fail and strange problems may arise in '
                   'symmetrization. IMPORTANT NOTICE 2: do not use celldm(1) '
                   'or A as a.u. to Ang conversion factor, use the true '
                   'lattice parameters or nothing, specify units in @ref '
                   'CELL_PARAMETERS and @ref ATOMIC_POSITIONS  ibrav      '
                   'structure                   celldm(2)-celldm(6) or: '
                   'b,c,cosbc,cosac,cosab 0          free crystal axis '
                   'provided in input: see card @ref CELL_PARAMETERS  '
                   '1          cubic P (sc) v1 = a(1,0,0),  v2 = a(0,1,0),  v3 '
                   '= a(0,0,1)  2          cubic F (fcc) v1 = (a/2)(-1,0,1),  '
                   'v2 = (a/2)(0,1,1), v3 = (a/2)(-1,1,0)  3          cubic I '
                   '(bcc) v1 = (a/2)(1,1,1),  v2 = (a/2)(-1,1,1),  v3 = '
                   '(a/2)(-1,-1,1) -3          cubic I (bcc), more symmetric '
                   'axis: v1 = (a/2)(-1,1,1), v2 = (a/2)(1,-1,1),  v3 = '
                   '(a/2)(1,1,-1)  4          Hexagonal and Trigonal P        '
                   'celldm(3)=c/a v1 = a(1,0,0),  v2 = a(-1/2,sqrt(3)/2,0),  '
                   'v3 = a(0,0,c/a)  5          Trigonal R, 3fold axis '
                   'c        celldm(4)=cos(gamma) The crystallographic vectors '
                   'form a three-fold star around the z-axis, the primitive '
                   'cell is a simple rhombohedron: v1 = a(tx,-ty,tz),   v2 = '
                   'a(0,2ty,tz),   v3 = a(-tx,-ty,tz) where c=cos(gamma) is '
                   'the cosine of the angle gamma between any pair of '
                   'crystallographic vectors, tx, ty, tz are: '
                   'tx=sqrt((1-c)/2), ty=sqrt((1-c)/6), tz=sqrt((1+2c)/3) '
                   '-5          Trigonal R, 3fold axis <111>    '
                   'celldm(4)=cos(gamma) The crystallographic vectors form a '
                   "three-fold star around <111>. Defining a' = a/sqrt(3) : v1 "
                   "= a' (u,v,v),   v2 = a' (v,u,v),   v3 = a' (v,v,u) where u "
                   'and v are defined as u = tz - 2*sqrt(2)*ty,  v = tz + '
                   'sqrt(2)*ty and tx, ty, tz as for case ibrav=5 Note: if you '
                   'prefer x,y,z as axis in the cubic limit, set  u = tz + '
                   '2*sqrt(2)*ty,  v = tz - sqrt(2)*ty See also the note in '
                   'Modules/latgen.f90  6          Tetragonal P '
                   '(st)               celldm(3)=c/a v1 = a(1,0,0),  v2 = '
                   'a(0,1,0),  v3 = a(0,0,c/a)  7          Tetragonal I '
                   '(bct)              celldm(3)=c/a v1=(a/2)(1,-1,c/a),  '
                   'v2=(a/2)(1,1,c/a),  v3=(a/2)(-1,-1,c/a)  8          '
                   'Orthorhombic P                  celldm(2)=b/a '
                   'celldm(3)=c/a v1 = (a,0,0),  v2 = (0,b,0), v3 = (0,0,c)  '
                   '9          Orthorhombic base-centered(bco) celldm(2)=b/a '
                   'celldm(3)=c/a v1 = (a/2, b/2,0),  v2 = (-a/2,b/2,0),  v3 = '
                   '(0,0,c) -9          as 9, alternate description v1 = '
                   '(a/2,-b/2,0),  v2 = (a/2, b/2,0),  v3 = (0,0,c) '
                   '91          Orthorhombic one-face base-centered A-type '
                   'celldm(2)=b/a celldm(3)=c/a v1 = (a, 0, 0),  v2 = '
                   '(0,b/2,-c/2),  v3 = (0,b/2,c/2)  10          Orthorhombic '
                   'face-centered      celldm(2)=b/a celldm(3)=c/a v1 = '
                   '(a/2,0,c/2),  v2 = (a/2,b/2,0),  v3 = (0,b/2,c/2)  '
                   '11          Orthorhombic body-centered      celldm(2)=b/a '
                   'celldm(3)=c/a v1=(a/2,b/2,c/2),  v2=(-a/2,b/2,c/2),  '
                   'v3=(-a/2,-b/2,c/2)  12          Monoclinic P, unique axis '
                   'c     celldm(2)=b/a celldm(3)=c/a, celldm(4)=cos(ab) '
                   'v1=(a,0,0), v2=(b*cos(gamma),b*sin(gamma),0),  v3 = '
                   '(0,0,c) where gamma is the angle between axis a and b. '
                   '-12          Monoclinic P, unique axis b     celldm(2)=b/a '
                   'celldm(3)=c/a, celldm(5)=cos(ac) v1 = (a,0,0), v2 = '
                   '(0,b,0), v3 = (c*cos(beta),0,c*sin(beta)) where beta is '
                   'the angle between axis a and c  13          Monoclinic '
                   'base-centered        celldm(2)=b/a (unique axis '
                   'c)                 celldm(3)=c/a, celldm(4)=cos(gamma) v1 '
                   '= (  a/2,         0,          -c/2), v2 = (b*cos(gamma), '
                   'b*sin(gamma), 0  ), v3 = (  a/2,         0,           '
                   'c/2), where gamma=angle between axis a and b projected on '
                   'xy plane  -13          Monoclinic base-centered        '
                   'celldm(2)=b/a (unique axis b)                 '
                   'celldm(3)=c/a, celldm(5)=cos(beta) v1 = (  a/2,       '
                   'b/2,             0), v2 = ( -a/2,       b/2,             '
                   '0), v3 = (c*cos(beta),   0,   c*sin(beta)), where '
                   'beta=angle between axis a and c projected on xz plane '
                   'IMPORTANT NOTICE: until QE v.6.4.1, axis for ibrav=-13 had '
                   'a different definition: v1(old) =-v2(now), v2(old) = '
                   'v1(now)  14          Triclinic                       '
                   'celldm(2)= b/a, celldm(3)= c/a, celldm(4)= cos(bc), '
                   'celldm(5)= cos(ac), celldm(6)= cos(ab) v1 = (a, 0, 0), v2 '
                   '= (b*cos(gamma), b*sin(gamma), 0) v3 = (c*cos(beta),  '
                   'c*(cos(alpha)-cos(beta)cos(gamma))/sin(gamma), c*sqrt( 1 + '
                   '2*cos(alpha)cos(beta)cos(gamma) - '
                   'cos(alpha)^2-cos(beta)^2-cos(gamma)^2 )/sin(gamma) ) where '
                   'alpha is the angle between axis b and c beta is the angle '
                   'between axis a and c gamma is the angle between axis a and '
                   'b }',
           'status': '',
           'input_type': 'select_multiple',
           'options': ['0', '1', '2', '3', '-3', '4', '5', '-5', '6', '7', '8', '9', '-9', '91', '10', '11', '12', '-12', '-13', '14'],
           'type': 'INTEGER'},
 'celldm(1)' : {
                'description' : 'Specifies the length scale of the simulation cell in Bohr units or angstroms',
                'info' : 'Crystallographic constants - see the ibrav variable.'
                        'Specify either these OR A,B,C,cosAB,cosBC,cosAC NOT both.'
                        'Only needed values (depending on "ibrav") must be specified'
                        'alat = celldm(1) is the lattice parameter "a" (in BOHR)'
                        'If ibrav==0, only celldm(1) is used if present;'
                        'cell vectors are read from card CELL_PARAMETERS',
                'type': 'REAL'},
 'celldm(2)' : {'description' : 'Specifies the ratio b/a of lattice vectors in orthorhombic cell (celldm(2) = b/a)',
                'info' : 'Crystallographic constants - see the ibrav variable.'
                        'Specify either these OR A,B,C,cosAB,cosBC,cosAC NOT both.'
                        'Only needed values (depending on "ibrav") must be specified'
                        'alat = celldm(1) is the lattice parameter "a" (in BOHR)'
                        'If ibrav==0, only celldm(1) is used if present;'
                        'cell vectors are read from card CELL_PARAMETERS',
                'type': 'REAL'},
 'celldm(3)' : {'description' : 'Specifies the ratio c/a of lattice vectors in orthorhombic cell (celldm(3) = c/a)',
                'info' : 'Crystallographic constants - see the ibrav variable.'
                        'Specify either these OR A,B,C,cosAB,cosBC,cosAC NOT both.'
                        'Only needed values (depending on "ibrav") must be specified'
                        'alat = celldm(1) is the lattice parameter "a" (in BOHR)'
                        'If ibrav==0, only celldm(1) is used if present;'
                        'cell vectors are read from card CELL_PARAMETERS',
                'type': 'REAL'},
 'celldm(4)' : {'description' : 'Specifies the value of cos(α) in monoclinic cell (celldm(4) = cos(α))',
                'info' : 'Crystallographic constants - see the ibrav variable.'
                        'Specify either these OR A,B,C,cosAB,cosBC,cosAC NOT both.'
                        'Only needed values (depending on "ibrav") must be specified'
                        'alat = celldm(1) is the lattice parameter "a" (in BOHR)'
                        'If ibrav==0, only celldm(1) is used if present;'
                        'cell vectors are read from card CELL_PARAMETERS',
                'type': 'REAL'},
 'celldm(5)' : {'description' : 'Specifies the value of cos(β) in monoclinic cell (celldm(5) = cos(β))',
                'info' : 'Crystallographic constants - see the ibrav variable.'
                        'Specify either these OR A,B,C,cosAB,cosBC,cosAC NOT both.'
                        'Only needed values (depending on "ibrav") must be specified'
                        'alat = celldm(1) is the lattice parameter "a" (in BOHR)'
                        'If ibrav==0, only celldm(1) is used if present;'
                        'cell vectors are read from card CELL_PARAMETERS',
                'type': 'REAL'},
 'celldm(6)' : {'description' : 'Specifies the value of cos(γ) in monoclinic cell (celldm(6) = cos(γ))',
                'info' : 'Crystallographic constants - see the ibrav variable.'
                        'Specify either these OR A,B,C,cosAB,cosBC,cosAC NOT both.'
                        'Only needed values (depending on "ibrav") must be specified'
                        'alat = celldm(1) is the lattice parameter "a" (in BOHR)'
                        'If ibrav==0, only celldm(1) is used if present;'
                        'cell vectors are read from card CELL_PARAMETERS',
                'type': 'REAL'},
'A' : {'description' : 'Specifies the length of the first lattice vector',
       'info' : 'Traditional crystallographic constants:'
                'A,B,C in ANGSTROM'
                'cosAB = cosine of the angle between axis a and b (gamma)'
                'cosAC = cosine of the angle between axis a and c (beta)'
                'cosBC = cosine of the angle between axis b and c (alpha)'
                'The axis are chosen according to the value of ibrav.'
                'Specify either these OR celldm but NOT both.'
                'Only needed values (depending on ibrav) must be specified.'
                'The lattice parameter alat = A (in ANGSTROM ).'
                'If ibrav == 0, only A is used if present, and'
                'cell vectors are read from card CELL_PARAMETERS.',
                'type': 'REAL'},
'B' : {'description' : 'Specifies the length of the second lattice vector',
       'info' : 'Traditional crystallographic constants:'
                'A,B,C in ANGSTROM'
                'cosAB = cosine of the angle between axis a and b (gamma)'
                'cosAC = cosine of the angle between axis a and c (beta)'
                'cosBC = cosine of the angle between axis b and c (alpha)'
                'The axis are chosen according to the value of ibrav.'
                'Specify either these OR celldm but NOT both.'
                'Only needed values (depending on ibrav) must be specified.'
                'The lattice parameter alat = A (in ANGSTROM ).'
                'If ibrav == 0, only A is used if present, and'
                'cell vectors are read from card CELL_PARAMETERS.',
                'type': 'REAL'},
'C' : {'description' : 'Specifies the length of the third lattice vector',
       'info' : 'Traditional crystallographic constants:'
                'A,B,C in ANGSTROM'
                'cosAB = cosine of the angle between axis a and b (gamma)'
                'cosAC = cosine of the angle between axis a and c (beta)'
                'cosBC = cosine of the angle between axis b and c (alpha)'
                'The axis are chosen according to the value of ibrav.'
                'Specify either these OR celldm but NOT both.'
                'Only needed values (depending on ibrav) must be specified.'
                'The lattice parameter alat = A (in ANGSTROM ).'
                'If ibrav == 0, only A is used if present, and'
                'cell vectors are read from card CELL_PARAMETERS.',
                'type': 'REAL'},
'cosAB' : {'description' : 'Specifies the cosine of the angle between lattice vectors A and B',
       'info' : 'Traditional crystallographic constants:'
                'A,B,C in ANGSTROM'
                'cosAB = cosine of the angle between axis a and b (gamma)'
                'cosAC = cosine of the angle between axis a and c (beta)'
                'cosBC = cosine of the angle between axis b and c (alpha)'
                'The axis are chosen according to the value of ibrav.'
                'Specify either these OR celldm but NOT both.'
                'Only needed values (depending on ibrav) must be specified.'
                'The lattice parameter alat = A (in ANGSTROM ).'
                'If ibrav == 0, only A is used if present, and'
                'cell vectors are read from card CELL_PARAMETERS.',
                'type': 'REAL'},
'cosAC' : {'description' : 'Specifies the cosine of the angle between lattice vectors A and C',
       'info' : 'Traditional crystallographic constants:'
                'A,B,C in ANGSTROM'
                'cosAB = cosine of the angle between axis a and b (gamma)'
                'cosAC = cosine of the angle between axis a and c (beta)'
                'cosBC = cosine of the angle between axis b and c (alpha)'
                'The axis are chosen according to the value of ibrav.'
                'Specify either these OR celldm but NOT both.'
                'Only needed values (depending on ibrav) must be specified.'
                'The lattice parameter alat = A (in ANGSTROM ).'
                'If ibrav == 0, only A is used if present, and'
                'cell vectors are read from card CELL_PARAMETERS.',
                'type': 'REAL'},
'cosBC' : {'description' : 'Specifies the cosine of the angle between lattice vectors B and C',
       'info' : 'Traditional crystallographic constants:'
                'A,B,C in ANGSTROM'
                'cosAB = cosine of the angle between axis a and b (gamma)'
                'cosAC = cosine of the angle between axis a and c (beta)'
                'cosBC = cosine of the angle between axis b and c (alpha)'
                'The axis are chosen according to the value of ibrav.'
                'Specify either these OR celldm but NOT both.'
                'Only needed values (depending on ibrav) must be specified.'
                'The lattice parameter alat = A (in ANGSTROM ).'
                'If ibrav == 0, only A is used if present, and'
                'cell vectors are read from card CELL_PARAMETERS.',
                'type': 'REAL'},
 'ace': {'default': 'true',
         'description': 'Adaptive coordinate exchange',
         'info': ' Use Adaptively Compressed Exchange operator as in Lin Lin, '
                 'J. Chem. Theory Comput. 2016, 12, 2242--2249, '
                 'doi:10.1021/acs.jctc.6b00092  Set to false to use standard '
                 'Exchange (much slower) }',
         'input_type': 'select_multiple',
         'options': ['.FALSE.', '.TRUE.'],
         'type': 'LOGICAL'},
 'assume_isolated': {'default': 'none',
                     'description': 'Assume isolated option',
                     'info': ' Used to perform calculation assuming the system '
                             'to be isolated (a molecule or a cluster in a 3D '
                             'supercell).  Currently available choices: }',
                     'input_type': 'select_multiple',
                     'options': ['none', 'esm', '2D'],
                     'type': 'CHARACTER'},
 'constrained_magnetization': {'default': 'none',
                               'description': 'Constrained magnetization',
                               'info': ' Used to perform constrained '
                                       'calculations in magnetic systems. '
                                       'Currently available choices: } N.B.: '
                                       'symmetrization may prevent to reach '
                                       'the desired orientation of the '
                                       'magnetization. Try not to start with '
                                       'very highly symmetric configurations '
                                       'or use the nosym flag (only as a last '
                                       'remedy) }',
                               'input_type': 'select_multiple',
                               'options': ['none', 'total', 'atomic'],
                               'type': 'CHARACTER'},
 'degauss': {'default': '0.D0 Ry',
             'description': 'Gaussian broadening',
             'info': ' value of the gaussian spreading (Ry) for brillouin-zone '
                     'integration in metals. }',
             'type': 'REAL'},
 'degauss_cond': {'default': '0.D0 Ry',
                  'description': 'Conditional degauss',
                  'info': ' value of the gaussian spreading (Ry) for '
                          'brillouin-zone integration in the conduction '
                          'manifold in a two-chemical potential calculation '
                          '(@ref twochem=.true.). }',
                  'type': 'REAL'},
 'dftd3_threebody': {'default': 'TRUE',
                     'description': 'DFT-D3 three-body dispersion correction',
                     'info': ' Turn three-body terms in Grimme-D3 on. If '
                             '.false. two-body contributions only are '
                             'computed, using two-body parameters of '
                             'Grimme-D3. If dftd3_version=2, three-body '
                             'contribution is always disabled. }',
                     'input_type': 'select_multiple',
                     'options': ['.FALSE.', '.TRUE.'],
                     'type': 'LOGICAL'},
 'dftd3_version': {'default': '3',
                   'description': 'DFT-D3 version',
                   'info': 'Version of Grimme implementation of Grimme-D3:'
                        'dftd3_version = 2 : Original Grimme-D2 parametrization'
                        'dftd3_version = 3 : Grimme-D3 (zero damping)'
                        'dftd3_version = 4 : Grimme-D3 (BJ damping)'
                        'dftd3_version = 5 : Grimme-D3M (zero damping)'
                        'dftd3_version = 6 : Grimme-D3M (BJ damping)',
                   'input_type': 'select_multiple',
                   'options': ['2', '3', '4', '5', '6'],
                   'type': 'INTEGER'},
 'dmft': {'default': '.FALSE.',
          'description': 'Dynamical mean field theory',
          'info': ' If true, nscf calculation will exit in restart mode, scf '
                  'calculation will restart from there if DMFT updates are '
                  'provided as hdf5 archive. Scf calculation should be used '
                  'only with @ref electron_maxstep = 1. @ref K_POINTS have to '
                  'be identical and given explicitly with @ref nosym. }',
          'input_type': 'select_multiple',
          'options': ['.FALSE.', '.TRUE.'],
          'status': ' Requires compilation with hdf5 support }',
          'type': 'LOGICAL'},
 'dmft_prefix': {'default': '@ref prefix',
                 'description': 'DMFT prefix',
                 'info': ' prepended to hdf5 archive: dmft_prefix.h5  DMFT '
                         'update should be provided in group/dataset as: - '
                         'dft_misc_input/band_window with dimension [1, number '
                         'of k-points, 2 (real + complex)] - '
                         'dft_update/delta_N with dimension [number of '
                         'k-points, number of correlated orbitals, number of '
                         'correlated orbitals, 2 (real + complex)] }',
                 'type': 'CHARACTER'},
 'eamp': {'default': '0.001 a.u.',
          'description': 'Electric field amplitude',
          'info': ' Amplitude of the electric field, in ***Hartree*** a.u.; 1 '
                  'a.u. = 51.4220632*10^10 V/m. Used only if @ref '
                  'tefield==.TRUE. The saw-like potential increases with slope '
                  '@ref eamp in the region from (@ref emaxpos+@ref eopreg-1) '
                  'to (@ref emaxpos), then decreases to 0 until (@ref '
                  'emaxpos+@ref eopreg), in units of the crystal vector @ref '
                  'edir. Important: the change of slope of this potential must '
                  'be located in the empty region, or else unphysical forces '
                  'will result. }',
          'type': 'REAL'},
 'ecfixed': {'default': '0.0',
             'description': 'Fixed total energy',
             'type': 'REAL'},
 'ecutfock': {'default': 'ecutrho',
              'description': 'Fock operator cutoff energy',
              'info': ' Kinetic energy cutoff (Ry) for the exact exchange '
                      'operator in EXX type calculations. By default this is '
                      'the same as @ref ecutrho but in some EXX calculations, '
                      'a significant speed-up can be obtained by reducing '
                      'ecutfock, at the expense of some loss in accuracy. Must '
                      'be .gt. @ref ecutwfc. Not implemented for stress '
                      'calculation and for US-PP and PAW pseudopotentials. Use '
                      'with care, especially in metals where it may give raise '
                      'to instabilities. }',
              'type': 'REAL'},
 'ecutrho': {'default': '4 * @ref ecutwfc',
             'description': 'Charge density cutoff energy',
             'info': ' Kinetic energy cutoff (Ry) for charge density and '
                     'potential For norm-conserving pseudopotential you should '
                     'stick to the default value, you can reduce it by a '
                     'little but it will introduce noise especially on forces '
                     'and stress. If there are ultrasoft PP, a larger value '
                     'than the default is often desirable (ecutrho = 8 to 12 '
                     'times @ref ecutwfc, typically). PAW datasets can often '
                     'be used at 4*@ref ecutwfc, but it depends on the shape '
                     'of augmentation charge: testing is mandatory. The use of '
                     'gradient-corrected functional, especially in cells with '
                     'vacuum, or for pseudopotential without non-linear core '
                     'correction, usually requires an higher values of ecutrho '
                     'to be accurately converged. }',
             'type': 'REAL'},
 'ecutvcut': {'default': '0.0 Ry',
              'description': 'Cut-off for vacuum',
              'info': ' Reciprocal space cutoff for correcting Coulomb '
                      'potential divergencies at small q vectors. }',
              'type': 'REAL'},
 'ecutwfc': {'description': 'Wavefunction cutoff energy',
             'info': ' kinetic energy cutoff (Ry) for wavefunctions }',
             'status': 'REQUIRED',
             'type': 'REAL'},
 'edir': {'description': 'Electric field direction',
          'info': ' The direction of the electric field or dipole correction '
                  'is parallel to the bg(:,edir) reciprocal lattice vector, so '
                  'the potential is constant in planes defined by FFT grid '
                  'points; @ref edir = 1, 2 or 3. Used only if @ref tefield is '
                  '.TRUE. }',
          'input_type': 'select_multiple',
          'options': ['1', '2', '3'],
          'type': 'INTEGER'},
 'emaxpos': {'default': '0.5D0',
             'description': 'Maximum positions',
             'info': ' Position of the maximum of the saw-like potential along '
                     'crystal axis @ref edir, within the  unit cell (see '
                     'below), 0 < emaxpos < 1 Used only if @ref tefield is '
                     '.TRUE. }',
             'type': 'REAL'},
 'ensemble_energies': {'default': '.false.',
                       'description': 'Ensemble energies',
                       'info': ' If @ref ensemble_energies = .true., an '
                               'ensemble of xc energies is calculated '
                               'non-selfconsistently for perturbed '
                               'exchange-enhancement factors and LDA vs. PBE '
                               'correlation ratios after each converged '
                               'electronic ground state calculation.  Ensemble '
                               "energies can be analyzed with the 'bee' "
                               'utility included with libbeef.  Requires '
                               'linking against libbeef. @ref input_dft must '
                               'be set to a BEEF-type functional (e.g. '
                               "input_dft = 'BEEF-vdW') }",
                       'input_type': 'select_multiple',
                       'options': ['.FALSE.', '.TRUE.'],
                       'type': 'LOGICAL'},
 'eopreg': {'default': '0.1D0',
            'description': 'External pressure regularization',
            'info': ' Zone in the unit cell where the saw-like potential '
                    'decreases. ( see below, 0 < eopreg < 1 ). Used only if '
                    '@ref tefield is .TRUE. }',
            'type': 'REAL'},
 'esm_bc': {'default': 'bc1',
            'description': 'ESM boundary conditions',
            'info': " If @ref assume_isolated = 'esm', determines the boundary "
                    'conditions used for either side of the slab.  Currently '
                    'available choices: }',
            'input_type': 'select_multiple',
            'options': ['bc1', 'bc2', 'bc3'],
            'type': 'CHARACTER'},
 'esm_efield': {'default': '0.d0',
                'description': 'ESM electric field',
                'info': " If @ref assume_isolated = 'esm' and @ref esm_bc = "
                        "'bc2', gives the magnitude of the electric field "
                        '[Ry/a.u.] to be applied between semi-infinite ESM '
                        'electrodes. }',
                'type': 'REAL'},
 'esm_nfit': {'default': '4',
              'description': 'ESM number of fitting coefficients',
              'info': " If @ref assume_isolated = 'esm', gives the number of "
                      'z-grid points for the polynomial fit along the cell '
                      'edge. }',
              'type': 'INTEGER'},
 'esm_w': {'default': '0.d0',
           'description': 'ESM weight',
           'info': " If @ref assume_isolated = 'esm', determines the position "
                   'offset [in a.u.] of the start of the effective screening '
                   'region, measured relative to the cell edge. (ESM region '
                   'begins at z = +/- [L_z/2 + esm_w] ). }',
           'type': 'REAL'},
 'exx_fraction': {'default': 'it depends on the specified functional',
                  'description': 'Fraction of exact exchange',
                  'info': ' Fraction of EXX for hybrid functional '
                          "calculations. In the case of @ref input_dft='PBE0', "
                          'the default value is 0.25, while for @ref '
                          "input_dft='B3LYP' the @ref exx_fraction default "
                          'value is 0.20. }',
                  'type': 'REAL'},
 'exxdiv_treatment': {'default': 'gygi-baldereschi',
                      'description': 'Exchange-divergent treatment',
                      'info': ' Specific for EXX. It selects the kind of '
                              'approach to be used for treating the Coulomb '
                              'potential divergencies at small q vectors. }',
                      'input_type': 'select_multiple',
                      'options': ['gygi-baldereschi',
                                  'vcut_spherical',
                                  'vcut_ws',
                                  'none'],
                      'type': 'CHARACTER'},
 'force_symmorphic': {'default': '.FALSE.',
                      'description': 'Force symmorphic option',
                      'info': ' if (.TRUE.) force the symmetry group to be '
                              'symmorphic by disabling symmetry operations '
                              'having an associated fractionary translation }',
                      'input_type': 'select_multiple',
                      'options': ['.FALSE.', '.TRUE.'],
                      'type': 'LOGICAL'},
 'gcscf_beta': {'default': '0.05D0',
                'description': 'Beta for gc scf',
                'info': ' Mixing factor for GC-SCF. Larger values are '
                        'recommended, if systems with small DOS on Fermi '
                        'surface as graphite. }',
                'type': 'REAL'},
 'gcscf_conv_thr': {'default': '1.D-2',
                    'description': 'GC-SCF convergence threshold',
                    'info': ' Convergence threshold of Fermi energy (eV) for '
                            'GC-SCF. }',
                    'type': 'REAL'},
 'gcscf_mu': {'description': 'Mu for gc scf',
              'info': ' The target Fermi energy (eV) of GC-SCF. One can start '
                      'with appropriate total charge of the system by giving '
                      '@ref tot_charge . }',
              'status': 'REQUIRED',
              'type': 'REAL'},
 'input_dft': {'default': 'read from pseudopotential files',
               'description': 'DFT type',
               'info': " Exchange-correlation functional: eg 'PBE', 'BLYP' etc "
                       'See Modules/funct.f90 for allowed values. Overrides '
                       'the value read from pseudopotential files. Use with '
                       'care and if you know what you are doing! }',
               'type': 'CHARACTER'},
 'lambda': {'default': '1.d0',
            'description': 'Lambda',
            'info': ' parameter used for constrained_magnetization '
                    'calculations N.B.: if the scf calculation does not '
                    'converge, try to reduce lambda to obtain convergence, '
                    'then restart the run with a larger lambda }',
            'type': 'REAL'},
 'lforcet': {'description': 'Local forces treatment',
             'info': ' When starting a non collinear calculation using an '
                     'existing density file from a collinear lsda calculation '
                     'assumes previous density points in @i z direction and '
                     'rotates it in the direction described by @ref angle1 and '
                     '@ref angle2 variables for atomic type 1 }',
             'input_type': 'select_multiple',
             'options': ['.FALSE.', '.TRUE.'],
             'type': 'LOGICAL'},
 'lgcscf': {'default': '.FALSE.',
            'description': 'Use of GC-SCF',
            'info': ' If .TRUE. perform a constant bias potential '
                    '(constant-mu) calculation with Grand-Canonical SCF. (JCP '
                    '146, 114104 (2017), R.Sundararaman, et al.)  NB: - The '
                    'total energy displayed in output includes the '
                    'potentiostat contribution (-mu*N). - @ref assume_isolated '
                    "= 'esm' and @ref esm_bc = 'bc2' or 'bc3' must be set in "
                    '@ref SYSTEM namelist. - ESM-RISM is also supported (@ref '
                    "assume_isolated = 'esm' and @ref esm_bc = 'bc1' and @ref "
                    "trism = .TRUE.). - @ref mixing_mode has to be 'TF' or "
                    "'local-TF', also its default is 'TF.' - The default of "
                    '@ref mixing_beta is 0.1 with ESM-RISM, 0.2 without '
                    'ESM-RISM. - The default of @ref diago_thr_init is 1.D-5. '
                    '- @ref diago_full_acc is always .TRUE. . - @ref '
                    'diago_rmm_conv is always .TRUE. . }',
            'input_type': 'select_multiple',
            'options': ['.FALSE.', '.TRUE.'],
            'type': 'LOGICAL'},
 'localization_thr': {'description': 'Localization threshold',
                      'info': ' Overlap threshold over which the exchange '
                              'integral over a pair of localized orbitals is '
                              'included in the evaluation of EXX operator. Any '
                              'value greater than 0.0 triggers the SCDM '
                              'localization and the evaluation on EXX using '
                              'the localized orbitals. Very small value of the '
                              'threshold should yield the same result as the '
                              'default EXX evaluation }',
                      'type': 'REAL'},
 'london': {'default': '.FALSE.',
            'description': 'London dispersion correction',
            'input_type': 'select_multiple',
            'options': ['.FALSE.', '.TRUE.'],
            'status': " OBSOLESCENT, same as @ref vdw_corr='DFT-D' }",
            'type': 'LOGICAL'},
 'london_rcut': {'default': '200',
                 'description': 'London dispersion correction cutoff radius',
                 'info': ' cutoff radius (a.u.) for dispersion interactions }',
                 'type': 'REAL'},
 'london_s6': {'default': '0.75',
               'description': 'London dispersion correction s6 coefficient',
               'info': ' global scaling parameter for DFT-D. Default is good '
                       'for PBE. }',
               'type': 'REAL'},
 'lspinorb': {'description': 'Spin-orbit coupling',
              'info': ' if .TRUE. the noncollinear code can use a '
                      'pseudopotential with spin-orbit. }',
              'input_type': 'select_multiple',
              'options': ['.FALSE.', '.TRUE.'],
              'type': 'LOGICAL'},
 'nat': {'description': 'Number of atoms',
         'info': ' number of atoms in the unit cell (ALL atoms, except if '
                 'space_group is set, in which case, INEQUIVALENT atoms) }',
         'status': 'REQUIRED',
         'type': 'INTEGER'},
 'nbnd': {'description': 'Number of Kohn-Sham bands',
          'info': ' Number of electronic states (bands) to be calculated. Note '
                  'that in spin-polarized calculations the number of k-point, '
                  'not the number of bands per k-point, is doubled }',
          'type': 'INTEGER'},
 'nbnd_cond': {'description': 'Number of bands for conductivity',
               'info': ' Number of electronic states in the conduction '
                       'manifold for a two chemical-potential calculation '
                       '(@ref twochem=.true.). }',
               'type': 'INTEGER'},
 'nelec_cond': {'default': '0.D0',
                'description': 'Number of electrons for conductivity',
                'info': ' Number of electrons placed in the conduction '
                        'manifold in a two-chemical potential calculation '
                        '(@ref twochem=.true.). Of the total # of electrons '
                        'nelec, nelec-nelec_cond will occupy the valence '
                        'manifold and nelec_cond will be constrained in the '
                        'conduction manifold. }',
                'type': 'REAL'},
 'nextffield': {'default': '0',
                'description': 'Next field',
                'info': ' Number of activated external ionic force fields. See '
                        'Doc/ExternalForceFields.tex for further explanation '
                        'and parameterizations }',
                'type': 'INTEGER'},
 'no_t_rev': {'default': '.FALSE.',
              'description': 'Time reversal symmetry suppression',
              'info': ' if (.TRUE.) disable the usage of magnetic symmetry '
                      'operations that consist in a rotation + time reversal. '
                      '}',
              'input_type': 'select_multiple',
              'options': ['.FALSE.', '.TRUE.'],
              'type': 'LOGICAL'},
 'noinv': {'default': '.FALSE.',
           'description': 'Inversion symmetry suppression',
           'info': ' if (.TRUE.) disable the usage of k => -k symmetry (time '
                   'reversal) in k-point generation }',
           'input_type': 'select_multiple',
           'options': ['.FALSE.', '.TRUE.'],
           'type': 'LOGICAL'},
 'noncolin': {'default': '.false.',
              'description': 'Non-collinear magnetization',
              'info': ' if .true. the program will perform a noncollinear '
                      'calculation. }',
              'input_type': 'select_multiple',
              'options': ['.FALSE.', '.TRUE.'],
              'type': 'LOGICAL'},
 'nosym': {'default': '.FALSE.',
           'description': 'Symmetry operation suppression',
           'info': ' if (.TRUE.) symmetry is not used. Consequences:  - if a '
                   'list of k points is provided in input, it is used "as is": '
                   'symmetry-inequivalent k-points are not generated, and the '
                   'charge density is not symmetrized;  - if a uniform '
                   '(Monkhorst-Pack) k-point grid is provided in input, it is '
                   'expanded to cover the entire Brillouin Zone, irrespective '
                   'of the crystal symmetry. Time reversal symmetry is assumed '
                   'so k and -k are considered as equivalent unless @ref '
                   'noinv=.true. is specified.  Do not use this option unless '
                   'you know exactly what you want and what you get. May be '
                   'useful in the following cases: - in low-symmetry large '
                   'cells, if you cannot afford a k-point grid with the '
                   'correct symmetry - in MD simulations - in calculations for '
                   'isolated atoms }',
           'input_type': 'select_multiple',
           'options': ['.FALSE.', '.TRUE.'],
           'type': 'LOGICAL'},
 'nosym_evc': {'default': '.FALSE.',
               'description': 'Symmetry operation suppression for '
                              'wavefunctions',
               'info': ' if (.TRUE.) symmetry is not used, and k points are '
                       'forced to have the symmetry of the Bravais lattice; an '
                       'automatically generated Monkhorst-Pack grid will '
                       'contain all points of the grid over the entire '
                       'Brillouin Zone, plus the points rotated by the '
                       'symmetries of the Bravais lattice which were not in '
                       'the original grid. The same applies if a k-point list '
                       'is provided in input instead of a Monkhorst-Pack grid. '
                       'Time reversal symmetry is assumed so k and -k are '
                       'equivalent unless @ref noinv=.true. is specified. This '
                       'option differs from @ref nosym because it forces '
                       'k-points in all cases to have the full symmetry of the '
                       'Bravais lattice (not all uniform grids have such '
                       'property!) }',
               'input_type': 'select_multiple',
               'options': ['.FALSE.', '.TRUE.'],
               'type': 'LOGICAL'},
 'nqx1' : { 'description' : 'Number of grid points for charge density in the first reciprocal lattice vector',
             'info' : 'Three-dimensional mesh for q (k1-k2) sampling of '
                'the Fock operator (EXX). Can be smaller than '
                'the number of k-points.'
                'Currently this defaults to the size of the k-point mesh used.'
                'In QE =< 5.0.2 it defaulted to nqx1=nqx2=nqx3=1.',
             'type' : 'INTEGER'},
 'nqx2' : { 'description' : 'Number of grid points for charge density in the second reciprocal lattice vector',
             'info' : 'Three-dimensional mesh for q (k1-k2) sampling of '
                'the Fock operator (EXX). Can be smaller than '
                'the number of k-points.'
                'Currently this defaults to the size of the k-point mesh used.'
                'In QE =< 5.0.2 it defaulted to nqx1=nqx2=nqx3=1.',
             'type' : 'INTEGER'},
 'nqx3' : { 'description' : 'Number of grid points for charge density in the third reciprocal lattice vector',
             'info' : 'Three-dimensional mesh for q (k1-k2) sampling of '
                'the Fock operator (EXX). Can be smaller than '
                'the number of k-points.'
                'Currently this defaults to the size of the k-point mesh used.'
                'In QE =< 5.0.2 it defaulted to nqx1=nqx2=nqx3=1.',
             'type' : 'INTEGER'},
 'nr1s' : { 'description' : 'Number of grid points in the first reciprocal lattice vector',
             'info' : 'Three-dimensional mesh for wavefunction FFT and for the smooth'
                'part of charge density ( smooth grid ).'
                'Coincides with @ref nr1, @ref nr2, @ref nr3 if @ref ecutrho = 4 * ecutwfc',
             'type' : 'INTEGER'},
 'nr2s' : { 'description' : 'Number of grid points in the first reciprocal lattice vector',
             'info' : 'Three-dimensional mesh for wavefunction FFT and for the smooth'
                'part of charge density ( smooth grid ).'
                'Coincides with @ref nr1, @ref nr2, @ref nr3 if @ref ecutrho = 4 * ecutwfc',
             'type' : 'INTEGER'},
 'nr3s' : { 'description' : 'Number of grid points in the first reciprocal lattice vector',
             'info' : 'Three-dimensional mesh for wavefunction FFT and for the smooth'
                'part of charge density ( smooth grid ).'
                'Coincides with @ref nr1, @ref nr2, @ref nr3 if @ref ecutrho = 4 * ecutwfc',
             'type' : 'INTEGER'},
 'nspin': {'default': '1',
           'description': 'Spin polarized calculation',
           'info': ' nspin = 1 :  non-polarized calculation (default)  nspin = '
                   '2 :  spin-polarized calculation, LSDA (magnetization along '
                   'z axis)  nspin = 4 :  spin-polarized calculation, '
                   'noncollinear (magnetization in generic direction) DO NOT '
                   'specify @ref nspin in this case; specify @ref '
                   'noncolin=.TRUE. instead }',
           'input_type': 'select_multiple',
           'options' : ['1','2','4'],
           'type': 'INTEGER'},
 'ntyp': {'description': 'Number of atom types',
          'info': ' number of types of atoms in the unit cell }',
          'status': 'REQUIRED',
          'type': 'INTEGER'},
 'occupations': {'description': 'Occupations',
                 'info': '',
                 'input_type': 'select_multiple',
                 'options': ['smearing',
                             'tetrahedra',
                             'tetrahedra_lin',
                             'tetrahedra_opt',
                             'fixed',
                             'from_input'],
                 'type': 'CHARACTER'},
 'one_atom_occupations': {'default': '.FALSE.',
                          'description': 'One atom occupations',
                          'info': ' This flag is used for isolated atoms (@ref '
                                  'nat=1) together with @ref '
                                  "occupations='from_input'. If it is .TRUE., "
                                  'the wavefunctions are ordered as the atomic '
                                  'starting wavefunctions, independently from '
                                  'their eigenvalue. The occupations indicate '
                                  'which atomic states are filled.  The order '
                                  'of the states is written inside the UPF '
                                  'pseudopotential file. In the scalar '
                                  'relativistic case: S -> l=0, m=0 P -> l=1, '
                                  'z, x, y D -> l=2, r^2-3z^2, xz, yz, xy, '
                                  'x^2-y^2  In the noncollinear magnetic case '
                                  '(with or without spin-orbit), each group of '
                                  'states is doubled. For instance: P -> l=1, '
                                  'z, x, y for spin up, l=1, z, x, y for spin '
                                  'down. Up and down is relative to the '
                                  'direction of the starting magnetization.  '
                                  'In the case with spin-orbit and '
                                  'time-reversal (@ref '
                                  'starting_magnetization=0.0) the atomic '
                                  'wavefunctions are radial functions '
                                  'multiplied by spin-angle functions. For '
                                  'instance: P -> l=1, j=1/2, m_j=-1/2,1/2. '
                                  'l=1, j=3/2, m_j=-3/2, -1/2, 1/2, 3/2.  In '
                                  'the magnetic case with spin-orbit the '
                                  'atomic wavefunctions can be forced to be '
                                  'spin-angle functions by setting @ref '
                                  'starting_spin_angle to .TRUE.. }',
                          'input_type': 'select_multiple',
                          'options': ['.FALSE.', '.TRUE.'],
                          'type': 'LOGICAL'},
 'origin_choice': {'default': '1',
                   'description': 'Origin choice',
                   'info': 'Used only for space groups that in the ITA allow '
                           'the use of two different origins. @ref '
                           'origin_choice=1, means the first origin, while '
                           '@ref origin_choice=2 is the second origin. }',
                   'input_type': 'select_multiple',
                   'options': ['1', '2'],
                   'type': 'INTEGER'},
 'pol_type': {'description': 'Polarization type',
              'info': ' Type of polaron in gammaDFT. }',
              'input_type': 'select_multiple',
              'options': ['e', 'h'],
              'type': 'CHARACTER'},
 'q2sigma': {'default': '0.1',
             'description': 'Sigma for broadening',
             'info': ' ecfixed, qcutz, q2sigma:  parameters for modified '
                     'functional to be used in variable-cell molecular '
                     'dynamics (or in stress calculation). "ecfixed" is the '
                     'value (in Rydberg) of the constant-cutoff; "qcutz" and '
                     '"q2sigma" are the height and the width (in Rydberg) of '
                     'the energy step for reciprocal vectors whose square '
                     'modulus is greater than "ecfixed". In the kinetic '
                     'energy, G^2 is replaced by G^2 + qcutz * (1 + erf ( (G^2 '
                     '- ecfixed)/q2sigma) ) See: M. Bernasconi et al, J. Phys. '
                     'Chem. Solids 56, 501 (1995), '
                     'doi:10.1016/0022-3697(94)00228-2 }',
             'type': 'REAL'},
 'qcutz': {'default': '0.0',
           'description': 'Zone sampling cut-off',
           'type': 'REAL'},
 'report': {'default': '-1',
            'description': 'Report',
            'info': ' determines when atomic magnetic moments are printed on '
                    'output: @b {report = 0}  never',
            'type': 'INTEGER'},
 'rhombohedral': {'default': '.TRUE.',
                  'description': 'Rhombohedral',
                  'info': ' Used only for rhombohedral space groups. When '
                          '.TRUE. the coordinates of the inequivalent atoms '
                          'are given with respect to the rhombohedral axes, '
                          'when .FALSE. the coordinates of the inequivalent '
                          'atoms are given with respect to the hexagonal axes. '
                          'They are converted internally to the rhombohedral '
                          'axes and @ref ibrav=5 is used in both cases. }',
                  'input_type': 'select_multiple',
                  'options': ['.FALSE.', '.TRUE.'],
                  'type': 'LOGICAL'},
 'sci_cb': {'default': '0',
            'description': 'Conduction band spin',
            'info': ' Conduction band band shift (in eV) through '
                    'self-consistent scissor operator. When performing '
                    'gammaDFT calculations of polarons, the polaron level is '
                    'not shifted. }',
            'type': 'REAL'},
 'sci_vb': {'default': '0',
            'description': 'Valence band spin',
            'info': ' Valence band shift (in eV) through self-consistent '
                    'scissor operator. When performing gammaDFT calculations '
                    'of polarons, the polaron level is not shifted. }',
            'type': 'REAL'},
 'screening_parameter': {'default': '0.106',
                         'description': 'Screening parameter',
                         'info': ' screening_parameter for HSE like hybrid '
                                 'functionals. For more information, see: J. '
                                 'Chem. Phys. 118, 8207 (2003), '
                                 'doi:10.1063/1.1564060 J. Chem. Phys. 124, '
                                 '219906 (2006), doi:10.1063/1.2204597 }',
                         'type': 'REAL'},
 'sic_energy': {'default': '.false.',
                'description': 'SIC energy',
                'info': ' Enable the calculation of the total energy in '
                        'gammaDFT. When .true., a preliminary calculation is '
                        'performed to calculate the electron density in the '
                        'absence of the polaron. When .false., the total '
                        'energy printed in output should not be considered. '
                        'For structural relaxations, it is recommended to use '
                        '.false. to avoid doubling the computational cost. }',
                'input_type': 'select_multiple',
                'options': ['.FALSE.', '.TRUE.'],
                'type': 'LOGICAL'},
 'sic_gamma': {'default': '0',
               'description': 'SIC gamma',
               'info': ' Strength of the gammaDFT potential. }',
               'type': 'REAL'},
 'smearing': {'default': 'gaussian',
              'description': 'Smearing technique',
              'info': ' Available options are: }',
              'type': 'CHARACTER'},
 'space_group': {'default': '0',
                 'description': 'Space group',
                 'info': ' The number of the space group of the crystal, as '
                         'given in the International Tables of Crystallography '
                         'A (ITA). This allows to give in input only the '
                         'inequivalent atomic positions. The positions of all '
                         'the symmetry equivalent atoms are calculated by the '
                         'code. Used only when the atomic positions are of '
                         'type crystal_sg. See also @ref uniqueb, @ref '
                         'origin_choice, @ref rhombohedral }',
                 'type': 'INTEGER'},
 'starting_spin_angle': {'default': '.FALSE.',
                         'description': 'Starting spin angle',
                         'info': ' In the spin-orbit case when @ref '
                                 'domag=.TRUE., by default, the starting '
                                 'wavefunctions are initialized as in scalar '
                                 'relativistic noncollinear case without '
                                 'spin-orbit.  By setting @ref '
                                 'starting_spin_angle=.TRUE. this behaviour '
                                 'can be changed and the initial wavefunctions '
                                 'are radial functions multiplied by '
                                 'spin-angle functions.  When @ref '
                                 'domag=.FALSE. the initial wavefunctions are '
                                 'always radial functions multiplied by '
                                 'spin-angle functions independently from this '
                                 'flag.  When @ref lspinorb is .FALSE. this '
                                 'flag is not used. }',
                         'input_type': 'select_multiple',
                         'options': ['.FALSE.', '.TRUE.'],
                         'type': 'LOGICAL'},
 'tot_charge': {'default': '0.0',
                'description': 'Total charge',
                'info': ' Total charge of the system. Useful for simulations '
                        'with charged cells. By default the unit cell is '
                        'assumed to be neutral (tot_charge=0). tot_charge=+1 '
                        'means one electron missing from the system, '
                        'tot_charge=-1 means one additional electron, and so '
                        'on.  In a periodic calculation a compensating jellium '
                        'background is inserted to remove divergences if the '
                        'cell is not neutral. }',
                'type': 'REAL'},
 'tot_magnetization': {'default': '-10000 [unspecified]',
                       'description': 'Total magnetization',
                       'info': ' Total majority spin charge - minority spin '
                               'charge. Used to impose a specific total '
                               'electronic magnetization. If unspecified then '
                               'tot_magnetization variable is ignored and the '
                               'amount of electronic magnetization is '
                               'determined during the self-consistent cycle. }',
                       'type': 'REAL'},
 'ts_vdw_econv_thr': {'default': '1.D-6',
                      'description': 'TS-vdW energy convergence threshold',
                      'info': ' Optional: controls the convergence of the vdW '
                              'energy (and forces). The default value is a '
                              'safe choice, likely too safe, but you do not '
                              'gain much in increasing it }',
                      'type': 'REAL'},
 'ts_vdw_isolated': {'default': '.FALSE.',
                     'description': 'TS-vdW isolated',
                     'info': ' Optional: set it to .TRUE. when computing the '
                             'Tkatchenko-Scheffler vdW energy or the Many-Body '
                             'dispersion (MBD) energy for an isolated '
                             '(non-periodic) system. }',
                     'input_type': 'select_multiple',
                     'options': ['.FALSE.', '.TRUE.'],
                     'type': 'LOGICAL'},
 'uniqueb': {'default': '.FALSE.',
             'description': 'Unique B',
             'info': ' Used only for monoclinic lattices. If .TRUE. the b '
                     'unique @ref ibrav (-12 or -13) are used, and symmetry '
                     'equivalent positions are chosen assuming that the '
                     'twofold axis or the mirror normal is parallel to the b '
                     'axis. If .FALSE. it is parallel to the c axis. }',
             'input_type': 'select_multiple',
             'options': ['.FALSE.', '.TRUE.'],
             'type': 'LOGICAL'},
 'use_all_frac': {'default': '.FALSE.',
                  'description': 'Use of all fractional coordinates',
                  'info': ' if (.FALSE.) force real-space FFT grids to be '
                          'commensurate with fractionary translations of '
                          'non-symmorphic symmetry operations, if present '
                          '(e.g.: if a fractional translation (0,0,c/4) '
                          'exists, the FFT dimension along the c axis must be '
                          'multiple of 4). if (.TRUE.) do not impose any '
                          'constraints to FFT grids, even in the presence of '
                          'non-symmorphic symmetry operations. BEWARE: '
                          'use_all_frac=.TRUE. may lead to wrong results for '
                          'hybrid functionals and phonon calculations. Both '
                          'cases use symmetrization in real space that works '
                          'for non-symmorphic operations only if the '
                          'real-space FFT grids are commensurate. }',
                  'input_type': 'select_multiple',
                  'options': ['.FALSE.', '.TRUE.'],
                  'type': 'LOGICAL'},
 'vdw_corr': {'default': 'none',
              'description': 'Van der Waals correction',
              'info': '\'grimme-d2\', \'Grimme-D2\', \'DFT-D\', \'dft-d\' : Semiempirical Grimme\'s DFT-D2.'
                        '\'grimme-d3\', \'Grimme-D3\', \'DFT-D3\', \'dft-d3\' : Semiempirical Grimme\'s DFT-D3.'
                        '\'TS\', \'ts\', \'ts-vdw\', \'ts-vdW\', \'tkatchenko-scheffler\' : Tkatchenko-Scheffler dispersion corrections with first-principle derived C6 coefficients.'
                        '\'MBD\', \'mbd\', \'many-body-dispersion\', \'mbd_vdw\' : Many-body dipersion (MBD) correction to long-range interactions.'
                        '\'XDM\', \'xdm\' : Exchange-hole dipole-moment model.',
              'input_type': 'select_multiple',
              'options': ['grimme-d2', 'Grimme-D2', 'DFT-D','dft-d','grimme-d3', 
                          'Grimme-D3', 'DFT-D3','dft-d3', 'TS','ts','ts-vdw', 
                          'ts-vdW','tkatchenko-scheffler','MBD', 'mbd', 
                          'many-body-dispersion', 'mbd_vdw','XDM', 'xdm'],
              'type': 'CHARACTER'},
 'x_gamma_extrapolation': {'description': 'Extrapolation in Gamma',
                           'info': ' Specific for EXX. If .true., extrapolate '
                                   'the G=0 term of the potential (see README '
                                   'in examples/EXX_example for more) Set this '
                                   'to .false. for GAU-PBE. }',
                           'input_type': 'select_multiple',
                           'options': ['.FALSE.', '.TRUE.'],
                           'type': 'LOGICAL'},
 'xdm_a1': {'default': '0.6836',
            'description': 'Exchange dispersion correction parameter',
            'info': 'Damping function parameter a1 (adimensional). It is NOT '
                    'necessary to give a value if the functional is one of '
                    'B86bPBE, PW86PBE, PBE, BLYP. For functionals in this '
                    'list, the coefficients are given in: '
                    'http://schooner.chem.dal.ca/wiki/XDM A. Otero de la Roza, '
                    'E. R. Johnson, J. Chem. Phys. 138, 204109 (2013), '
                    'doi:10.1063/1.4705760 }',
            'type': 'REAL'},
 'xdm_a2': {'default': '1.5045',
            'description': 'Exchange dispersion correction parameter A2',
            'info': ' Damping function parameter a2 (angstrom). It is NOT '
                    'necessary to give a value if the functional is one of '
                    'B86bPBE, PW86PBE, PBE, BLYP. For functionals in this '
                    'list, the coefficients are given in: '
                    'http://schooner.chem.dal.ca/wiki/XDM A. Otero de la Roza, '
                    'E. R. Johnson, J. Chem. Phys. 138, 204109 (2013), '
                    'doi:10.1063/1.4705760 }',
            'type': 'REAL'},
 'zgate' : {'default' : '0.5',
            'description' : 'Z-gate parameter',
            'info' : 'used only if @ref gate = .TRUE.'
                    'Specifies the position of the charged plate which represents'
                    'the counter charge in doped systems (@ref tot_charge .ne. 0).'
                    'In units of the unit cell length in @i z direction, @ref zgate in ]0,1['
                    'Details of the gate potential can be found in'
                    'T. Brumme, M. Calandra, F. Mauri; PRB 89, 245406 (2014).',
            'type' : 'REAL'},
 'relaxz' : {'default' : '.FALSE.',
            'description' : 'Relaxation of the z coordinate',
            'input_type': 'select_multiple',
            'options': ['.FALSE.', '.TRUE.'],
            'info' : 'used only if @ref gate = .TRUE.'
                    'Allows the relaxation of the system towards the charged plate.'
                    'Use carefully and utilize either a layer of fixed atoms or a'
                    'potential barrier (@ref block=.TRUE.) to avoid the atoms moving to'
                    'the position of the plate or the dipole of the dipole'
                    'correction (@ref dipfield=.TRUE.).',
            'type' : 'LOGICAL'},
 'block' : {'default' : '.FALSE.',
            'description' : 'Block parameter',
            'input_type': 'select_multiple',
            'options': ['.FALSE.', '.TRUE.'],
            'info' : 'used only if @ref gate = .TRUE.'
                    'Adds a potential barrier to the total potential seen by the'
                    'electrons to mimic a dielectric in field effect configuration'
                    'and/or to avoid electrons spilling into the vacuum region for'
                    'electron dopediring. Potential barrier is from @ref block_1 to @ref block_2 and'
                    'has a height of block_height.'
                    'If @ref dipfield = .TRUE. then @ref eopreg is used for a smooth increase and'
                    'decrease of the potential barrier.',
            'type' : 'LOGICAL'},
 'block_1' : {'default' : '0.45',
            'description' : 'Block parameter 1',
            'info' : 'used only if @ref gate = .TRUE. and @ref block = .TRUE.'
                    'lower beginning of the potential barrier, in units of the'
                    'unit cell size along @i z, @ref block_1 in ]0,1[',
            'type' : 'REAL'},
 'block_2' : {'default' : '0.55',
            'description' : 'Block parameter 2',
            'info' : 'used only if @ref gate = .TRUE. and @ref block = .TRUE.'
                    'upper beginning of the potential barrier, in units of the'
                    'unit cell size along @i z, @ref block_2 in ]0,1[',
            'type' : 'REAL'},
 'block_height' : {'default' : '0.1',
            'description' : 'Block height parameter',
            'info' : 'used only if @ref gate = .TRUE. and @ref block = .TRUE.'
                    'Height of the potential barrier in Rydberg.',
            'type' : 'REAL'}
}