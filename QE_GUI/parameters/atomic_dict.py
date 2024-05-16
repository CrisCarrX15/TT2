INPUT_DATA = {
    'ATOMIC_SPECIES': {
        'description' : 'Specifies the atomic species and their pseudopotentials',
        'input_type': 'atomic_species',
        'info' : 'Syntax:\n'
                'X(1)  	 Mass_X(1)  	 PseudoPot_X(1)\n'
                'X(2)  	 Mass_X(2)  	 PseudoPot_X(2)\n'
                '. . .\n'
                'X(ntyp)  	 Mass_X(ntyp)  	 PseudoPot_X(ntyp)\n\n'
                'X : label of the atom. Acceptable syntax:'
                'chemical symbol X (1 or 2 characters, case-insensitive)'
                'or chemical symbol plus a number or a letter, as in'
                '\"Xn\" (e.g. Fe1) or "X_*" or "X-*" (e.g. C1, C_h;'
                'max total length cannot exceed 3 characters)\n'
                'Mass_X : mass of the atomic species [amu: mass of C = 12]'
                'Used only when performing Molecular Dynamics run'
                'or structural optimization runs using Damped MD.'
                'Not actually used in all other cases (but stored'
                'in data files, so phonon calculations will use'
                'these values unless other values are provided)\n'
                'PseudoPot_X : File containing PP for this species.'
    },
    'ATOMIC_POSITIONS': {
        'description' : '{alat} Specifies the atomic positions in Cartesian or crystal coordinates',
        'input_type': 'atomic_positions',
        'info' : 'atomic positions are in cartesian coordinates, in'
                'units of the lattice parameter (either celldm(1)'
                'or A). If no option is specified, \'alat\' is assumed;'
                'not specifying units is DEPRECATED and will no'
                'longer be allowed in the future'
    },
    'K_POINTS': {
        'description' : '{automatic} Specifies the k-points sampling method for Brillouin zone integration',
        'input_type': 'k_points',
        'info' : 'Syntax:\n'
                'K_POINTS automatic\n'
                'nk1  nk2  nk3  sk1  sk2  sk3\n\n'
                'Special k-points (xk_x/y/z) in the irreducible Brillouin Zone'
                '(IBZ) of the lattice (with all symmetries) and weights (wk)'
                'See the literature for lists of special points and'
                'the corresponding weights.'
                'If the symmetry is lower than the full symmetry'
                'of the lattice, additional points with appropriate'
                'weights are generated. Notice that such procedure'
                'assumes that ONLY k-points in the IBZ are provided in input'
                'In a non-scf calculation, weights do not affect the results.'
                'If you just need eigenvalues and eigenvectors (for instance,'
                'for a band-structure plot), weights can be set to any value'
                '(for instance all equal to 1).'
    },
    'CELL_PARAMETERS': {
        'description' : '{alat} Specifies the lattice vectors of the simulation cell in Cartesian coordinates',
        'input_type': 'cell_parameters',
        'info' : 'Syntax:\n'
                'CELL_PARAMETERS alat\n'
                'v1(1)  	 v1(2)  	 v1(3)\n'
                'v2(1)  	 v2(2)  	 v2(3)\n'
                'v3(1)  	 v3(2)  	 v3(3)\n\n'
                'v1, v2, v3\n'
                'Crystal lattice vectors (in cartesian axis):\n'
                'v1(1)  v1(2)  v1(3)    ... 1st lattice vector\n'
                'v2(1)  v2(2)  v2(3)    ... 2nd lattice vector\n'
                'v3(1)  v3(2)  v3(3)    ... 3rd lattice vector\n'
    }
}

