INPUT_DATA = {
    'ATOMIC_SPECIES': {
        'description' : 'Specifies the atomic species and their pseudopotentials',
        'input_type': 'atomic_species',
        'data': {
            'X': {'mass': 'Mass_X', 'pseudo_potential': 'PseudoPot_X'},
            'Y': {'mass': 'Mass_Y', 'pseudo_potential': 'PseudoPot_Y'},
            'Z': {'mass': 'Mass_Z', 'pseudo_potential': 'PseudoPot_Z'}
        }
    },
    'ATOMIC_POSITIONS': {
        'description' : 'Specifies the atomic positions in Cartesian or crystal coordinates',
        'input_type': 'atomic_positions',
        'data': {
            'units': 'angstrom',
            'atoms': [
                {'atom': 'X', 'position': (0.0, 0.0, 0.0)},
                {'atom': 'Y', 'position': (0.5, 0.0, 0.0)},
                {'atom': 'Z', 'position': (0.0, 0.2, 0.2)}
            ]
        }
    },
    'K_POINTS': {
        'description' : 'Specifies the k-points sampling method for Brillouin zone integration',
        'input_type': 'k_points',
        'data': {
            'type': 'manual',
            'matrix': [
                [None, None, None],
                [None, None, None]
            ]
        }
    },
    'CELL_PARAMETERS': {
        'description' : 'Specifies the lattice vectors of the simulation cell in Cartesian coordinates',
        'input_type': 'cell_parameters',
    }
}

