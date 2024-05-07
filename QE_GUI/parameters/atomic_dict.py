INPUT_DATA = {
    'ATOMIC_SPECIES': {
        'input_type': 'text',
        'data': {
            'X': {'mass': 'Mass_X', 'pseudo_potential': 'PseudoPot_X'},
            'Y': {'mass': 'Mass_Y', 'pseudo_potential': 'PseudoPot_Y'},
            'Z': {'mass': 'Mass_Z', 'pseudo_potential': 'PseudoPot_Z'}
        }
    },
    'ATOMIC_POSITIONS': {
        'input_type': 'text',
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
        'input_type': 'matrix',
        'data': {
            'type': 'manual',
            'matrix': [
                [None, None, None],
                [None, None, None]
            ]
        }
    }
}

