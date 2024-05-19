CONFIG_DICT = { 
        'performance': {
            'description' : 'Select the performance dedicated to Quantum ESPRESSO calculation',
            'info' : 'Standard: The calculation will be executed with one core\n'
                     'High: Will run with maximum cores depending on RAM',
            'input_type': 'select_multiple',
            'options': ['High', 
                        'Standard'],
            },
        'graph_3D': {
            'description' : 'Enable 3D graphics for .in and .out files',
            'info' : '.xyz files will be generated based on the atomic positions of the .in and .out file',
            'input_type': 'select_multiple',
            'options': ['Yes', 
                        'No'],
        }
}