import re
from control_dict import CONTROL_DICT
from system_dict import SYSTEM_DICT

def validate_all_entries(in_parameters):
    error_message = validate_types(in_parameters, {'CONTROL' : CONTROL_DICT, 'SYSTEM' : SYSTEM_DICT})
    validate_rules()

    return error_message


def validate_types(in_parameters, parameters_dict):
    error_message = ''
    for section, params in in_parameters.items():
        for key, value in params.items():
            if value != '' and value is not None:
                input_type = parameters_dict[section][key]['type']

                if input_type == 'REAL' and not is_real_number(value):
                    error_message += f'The parameter {key} needs to be a real number\n'
                elif input_type == 'INTEGER' and not is_integer(value):
                    error_message += f'The parameter {key} needs to be an integer number\n'
                elif input_type == 'LOGICAL' and not is_logical(value):
                    error_message += f'The parameter {key} needs to be .TRUE. or .FALSE.\n'
    return error_message

# En esta función se verificarán las reglas en los parámetros que aplique
def validate_rules():
    return



def is_real_number(string):
    # Expresión regular para un número real
    pattern = r'^[-+]?[0-9]*\.?[0-9]+(?:[eEdD][-+]?[0-9]+)?$'
    # Verificar si la cadena coincide con el patrón
    if re.match(pattern, string):
        return True
    else:
        return False

def is_integer(string):
    # Verificar si la cadena es un número entero
    return string.isdigit()

def is_logical(string):
    if string.lower() == '.false.' or string.lower() == '.true.':
        return True
    else:
        return False