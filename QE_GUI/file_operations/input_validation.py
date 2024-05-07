import re
import os
from parameters.control_dict import CONTROL_DICT
from parameters.system_dict import SYSTEM_DICT
from parameters.electrons_dict import ELECTRONS_DICT
from parameters.ions_dict import IONS_DICT
from parameters.rism_dict import RISM_DICT

##### NOTA: AL INICIAR LA GUI POR PRIMERA VEZ DEBEMOS PREGUNTAR DONDE ESTA GUARDADA LA CARPETA DE qe-7.3


def validate_all_entries(in_parameters):
    error_message = ''
    error_message += validate_required_fields(in_parameters, {'CONTROL' : CONTROL_DICT, 'SYSTEM' : SYSTEM_DICT, 'ELECTRONS' : ELECTRONS_DICT, 'IONS' : IONS_DICT, 'RISM' : RISM_DICT})
    error_message += validate_types(in_parameters, {'CONTROL' : CONTROL_DICT, 'SYSTEM' : SYSTEM_DICT, 'ELECTRONS' : ELECTRONS_DICT, 'IONS' : IONS_DICT, 'RISM' : RISM_DICT})
    #error_message +=validate_rules(in_parameters)

    return error_message



# ===================================================================
# ================== VALIDATION OF REQUIRED FIELDS ==================
# ===================================================================
def validate_required_fields(in_parameters, parameters_dict):
    error_message = ''
    for section, params in parameters_dict.items():
        for key, param_info in params.items():
            # Check if 'status' is present in param_info, if not, by default it is considered ''
            status = param_info.get('status', '')
            if status == 'REQUIRED' and key in in_parameters.get(section, {}):
                value = in_parameters[section][key]
                if value == '' or value is None:
                    error_message += f'The parameter {key} is required but not provided\n'
    return error_message



# ===================================================================
# ======================= VALIDATION OF TYPES =======================
# ===================================================================
def validate_types(in_parameters, parameters_info):
    error_message = ''
    for section, params in in_parameters.items():
        for key, value in params.items():
            if value != '' and value is not None:
                input_type = parameters_info[section][key].get('type', 'UNKNOWN')

                if input_type == 'REAL' and not is_real_number(value):
                    error_message += f'The parameter {key} needs to be a real number\n'
                elif input_type == 'INTEGER' and not is_integer(value):
                    error_message += f'The parameter {key} needs to be an integer number\n'
                elif input_type == 'LOGICAL' and not is_logical(value):
                    error_message += f'The parameter {key} needs to be .TRUE. or .FALSE.\n'
                elif input_type == 'UNKNOWN':
                    error_message += f'The parameter {key} has an unknown data type\n'
    return error_message


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


# ===================================================================
# ======================= VALIDATION OF RULES =======================
# ===================================================================
# En esta función se verificarán las reglas en los parámetros que aplique
def validate_rules(in_parameters):
    control_validation(in_parameters)
    return



# ==== VALIDATION OF RULES ====
def control_validation(control):
    error_message = ''

    if (control['lfcp'] == '.TRUE.'):
        if (control['calulation']!='relax') or (control['assume_isolated']!='esm') or (control['esm_bc']!='bc2' and control['esm_bc']!='bc3'):
            error_message+='So that \'lfcp\' is .TRUE. does not meet the necessary requirements, see the More Info section for more information'

    return error_message