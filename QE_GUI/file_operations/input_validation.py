import re
import os
from parameters.control_dict import CONTROL_DICT
from parameters.system_dict import SYSTEM_DICT
from parameters.electrons_dict import ELECTRONS_DICT
from parameters.ions_dict import IONS_DICT
from parameters.rism_dict import RISM_DICT



class QEInputValidator:
    def __init__(self, data):
        # Data is dict of dict
        self.input_data = data

    # ===================================================================
    # ================== validation of REQUIRED fields ==================
    # ===================================================================
    def validate_required_fields(self, in_parameters, parameters_dict):
        error_message = ''
        for key, param_info in parameters_dict.items():
            # Check if 'status' is present in param_info, if not, by default it is considered ''
            status = parameters_dict[key].get('status', '')
            #if status == 'REQUIRED' and key in in_parameters:
            if status == 'REQUIRED':
                value = in_parameters[key]
                if value == '' or value is None:
                    error_message += f'The parameter {key} is required but not provided\n'
        return error_message

    # ===================================================================
    # ======================= VALIDATION OF TYPES =======================
    # ===================================================================
    def validate_types(self, in_parameters, parameters_info):
        error_message = ''
        for key, value in in_parameters.items():
            if value != '' and value is not None:
                input_type = parameters_info[key].get('type', '')

                if input_type == 'REAL' and not is_real_number(value):
                    error_message += f'The parameter {key} needs to be a real number\n'
                elif input_type == 'INTEGER' and not is_integer(value):
                    error_message += f'The parameter {key} needs to be an integer number\n'
                elif input_type == 'LOGICAL' and not is_logical(value):
                    error_message += f'The parameter {key} needs to be .TRUE. or .FALSE.\n'

        return error_message

    # ===================================================================
    # ================== validation of CONTROL fields ===================
    # ===================================================================
    def validate_control(self):
        control_input = self.input_data.get('CONTROL', {})
        error_message = ''

        error_message += self.validate_required_fields(control_input, CONTROL_DICT)
        error_message += self.validate_types(control_input, CONTROL_DICT)
        
        # ========== VALIDATE lfcp ==========
        if (control_input.get('lfcp', '') == '.TRUE.'):
            if (control_input.get('calulation', '')!='relax') or (control_input.get('assume_isolated', '')!='esm') or (control_input.get('esm', '')!='bc2' and control_input.get('esm', '')!='bc3'):
                error_message+='So that \'lfcp\' is .TRUE. does not meet the necessary requirements, see the More Info section for more information\n'

        return error_message


    # ===================================================================
    # ================== validation of SYSTEM fields ===================
    # ===================================================================
    def validate_system(self):
        system_input = self.input_data.get('SYSTEM', {})
        error_message = ''

        error_message += self.validate_required_fields(system_input, SYSTEM_DICT)
        error_message += self.validate_types(system_input, SYSTEM_DICT)

        if error_message == '':
            # Special Parameters
            error_message += self.validate_celldm(system_input)
            error_message += self.validate_ecutwfc(system_input.get('ecutwfc', ''))
            error_message += self.validate_nat_ntyp(system_input.get('nat', ''), system_input.get('ntyp', ''))
            # The other parameters
            error_message += self.validate_system_other_parameters(system_input)

        return error_message

    def validate_celldm(self, system):
        error_message = ''
        
        a = system.get('A', '')
        b = system.get('B', '')
        c = system.get('C', '')
        cosAB = system.get('cosAB', '')
        cosAC = system.get('cosAC', '')
        cosBC = system.get('cosBC', '')
        celldm_1 = system.get('celldm(1)', '')
        celldm_2 = system.get('celldm(2)', '')
        celldm_3 = system.get('celldm(3)', '')
        celldm_4 = system.get('celldm(4)', '')
        celldm_5 = system.get('celldm(5)', '')
        celldm_6 = system.get('celldm(6)', '')

        if (a!='' or b!='' or c!='' or cosAB!='' or cosAC!='' or cosBC!='') and (celldm_1!='' or celldm_2!='' or celldm_3!='' or celldm_4!='' or celldm_5!='' or celldm_6!=''):
            error_message += 'Specify either celldm(1...6) OR A,B,C,cosAB,cosBC,cosAC NOT both.\n'

        return error_message
    
    def validate_ecutwfc(self, ecutwfc):
        error_message = ''
        # Validate ecutwfc as a positive float and check if it's within a practical range
        try:
            ecutwfc = float(ecutwfc)
            if ecutwfc <= 0:
                error_message += ('ecutwfc must be a positive number.\n')
            # Check if ecutwfc is within a practical range
            if ecutwfc > 200:
                error_message += ('ecutwfc is unusually high. Please verify if this is intentional.\n')
        except ValueError:
            error_message += ('ecutwfc must be a number.\n')
        
        return error_message
    
    def validate_nat_ntyp(self, nat, ntyp):
        error_message = ''
        try:
            nat = int(nat)
            ntyp = int(ntyp)
            if nat <= 0 or ntyp <= 0:
                error_message += ('nat and ntyp must be positive integers.\n')
        except ValueError:
            error_message += ('nat and ntyp must be integers.\n')
        return error_message
    
    def validate_system_other_parameters(self, system):
        error_message = ''
        
        if (system.get('degauss_cond', '')!= '') and system.get('twochem', '').lower != '.true.':
            error_message += 'To use degauss_cond is necessary @twochem = .TRUE.\n'

        if system.get('gate', '').lower() != '.true.':
            # zgate, relaxz, block, block_1, block_2, block_height
            if (system.get('zgate', '')!='' or system.get('relaxz', '')!='' or system.get('block', '')!='' or system.get('block_1', '')!='' or system.get('block_2', '')!='' or system.get('block_height', '')!=''):
                error_message += 'To use zgate, relaxz, block, block_1, block_2 or block_height, is necessary use @gate = .TRUE.\n'

            # block, block_1, block_2, block_height
            if system.get('block', '').lower() != '.true.':
                if system.get('block_1', '')!='' or system.get('block_2', '')!='' or system.get('block_height', '')!='':
                    error_message += 'To use block_1, block_2 or block_height is necessary @gate = .TRUE. and @block = .TRUE.\n'
            else:
                if system.get('block_1', '')!='' or system.get('block_2', '')!='':
                    try:
                        if float(system.get('block_1', '')) < 0 and float(system.get('block_1', '')) > 1:
                            error_message += '@block_1 needs to be between 0 and 1\n'
                        if float(system.get('block_2', '')) < 0 and float(system.get('block_2', '')) > 1:
                            error_message += '@block_2 needs to be between 0 and 1\n'
                    except:
                        error_message += 'block_1 and block_2 must be float.\n'
            
            # one_atom_occupation
            if system.get('one_atom_occupations', '')!='':
                if system.get('nat', '')!='1' or system.get('occupations', '')!='from_input':
                    error_message += 'To use @one_atom_occupations is necessary @nat=1 and @occupations=\'from_input\'\n'
            
            # edir, emaxpos, eopreg
            if system.get('tefield', '').lower() != '.true.':
                if system.get('edir', '') != '':
                    error_message += 'To use edir is necessary @tefield = .TRUE.\n'
                if system.get('emaxpos', '') != '':
                    error_message += 'To use emaxpos is necessary @tefield = .TRUE.\n'
                if system.get('eopreg', '') != '':
                    error_message += 'To use eopreg is necessary @tefield = .TRUE.\n'
            else:
                if system.get('emaxpos', '')!='':
                    try:
                        if float(system.get('emaxpos', '')) < 0 and float(system.get('emaxpos', '')) > 1:
                                error_message += '@emaxpos needs to be between 0 and 1\n'   
                    except:
                        error_message += 'emaxpos must be float.\n'
                    if system.get('eopreg', '')!='':
                        try:
                            if float(system.get('eopreg', '')) < 0 and float(system.get('eopreg', '')) > 1:
                                    error_message += '@eopreg needs to be between 0 and 1\n'
                        except:
                            error_message += 'eopreg must be float.\n'

        return error_message
    
    #def validate_k_points(self, k_points):

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