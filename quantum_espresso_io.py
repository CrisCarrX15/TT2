###################################################
##                                               ##
##  quantum_espresso_io.py                       ##
##    Modules to create, modify and delete       ##
##    Quantum ESPRESSO files                     ##
##                                               ##
##  Authors:                                     ##
##    Marco Uriel Aguilar Lara                   ##
##    Cristian Eduardo Carrillo Soto             ##
##                                               ##
##  From: ESCOM, National Polytechnic Institute  ##
##                                               ##
###################################################

from input_validation import is_real_number, is_integer

def find_total_energy(filename):
    # Read the input file
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    total_energy = "0 Ry"
    # Find the line that contains the total energy
    for i, line in enumerate(lines):
        if '!    total energy' in line:
            # Find the position of = and extract the current value
            pos_equal = line.find('=')
            #line = line.replace('Ry','')
            total_energy = line[pos_equal+1:].strip()
            break  # Breaks the loop once it has been found
    return(total_energy)


def modify_k_points(filename, k_points):
    # Read the input file
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    print('Escribiendo puntos k...')
    #print('Los k_points son: ', k_points)
    for i, line in enumerate(lines):
        if 'K_POINTS automatic' in line:
            k_points_str = '  '
            for k_point in k_points:
                k_points_str += str(k_point[0]) + ' ' + str(k_point[1]) + ' ' + str(k_point[2]) + '   '
            lines[i+1] = k_points_str

    # Save the modified file
    with open(filename, 'w') as modified_file:
        modified_file.writelines(lines)


def modify_ecut(filename, ecut):
    # Read the input file
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Find the line that contains the cutting energy
    for i, line in enumerate(lines):
        if 'ecutwfc' in line:
            # Actualiza la línea con el nuevo valor
            lines[i] = f' ecutwfc = {float(ecut)},\n'
            break  # Breaks the loop once the line has been found and modified

    # Save the modified file
    with open(filename, 'w') as modified_file:
        modified_file.writelines(lines)   


def sum_ecut(filename, number_to_add):
    # Read the input file
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Find the line that contains the cutting energy
    for i, line in enumerate(lines):
        if 'ecutwfc' in line:
            # Find the position of = and extract the current value
            pos_equal = line.find('=')
            line = line.replace(',','')
            actual_value = float(line[pos_equal+1:].strip())
            
            # Add the desired amount
            new_value = actual_value + number_to_add

            # Update the line with the new value
            lines[i] = f'  ecutwfc = {new_value},\n'
            break  # Breaks the loop once the line has been found and modified

    # Save the modified file
    with open(filename, 'w') as modified_file:
        modified_file.writelines(lines)


def create_in_file(parameters):
    filename = 'test_in_file.in'
    with open(filename, 'w') as f:
        for section, params in parameters.items():
            f.write(f"&{section}\n")
            for param, value in params.items():
                if value != '' and value is not None:
                    if is_real_number(value) or is_integer(value):
                        f.write(f"\t{param} = {value}\n")
                    else:
                        f.write(f"\t{param} = \'{value}\'\n")
            f.write("\\\n\n")
