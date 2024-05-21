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

from file_operations.input_validation import is_real_number, is_integer

# ========== FIND TOTAL ENERGY ==========
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


# ========== MODIFY K_POINTS ==========
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


# ========== MODIFY ECUTWFC ==========
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


# ========== ADD CERTAIN AMOUNT TO ECUTWFC ==========
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


# ========== EXTRACT ATOMIC_POSITIONS IN .IN FILE ==========
def extract_atomic_positions_in(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    atomic_positions = []
    in_atomic_positions = False

    for line in lines:
        if 'ATOMIC_POSITIONS' in line:
            in_atomic_positions = True
            continue  # Saltar la línea de cabecera

        if in_atomic_positions:
            if line.strip() == '' or 'K_POINTS' in line or 'CELL_PARAMETERS' in line:
                break  # Salir si se encuentra una línea vacía
            atomic_positions.append(line.strip())

    return atomic_positions


# ========== EXTRACT ATOMIC_POSITIONS IN .OUT FILE ==========
def extract_atomic_positions_out(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    atomic_positions = []
    in_atomic_positions = False

    for line in lines:
        if 'positions (alat units)' in line:
            in_atomic_positions = True
            continue  # Saltar la línea de cabecera

        if in_atomic_positions:
            #print(f'Analizando la linea {line}')
            if 'tau(' in line:
                parts = line.split()
                atom = parts[1]
                x = parts[-4].strip('()')
                y = parts[-3].strip('()')
                z = parts[-2].strip('()')
                #print(f"{atom} {x} {y} {z}")
                atomic_positions.append(f"{atom} {x} {y} {z}")
            else:
                break  # Salir si no se encuentra 'tau('

    return atomic_positions

# ========== CREATE .XYZ FILE ==========
def create_xyz_file(atomic_positions, output_filename, description="Molecula"):
    num_atoms = len(atomic_positions)
    with open(output_filename, 'w') as file:
        file.write(f"{num_atoms}\n")
        file.write(f"{description}\n")
        for pos in atomic_positions:
            file.write(f"{pos}\n")



def create_in_file(file_path, filename, parameters):
    filename = f'{file_path}/{filename}.in'
    with open(filename, 'w') as f:
        for section, params in parameters.items():
            if section == 'ATOMIC_K_POINTS':

                # Verificación para ATOMIC_SPECIES
                write_atomic_species = any(param != '' for key, param in params.items() if key.startswith('atomic_species'))

                # Verificación para ATOMIC_POSITIONS
                write_atomic_positions = any(param != '' for key, param in params.items() if key.startswith('atomic_positions'))

                write_k_points = all(param != '' for param in params.values() if param.startswith('K_POINTS'))

                # Verificación para CELL_PARAMETERS
                write_cell_parameters = any(param != '' for key, param in params.items() if key.startswith('CELL_PARAMETERS'))

                if write_atomic_species:
                    max_row_atomic_species = max(int(param.split('_')[2]) for param in params if param.startswith('atomic_species'))
                    f.write("ATOMIC_SPECIES\n")
                    for i in range(max_row_atomic_species + 1):
                        f.write(f'    {params.get(f"atomic_species_{i}_0", "")}  {float(params.get(f"atomic_species_{i}_1", 0)):.4f}  {params.get(f"atomic_species_{i}_2", "")}\n')
                    f.write("\n\n")

                if write_atomic_positions:
                    max_row_atomic_positions = max(int(param.split('_')[2]) for param in params if param.startswith('atomic_positions'))
                    f.write("ATOMIC_POSITIONS alat\n")
                    for i in range(max_row_atomic_positions + 1):
                        f.write(f'    {params.get(f"atomic_positions_{i}_0", "")}  {float(params.get(f"atomic_positions_{i}_1", 0)):.6f}  {float(params.get(f"atomic_positions_{i}_2", 0)):.6f}  {float(params.get(f"atomic_positions_{i}_3", 0)):.6f}\n')
                    f.write("\n\n")

                if write_k_points:
                    f.write("K_POINTS automatic\n  ")
                    for i in range(6):  # Siempre hay 6 parámetros de K_POINTS
                        f.write(f'{params.get(f"K_POINTS_{i}", "")}')
                        if i < 2:
                            f.write(' ')
                        elif i == 2:
                            f.write('   ')
                        elif i < 5:
                            f.write(' ')
                    f.write("\n\n")
                
                if write_cell_parameters:
                    f.write("CELL_PARAMETERS\n")
                    for i in range(3):
                        for j in range(3):
                            f.write(f'    {float(params.get(f"CELL_PARAMETERS_{i}_{j}", 0)):.6f}')
                        f.write('\n')
                    f.write('\n')
            else:
                if any(value != '' and value is not None for value in params.values()):
                    f.write(f"&{section}\n")
                    for param, value in params.items():
                        if value != '' and value is not None:
                            if is_real_number(value) or is_integer(value):
                                f.write(f"\t{param} = {value},\n")
                            else:
                                f.write(f"\t{param} = \'{value}\',\n")
                    f.write("/\n\n")

# ========== FIND ERROR OR JOB DONE IN .OUT FILE ==========
def check_qe_output(file_path):
    message = ''
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            
        # Check for specific error messages
        error_messages = []
        error_detected = False
        job_done = False
        
        for line in lines:
            if "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" in line:
                error_detected = True
                continue
            if error_detected and "Error in routine" in line:
                error_messages.append(line.strip())
                error_detected = False
            
            if "JOB DONE" in line:
                job_done = True
        
        if job_done:
            message += ("JOB DONE")
        elif error_messages:
            message += ("Errors detected in the output file:")
            for msg in error_messages:
                message += (f" - {msg}")
        else:
            message += ("The job did not complete successfully, but no specific errors were found.")
    
    except FileNotFoundError:
        message += (f"File {file_path} not found.")
    except Exception as e:
        message += (f"An error occurred while reading the file: {e}")
    
    return message

