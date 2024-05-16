import re
import xml.etree.ElementTree as ET

def save_data(file_path, project_name, data):
    project = ET.Element(project_name)
    for section, parameters in data.items():
        section_element = ET.SubElement(project, section)
        for key, value in parameters.items():
            if str(value) != '' and value is not None:
                key = re.sub(r'celldm\((\d+)\)', r'celldm_\1', key)
                parameter_element = ET.SubElement(section_element, key)
                parameter_element.text = str(value)
    
    tree = ET.ElementTree(project)
    tree.write(f'{file_path}/{project_name}.qg')



def load_data(file):
    data = {}
    try:
        with open(file, 'r', encoding='utf-8') as f:
            tree = ET.parse(f)
        root = tree.getroot()
        for section in root:
            parameters = {}
            for parameter in section:
                key = re.sub(r'celldm_(\d+)', r'celldm(\1)', parameter.tag)
                parameters[key] = parameter.text
            data[section.tag] = parameters
    except ET.ParseError as e:
        print(f"Error de análisis XML: {e}")
    return data


# Find the number of rows created for atomic_species and atomic_positions
def find_max_rows(file):
    tree = ET.parse(file)
    root = tree.getroot()
    rows = {}

    # Find the elements within ATOMIC_K_POINTS
    atomic_k_points = root.find("ATOMIC_K_POINTS")
    if atomic_k_points is not None:
        for child in atomic_k_points:
            # In the .qg file these values ​​are in the following format:
            # <atomic_{species}/{positions}_{row}_{column}>
            if child.tag.startswith("atomic_species") or child.tag.startswith("atomic_positions"):
                key = 'atomic_' + str(child.tag.split('_')[1])
                #print('Found ' + str(child.tag))
                row_number = int(child.tag.split("_")[2])
                if key not in rows or row_number+1 > rows.get(key, -1):
                    rows[key] = row_number+1
    #print('rows', rows)
    return rows
