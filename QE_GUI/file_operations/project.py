import xml.etree.ElementTree as ET

def save_data(file_path, project_name, data):
    project = ET.Element(project_name)
    for section, parameters in data.items():
        section_element = ET.SubElement(project, section)
        for key, value in parameters.items():
            if str(value) != '' and value is not None:
                parameter_element = ET.SubElement(section_element, key)
                parameter_element.text = str(value)
    
    tree = ET.ElementTree(project)
    tree.write(f'{file_path}/{project_name}.qg')



def load_data(file):
    data = {}
    tree = ET.parse(file)
    root = tree.getroot()
    for section in root:
        parameters = {}
        for parameter in section:
            parameters[parameter.tag] = parameter.text
        data[section.tag] = parameters
    #print('data:')
    #print(data)
    return data

#loaded_data = load_data_from_xml(xml_file)

#print(loaded_data)
