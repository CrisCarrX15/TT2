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

from PySide2.QtWidgets import QHBoxLayout, QMainWindow, QComboBox, QVBoxLayout, QPushButton, QScrollArea, QSizePolicy, QTabWidget, QWidget, QLineEdit, QLabel, QDialog,QGridLayout
from PySide2.QtGui import QColor, QFont
from functools import partial
from PySide2.QtCore import Qt
from PySide2.QtWidgets import QTextEdit
from parameters.control_dict import CONTROL_DICT
from parameters.system_dict import SYSTEM_DICT
from parameters.electrons_dict import ELECTRONS_DICT
from parameters.ions_dict import IONS_DICT
from parameters.atomic_dict import INPUT_DATA
from parameters.config_dict import CONFIG_DICT
from parameters.rism_dict import RISM_DICT
from file_operations.quantum_espresso_io import create_in_file
from file_operations.project import save_data, load_data

class AppParametersStyle:
    @staticmethod
    def apply(window):
        window.setMinimumSize(800, 600)
        window.setMaximumSize(1600, 1200)
        color_fondo = QColor("#5f7eb2")
        window.setStyleSheet(f'background-color: {color_fondo.name()};')

        button_style = """
            QPushButton {
                background-color: #f0a500;
                color: white;
                font-size: 18px;
                min-width: 150px;
                min-height: 40px;
                border: 2px solid #f0a500;
                border-radius: 5px;
                padding: 10px 20px;
            }
            QPushButton:hover {
                background-color: #e69500;
                border: 2px solid #e69500;
            }
            """
        window.setStyleSheet(window.styleSheet() + button_style)


class ParametersWindow(QMainWindow):
    def __init__(self, project_name, file_path, load_qg_project):
        super().__init__()
        self.control_dict = CONTROL_DICT
        self.system_dict = SYSTEM_DICT
        self.electrons_dict = ELECTRONS_DICT
        self.ions_dict = IONS_DICT 
        self.rism_dict = RISM_DICT
        self.atomic_dict = INPUT_DATA
        self.config_dict = CONFIG_DICT
        self.project_name = project_name
        self.file_path = file_path
        self.setWindowTitle(f'Project: {project_name}')
        self.widget_dict = {}
        self.tab_widget_dict_widgets = {}
        self.tab_control_name = 'CONTROL'
        self.tab_system_name = 'SYSTEM'
        self.tab_electrons_name = 'ELECTRONS'  
        self.tab_ions_name ='IONS'
        self.tab_rism_name = 'RISM'
        self.tab_input_name = 'ATOMIC & K_POINTS'
        self.tab_config_name = 'SETTINGS'

        self.atomic_positions_counter = 0
        self.atomic_species_counter = 0

        central_widget = QWidget()

        self.tab_widget = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tab3 = QWidget()  
        self.tab4 = QWidget()
        self.tab5 = QWidget()
        self.tab6 = QWidget()
        self.tab7 = QWidget()

        self.tab_widget.addTab(self.tab1, self.tab_control_name)
        self.tab_widget.addTab(self.tab2, self.tab_system_name)
        self.tab_widget.addTab(self.tab3, self.tab_electrons_name)  
        self.tab_widget.addTab(self.tab4, self.tab_ions_name)
        self.tab_widget.addTab(self.tab5, self.tab_rism_name)
        self.tab_widget.addTab(self.tab6, self.tab_input_name)
        self.tab_widget.addTab(self.tab7, self.tab_config_name)


        self.setup_tab(self.tab1, self.control_dict, self.tab_control_name)
        self.setup_tab(self.tab2, self.system_dict, self.tab_system_name)
        self.setup_tab(self.tab3, self.electrons_dict, self.tab_electrons_name)  
        self.setup_tab(self.tab4, self.ions_dict, self.tab_ions_name)
        self.setup_tab(self.tab5, self.rism_dict, self.tab_rism_name)
        self.setup_tab(self.tab6, self.atomic_dict, self.tab_input_name)
        self.setup_tab(self.tab7, self.config_dict, self.tab_config_name)

        central_layout = QVBoxLayout(central_widget)
        central_layout.addWidget(self.tab_widget)

        self.setCentralWidget(central_widget)

        if load_qg_project:
            data = load_data(f'{file_path}/{self.project_name}.qg')
            self.load_data_into_interface(data)


    def setup_tab(self, tab, control_dict, tab_index):
        scroll_area = QScrollArea()
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll_area.setWidgetResizable(True)

        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        scroll_layout.setSpacing(20)

        self.create_widgets(control_dict, scroll_layout, tab_index)

        scroll_area.setWidget(scroll_widget)
        tab_layout = QVBoxLayout(tab)
        tab_layout.addWidget(scroll_area)

# ==================== NUEVO =========================
    # ============== CREATE THE WIDGET FOR ATOMIC SPECIES ===============
    def create_widgets_atomic_species(self, config, layout):
        widget_dict = self.tab_widget_dict_widgets.get(self.tab_input_name, {})

        text_layout = QVBoxLayout() 
        text_matrix_layout = QHBoxLayout()

        for i in range(3):  # Crear 1x3 entradas de texto para atomic_species
            text_edit = QLineEdit()
            text_matrix_layout.addWidget(text_edit)
            identifier = f"atomic_species_0_{i}"
            widget_dict[identifier] = text_edit

        text_layout.addLayout(text_matrix_layout)

        add_table_button = QPushButton("Add row")
        add_table_button.clicked.connect(lambda: self.add_additional_table_atomic_species(text_layout))
        text_layout.addWidget(add_table_button)

        widget = QWidget()
        widget.setLayout(text_layout)
        layout.addWidget(widget)

        return widget_dict

    # =========== ADD NEW ROW OF ENTRIES FOR ATOMIC SPECIES =============
    def add_additional_table_atomic_species(self, text_layout):
        widget_dict = self.tab_widget_dict_widgets.get(self.tab_input_name, {})
        # Obtener el número actual de filas
        current_rows = text_layout.count() - 1  # Excluyendo el botón "Agregar fila"

        # Crear una nueva fila de entradas de texto
        text_matrix_layout = QHBoxLayout()
        for i in range(3):  # 3 columnas para atomic_species
            text_edit = QLineEdit()
            text_matrix_layout.addWidget(text_edit)
            # Construir el identificador único en el formato requerido
            identifier = f"atomic_species_{current_rows}_{i}"
            widget_dict[identifier] = text_edit

        # Agregar la nueva fila al layout
        text_layout.insertLayout(current_rows, text_matrix_layout)

        # Actualizar el diccionario de widgets de la pestaña actual
        self.tab_widget_dict_widgets[self.tab_input_name] = widget_dict

    # ============== CREATE THE WIDGET FOR ATOMIC POSITIONS ===============
    def create_widgets_atomic_positions(self, config, layout):
        widget_dict = self.tab_widget_dict_widgets.get(self.tab_input_name, {})

        text_layout = QVBoxLayout() 
        text_matrix_layout = QHBoxLayout()

        for i in range(4):  # Crear 1x4 entradas de texto para atomic_positions
            text_edit = QLineEdit()
            text_matrix_layout.addWidget(text_edit)
            identifier = f"atomic_positions_0_{i}"
            widget_dict[identifier] = text_edit

        text_layout.addLayout(text_matrix_layout)

        add_table_button = QPushButton("Add row")
        add_table_button.clicked.connect(lambda: self.add_additional_table_atomic_positions(text_layout))
        text_layout.addWidget(add_table_button)

        widget = QWidget()
        widget.setLayout(text_layout)
        layout.addWidget(widget)

        return widget_dict

    # =========== ADD NEW ROW OF ENTRIES FOR ATOMIC POSITIONS =============
    def add_additional_table_atomic_positions(self, text_layout):
        widget_dict = self.tab_widget_dict_widgets.get(self.tab_input_name, {})
        # Obtener el número actual de filas
        current_rows = text_layout.count() - 1  # Excluyendo el botón "Agregar fila"

        # Crear una nueva fila de entradas de texto
        text_matrix_layout = QHBoxLayout()
        for i in range(4):  # 4 columnas para atomic_positions
            text_edit = QLineEdit()
            text_matrix_layout.addWidget(text_edit)
            # Construir el identificador único en el formato requerido
            identifier = f"atomic_positions_{current_rows}_{i}"
            widget_dict[identifier] = text_edit

        # Agregar la nueva fila al layout
        text_layout.insertLayout(current_rows, text_matrix_layout)

        # Actualizar el diccionario de widgets de la pestaña actual
        self.tab_widget_dict_widgets[self.tab_input_name] = widget_dict

# ==================== NUEVO =========================

    def create_widgets(self, control_dict, layout, tab_index):
        widget_dict = {}
        for key, config in control_dict.items():
            description_text = f"{key}: {config.get('description', 'No description available')} ({config.get('type', 'Any type')})"
            label = QLabel(description_text)
            layout.addWidget(label)

            full_info = config.get('info', 'No information available')
        
            button = QPushButton("More")
            button.setStyleSheet("background-color: #CCCCCC;")  # Cambiar color a gris
            button.setMaximumWidth(100)  # Establecer el ancho máximo del botón
            self.connect_button_clicked(button, full_info, key)  # Conecta el botón al hacer clic
            layout.addWidget(button)
        

            input_type = config.get('input_type', None)

            if input_type is None:
                text_edit = QLineEdit()
                layout.addWidget(text_edit)
                widget_dict[key] = text_edit
            elif input_type == 'select_multiple':
                combobox = QComboBox()
                combobox.setFont(QFont("Arial", 12))
                combobox.addItems([''] + config['options'])
                layout.addWidget(combobox)
                widget_dict[key] = combobox
            elif input_type == 'matrix':
                matrix_layout = QGridLayout()
                for i in range(2):
                    for j in range(3):
                        text_edit = QLineEdit()
                        matrix_layout.addWidget(text_edit, i, j)
                        widget_dict[f'{key}_{i}_{j}'] = text_edit
                matrix_widget = QWidget()
                matrix_widget.setLayout(matrix_layout)
                layout.addWidget(matrix_widget)
            elif input_type == 'atomic_species':
                widget_dict.update(self.create_widgets_atomic_species(config, layout))
            elif input_type == 'atomic_positions':
                widget_dict.update(self.create_widgets_atomic_positions(config, layout))
            else:
                text_edit = QLineEdit()
                layout.addWidget(text_edit)
                widget_dict[key] = text_edit

        self.tab_widget_dict_widgets[tab_index] = widget_dict

        button_layout = QVBoxLayout()
        button_layout.setSpacing(10)
        layout.addLayout(button_layout)

        save_project_button = QPushButton(f'Save {self.project_name} project')
        save_project_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        save_project_button.clicked.connect(self.save_project)
        button_layout.addWidget(save_project_button)

        save_in_button = QPushButton(f'Save {self.project_name}.in file')
        save_in_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        save_in_button.clicked.connect(self.create_in)
        button_layout.addWidget(save_in_button)

        cancel_button = QPushButton("Cancel")
        cancel_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        button_layout.addWidget(cancel_button)

    
    def show_full_text(self, text, parameter):
        dialog = QDialog(self)
        dialog.setWindowTitle(parameter)
        dialog.setWindowModality(Qt.ApplicationModal)
        dialog.setMinimumWidth(400)  # Establecer el ancho mínimo del diálogo

        layout = QVBoxLayout()

        text_edit = QTextEdit()
        text_edit.setPlainText(text)
        text_edit.setReadOnly(True)  # Hacer el QTextEdit de solo lectura
        text_edit.setFont(QFont("Arial", 12))  # Establecer un tamaño de fuente más grande
        layout.addWidget(text_edit)

        close_button = QPushButton("Close")
        close_button.clicked.connect(dialog.close)
        layout.addWidget(close_button)

        dialog.setLayout(layout)
        dialog.exec_()

    def save_project(self):
        control_info = self.get_tab_info(self.tab_control_name)
        system_info = self.get_tab_info(self.tab_system_name)
        electrons_info = self.get_tab_info(self.tab_electrons_name)
        ions_info = self.get_tab_info(self.tab_ions_name)
        rism_info = self.get_tab_info(self.tab_rism_name)
        input_info = self.get_tab_info(self.tab_input_name)
        input_info2 = self.get_tab_info(self.tab_input_name)
        input_info3 = self.get_tab_info(self.tab_input_name)
        config_info = self.get_tab_info(self.tab_config_name)


        info_dict = {}
        info_dict['CONTROL'] = control_info
        info_dict['SYSTEM'] = system_info
        info_dict['ELECTRONS'] = electrons_info
        info_dict['IONS'] = ions_info
        info_dict['RISM'] = rism_info
        info_dict['ATOMIC_SPECIES']= input_info
        info_dict['ATOMIC_POSISITIONS']= input_info2
        info_dict['K_POINTS']= input_info3
        info_dict['SETTINGS'] = config_info

        save_data(self.file_path, self.project_name, info_dict)

    def create_in(self):
        control_info = self.get_tab_info(self.tab_control_name)
        system_info = self.get_tab_info(self.tab_system_name)
        electrons_info = self.get_tab_info(self.tab_electrons_name)
        ions_info = self.get_tab_info(self.tab_ions_name)
        rism_info = self.get_tab_info(self.tab_rism_name)
        input_info = self.get_tab_info(self.tab_input_name)
        config_info = self.get_tab_info(self.tab_config_name)

        info_dict = {}
        info_dict['CONTROL'] = control_info
        info_dict['SYSTEM'] = system_info
        info_dict['ELECTRONS'] = electrons_info
        info_dict['IONS'] = ions_info
        info_dict['RISM'] = rism_info
        info_dict['ATOMIC'] = input_info
        info_dict['SETTINGS'] = config_info

        #create_in_file(self.file_path, self.project_name,info_dict)
    
        print("Información de Control Dict:")
        for key, value in control_info.items():
            print(f"{key}: {value}")

        print("\nInformación de System Dict:")
        for key, value in system_info.items():
            print(f"{key}: {value}")
        
        print("\nInformación de Electrons Dict:")
        for key, value in electrons_info.items():
            print(f"{key}: {value}")

        print("\nInformación de Ions Dict:")
        for key, value in ions_info.items():
            print(f"{key}: {value}")
        
        print("\nInformación de Rims Dict:")
        for key, value in rism_info.items():
            print(f"{key}: {value}")

        print("\nInformación de Input Data:")
        for key, value in input_info.items():
            print(f"{key}: {value}")
        
        print("\nInformación de Config Dict:")
        for key, value in config_info.items():
            print(f"{key}: {value}")
        
    
    def get_tab_info(self, tab_index):
        info = {}
        widget_dict = self.tab_widget_dict_widgets.get(tab_index, {})
        for label, widget in widget_dict.items():
            if isinstance(widget, QComboBox):
                value = widget.currentText()
            elif isinstance(widget, QLineEdit):
                value = widget.text()
            else:
                value = None

            #key = label.text().split(':')[0].strip()
            key = label
            info[key] = value
        return info

    def connect_button_clicked(self, button, text, parameter):
        def on_button_clicked():
            self.show_full_text(text, parameter)
        button.clicked.connect(on_button_clicked)
    

    def load_data_into_interface(self, data):
        for section, parameters in data.items():
            widget_dict = self.tab_widget_dict_widgets.get(section)
            if widget_dict:
                for parameter, value in parameters.items():
                    for label, widget in widget_dict.items():
                        # Obtener el texto de la etiqueta del widget y eliminar cualquier texto adicional
                        #label_text = label.text().split(':')[0].strip()
                        label_text = label
                        if label_text == parameter:
                            if isinstance(widget, QComboBox):
                                index = widget.findText(value)
                                if index != -1:
                                    widget.setCurrentIndex(index)
                            elif isinstance(widget, QLineEdit):
                                widget.setText(value)



def run_parameters(project_name, file_path, load_qg_project):
    parameters_window = ParametersWindow(project_name, file_path, load_qg_project)
    AppParametersStyle.apply(parameters_window)
    parameters_window.show()
"""
if __name__ == "__main__":
    app = QApplication(sys.argv)
    parametros_window = ParametersWIndow(CONTROL_DICT, SYSTEM_DICT)
    AppStyle.apply(parametros_window)
    parametros_window.show()
    sys.exit(app.exec_())
"""
