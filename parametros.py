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


import sys
from control_dict import CONTROL_DICT
from system_dict import SYSTEM_DICT
from quantum_espresso_io import create_in_file
from PySide2.QtWidgets import QApplication, QMainWindow, QComboBox, QVBoxLayout, QPushButton, QScrollArea, QSizePolicy, QTabWidget, QWidget, QLineEdit, QLabel
from PySide2.QtGui import QColor, QFont
from PySide2.QtCore import Qt
from PySide2.QtWidgets import QLineEdit

class AppStyle:
    @staticmethod
    def apply(window):
        # Establecer el tamaño mínimo y máximo de la ventana
        window.setMinimumSize(800, 600)
        window.setMaximumSize(1600, 1200)

        # Establecer el color de fondo
        color_fondo = QColor("#5f7eb2")  # Color azul
        window.setStyleSheet(f'background-color: {color_fondo.name()};')

        # Estilo para los botones
        button_style = """
            QPushButton {
                background-color: #f0a500; /* Color naranja */
                color: white;
                font-size: 18px;
                min-width: 150px;
                min-height: 40px;
                border: 2px solid #f0a500; /* Mismo color que el fondo */
                border-radius: 5px;
                padding: 10px 20px; /* Aumentar el espacio interior */
            }
            QPushButton:hover {
                background-color: #e69500; /* Tonos más oscuros de naranja */
                border: 2px solid #e69500;
            }
            """
        # Combinar el estilo general de la ventana con el estilo de los botones
        window.setStyleSheet(window.styleSheet() + button_style)


class ParametrosWindow(QMainWindow):
    def __init__(self, control_dict, system_dict):
        super().__init__()
        self.control_dict = control_dict
        self.system_dict = system_dict
        self.setWindowTitle("Parametros")
        self.widget_dict = {}  # Inicializar el diccionario de widgets como variable de instancia
        self.tab_widget_dict_widgets = {}  # Diccionario para almacenar los widgets de cada pestaña
        # ADD NAME OF TABS
        self.tab_control_name = 'CONTROL'
        self.tab_system_name = 'SYSTEM'

        central_widget = QWidget()

        # Crear un widget de pestañas
        self.tab_widget = QTabWidget()

        # Crear pestañas para los parámetros
        self.tab1 = QWidget()
        self.tab2 = QWidget()

        self.tab_widget.addTab(self.tab1, self.tab_control_name)
        self.tab_widget.addTab(self.tab2, self.tab_system_name)

        self.setup_tab(self.tab1, self.control_dict, self.tab_control_name)
        self.setup_tab(self.tab2, self.system_dict, self.tab_system_name)

        central_layout = QVBoxLayout(central_widget)
        central_layout.addWidget(self.tab_widget)

        self.setCentralWidget(central_widget)

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

    
    def create_widgets(self, control_dict, layout, tab_index):
        widget_dict = {}  # Diccionario para almacenar los widgets de esta pestaña específica
        for key, config in control_dict.items():
            label = QLabel(f"{key}: {config.get('info', 'No information available')}")
            layout.addWidget(label)
            
            input_type = config.get('input_type', None)  # Get input_type if present, otherwise None

            if input_type is not None and input_type == 'select_multiple':
                combobox = QComboBox()
                combobox.setFont(QFont("Arial", 12))  # Establecer un tamaño de fuente más grande
                combobox.addItems([''] + config['options'])  # Agregar una opción vacía al principio
                layout.addWidget(combobox)
                widget_dict[label] = combobox  # Agregar el QLabel y el QComboBox al diccionario
            else:
                text_edit = QLineEdit()  # Use un QLineEdit para que el usuario pueda escribir texto
                layout.addWidget(text_edit)
                widget_dict[label] = text_edit  # Agregar el QLabel y el QLineEdit al diccionario
        
        self.tab_widget_dict_widgets[tab_index] = widget_dict  # Almacenar el diccionario de widgets para esta pestaña específica


        # Agregar botones
        button_layout = QVBoxLayout()
        button_layout.setSpacing(10)
        layout.addLayout(button_layout)

        save_button = QPushButton("Guardar")
        save_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)  # Establecer política de tamaño para expandirse horizontalmente
        button_layout.addWidget(save_button)

        cancel_button = QPushButton("Cancelar")
        cancel_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)  # Establecer política de tamaño para expandirse horizontalmente
        button_layout.addWidget(cancel_button)

        # Conectar el botón "Guardar" a la función guardar_informacion
        save_button.clicked.connect(self.guardar_informacion)

    def guardar_informacion(self):
     # Recopilar información de ambas pestañas y hacerself.tab1, s algo con ella
        control_info = self.get_tab_info(self.tab_control_name)
        system_info = self.get_tab_info(self.tab_system_name)

        info_dict = {}
        info_dict['CONTROL'] = control_info
        info_dict['SYSTEM'] = system_info

        create_in_file(info_dict)
    
        print("Información de Control Dict:")
        for key, value in control_info.items():
            print(f"{key}: {value}")

        print("\nInformación de System Dict:")
        for key, value in system_info.items():
            print(f"{key}: {value}")

    def get_tab_info(self, tab_index):
        info = {}
        widget_dict = self.tab_widget_dict_widgets.get(tab_index, {})  # Obtener el diccionario de widgets para la pestaña actual
        for label, widget in widget_dict.items():
            if isinstance(widget, QComboBox):  # Verificar si el widget es un QComboBox
                value = widget.currentText()  # Obtener el texto seleccionado del QComboBox
            elif isinstance(widget, QLineEdit):  # Verificar si el widget es un QLineEdit
                value = widget.text()  # Obtener el texto escrito en el QLineEdit
            else:
                value = None
            key = label.text().split(':')[0]  # Obtener la clave del QLabel
            print(f'Encontrando en el parametro {key} el valor {value}')
            info[key] = value
        return info

if __name__ == "__main__":
    app = QApplication(sys.argv)
    parametros_window = ParametrosWindow(CONTROL_DICT, SYSTEM_DICT)
    AppStyle.apply(parametros_window)  # Aplicar el estilo a la ventana
    parametros_window.show()
    sys.exit(app.exec_())
