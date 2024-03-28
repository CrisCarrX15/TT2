import sys
from control_dict import CONTROL_DICT
from system_dict import SYSTEM_DICT
from PySide2.QtWidgets import QApplication, QMainWindow, QFormLayout, QLabel, QWidget, QComboBox, QScrollBar, QVBoxLayout, QPushButton, QScrollArea, QTabWidget
from PySide2.QtWidgets import QApplication, QMainWindow, QFormLayout, QLabel, QWidget, QComboBox, QScrollBar, QVBoxLayout, QPushButton, QScrollArea, QSizePolicy
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

        central_widget = QWidget()

        # Crear un widget de pestañas
        self.tab_widget = QTabWidget()

        # Crear pestañas para los parámetros
        self.tab1 = QWidget()
        self.tab2 = QWidget()

        self.tab_widget.addTab(self.tab1, "Control Dict")
        self.tab_widget.addTab(self.tab2, "System Dict")

        self.setup_tab(self.tab1, self.control_dict)
        self.setup_tab(self.tab2, self.system_dict)

        central_layout = QVBoxLayout(central_widget)
        central_layout.addWidget(self.tab_widget)

        self.setCentralWidget(central_widget)

    def setup_tab(self, tab, control_dict):
        scroll_area = QScrollArea()
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll_area.setWidgetResizable(True)

        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        scroll_layout.setSpacing(20)

        self.create_widgets(control_dict, scroll_layout)

        scroll_area.setWidget(scroll_widget)
        tab_layout = QVBoxLayout(tab)
        tab_layout.addWidget(scroll_area)

    
    def create_widgets(self, control_dict, layout):
        for key, config in control_dict.items():
            label = QLabel(f"{key}: {config.get('info', 'No information available')}")
            layout.addWidget(label)
            
            input_type = config.get('input_type', None)  # Get input_type if present, otherwise None

            if input_type is not None and input_type == 'select_multiple':
                combobox = QComboBox()
                combobox.setFont(QFont("Arial", 12))  # Establecer un tamaño de fuente más grande
                combobox.addItems([''] + config['options'])  # Agregar una opción vacía al principio
                layout.addWidget(combobox)
            else:
                text_edit = QLineEdit()  # Use un QLineEdit para que el usuario pueda escribir texto
                layout.addWidget(text_edit)

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

    def guardar_informacion(self):
        # Recopilar información de ambas pestañas y hacer algo con ella
        control_info = self.get_tab_info(self.tab1)
        system_info = self.get_tab_info(self.tab2)
        print("Información de Control Dict:", control_info)
        print("Información de System Dict:", system_info)

    def get_tab_info(self, tab):
        info = {}
        for widget in tab.findChildren(QComboBox):  # Buscar todos los ComboBox en la pestaña
            key = widget.parentWidget().layout().itemAt(0).widget().text().split(':')[0]  # Obtener la clave del QLabel
            value = widget.currentText()  # Obtener el texto seleccionado del ComboBox
            info[key] = value
        return info

if __name__ == "__main__":
    app = QApplication(sys.argv)
    parametros_window = ParametrosWindow(CONTROL_DICT, SYSTEM_DICT)
    AppStyle.apply(parametros_window)  # Aplicar el estilo a la ventana
    parametros_window.show()
    sys.exit(app.exec_())
