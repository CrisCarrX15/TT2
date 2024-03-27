import sys
import control_dict as control
CONTROL_DICT = control.CONTROL_DICT
from PySide2.QtWidgets import QApplication, QMainWindow, QFormLayout, QLabel, QWidget, QComboBox, QScrollBar, QVBoxLayout, QPushButton, QScrollArea
from PySide2.QtGui import QColor, QFont
from PySide2.QtCore import Qt

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
            }
            QPushButton:hover {
                background-color: #e69500; /* Tonos más oscuros de naranja */
                border: 2px solid #e69500;
            }
            """
        # Combinar el estilo general de la ventana con el estilo de los botones
        window.setStyleSheet(window.styleSheet() + button_style)

class ParametrosWindow(QMainWindow):
    def __init__(self, control_dict):
        super().__init__()
        self.control_dict = control_dict
        self.setWindowTitle("Parametros")

        central_widget = QWidget()
        scroll_area = QScrollArea()
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll_area.setWidgetResizable(True)

        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        scroll_layout.setSpacing(20)

        self.create_widgets(scroll_layout)

        scroll_area.setWidget(scroll_widget)
        central_layout = QVBoxLayout(central_widget)
        central_layout.addWidget(scroll_area)

        self.setCentralWidget(central_widget)

    def create_widgets(self, layout):
        for key, config in self.control_dict.items():
            label = QLabel(f"{key}: {config['info']}")
            layout.addWidget(label)
            
            input_type = config.get('input_type', None)  # Get input_type if present, otherwise None

            if input_type is not None and input_type == 'select_multiple':
                combobox = QComboBox()
                combobox.setFont(QFont("Arial", 12))  # Establecer un tamaño de fuente más grande
                combobox.addItems([''] + config['options'])  # Agregar una opción vacía al principio
                layout.addWidget(combobox)
        
        # Agregar botones
        button_layout = QVBoxLayout()
        button_layout.setSpacing(10)
        layout.addLayout(button_layout)

        save_button = QPushButton("Guardar")
        button_layout.addWidget(save_button)

        cancel_button = QPushButton("Cancelar")
        button_layout.addWidget(cancel_button)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    parametros_window = ParametrosWindow(CONTROL_DICT)
    AppStyle.apply(parametros_window)  # Aplicar el estilo a la ventana
    parametros_window.show()
    sys.exit(app.exec_())

