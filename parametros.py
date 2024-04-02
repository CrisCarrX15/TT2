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
from PySide2.QtWidgets import QApplication, QMainWindow, QComboBox, QVBoxLayout, QPushButton, QScrollArea, QSizePolicy, QTabWidget, QWidget, QLineEdit, QLabel, QDialog
from PySide2.QtGui import QColor, QFont
from functools import partial
from PySide2.QtCore import Qt
from PySide2.QtWidgets import QTextEdit
from control_dict import CONTROL_DICT
from system_dict import SYSTEM_DICT
from quantum_espresso_io import create_in_file

class AppStyle:
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


class ParametrosWindow(QMainWindow):
    def __init__(self, control_dict, system_dict):
        super().__init__()
        self.control_dict = control_dict
        self.system_dict = system_dict
        self.setWindowTitle("Parametros")
        self.widget_dict = {}
        self.tab_widget_dict_widgets = {}
        self.tab_control_name = 'CONTROL'
        self.tab_system_name = 'SYSTEM'

        central_widget = QWidget()

        self.tab_widget = QTabWidget()
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
        widget_dict = {}
        for key, config in control_dict.items():
            full_text = f"{key}: {config.get('info', 'No information available')}"
            if len(full_text) > 50:
                label = QLabel(full_text[:50] + "...")
                button = QPushButton("Ver más")
                button.setStyleSheet("background-color: #CCCCCC;")  # Cambiar color a gris
                button.setMaximumWidth(100)  # Establecer el ancho máximo del botón
                self.connect_button_clicked(button, full_text)  # Conecta el botón al hacer clic
                layout.addWidget(label)
                layout.addWidget(button)
            else:
                label = QLabel(full_text)
                layout.addWidget(label)

            input_type = config.get('input_type', None)

            if input_type is not None and input_type == 'select_multiple':
                combobox = QComboBox()
                combobox.setFont(QFont("Arial", 12))
                combobox.addItems([''] + config['options'])
                layout.addWidget(combobox)
                widget_dict[label] = combobox
            else:
                text_edit = QLineEdit()
                layout.addWidget(text_edit)
                widget_dict[label] = text_edit

        self.tab_widget_dict_widgets[tab_index] = widget_dict

        button_layout = QVBoxLayout()
        button_layout.setSpacing(10)
        layout.addLayout(button_layout)

        save_button = QPushButton("Guardar")
        save_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        button_layout.addWidget(save_button)

        cancel_button = QPushButton("Cancelar")
        cancel_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        button_layout.addWidget(cancel_button)

        save_button.clicked.connect(self.guardar_informacion)

    def show_full_text(self, text):
        dialog = QDialog(self)
        dialog.setWindowTitle("Texto completo")
        dialog.setWindowModality(Qt.ApplicationModal)
        dialog.setMinimumWidth(400)  # Establecer el ancho mínimo del diálogo

        layout = QVBoxLayout()

        text_edit = QTextEdit()
        text_edit.setPlainText(text)
        text_edit.setReadOnly(True)  # Hacer el QTextEdit de solo lectura
        text_edit.setFont(QFont("Arial", 12))  # Establecer un tamaño de fuente más grande
        layout.addWidget(text_edit)

        close_button = QPushButton("Cerrar")
        close_button.clicked.connect(dialog.close)
        layout.addWidget(close_button)

        dialog.setLayout(layout)
        dialog.exec_()


    def guardar_informacion(self):
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
        widget_dict = self.tab_widget_dict_widgets.get(tab_index, {})
        for label, widget in widget_dict.items():
            if isinstance(widget, QComboBox):
                value = widget.currentText()
            elif isinstance(widget, QLineEdit):
                value = widget.text()
            else:
                value = None
            key = label.text().split(':')[0]
            info[key] = value
        return info

    def connect_button_clicked(self, button, text):
        def on_button_clicked():
            self.show_full_text(text)
        button.clicked.connect(on_button_clicked)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    parametros_window = ParametrosWindow(CONTROL_DICT, SYSTEM_DICT)
    AppStyle.apply(parametros_window)
    parametros_window.show()
    sys.exit(app.exec_())

