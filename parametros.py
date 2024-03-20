import sys
from PySide2.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QComboBox, QLineEdit, QGroupBox, QCheckBox, QScrollArea
from PySide2.QtCore import Qt

import input_values

class Ventana(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Selección de parámetros")
        self.setMinimumSize(600, 400)  # Establecer el tamaño mínimo de la ventana
        self.setMaximumSize(800, 600)  # Establecer el tamaño máximo de la ventana

        # Crear un área de desplazamiento
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)  # Permitir que el contenido sea redimensionable
        self.layout = QVBoxLayout()

        # Iterar sobre los parámetros definidos en CONTROL_VALUES
        for parametro, detalles in input_values.CONTROL_VALUES.items():
            label = QLabel(parametro)
            self.layout.addWidget(label)

            # Dependiendo del input_type, crea el widget correspondiente
            if detalles['input_type'] == 'select_multiple':
                combo_box = QComboBox()
                for valor in detalles['values']:
                    combo_box.addItem(valor)
                combo_box.currentIndexChanged.connect(self.on_selection_change)
                self.layout.addWidget(combo_box)

            elif detalles['input_type'] == 'text':
                line_edit = QLineEdit()
                line_edit.textChanged.connect(self.on_text_change)
                self.layout.addWidget(line_edit)

            elif detalles['input_type'] == 'automatic':
                line_edit = QLineEdit()
                line_edit.setReadOnly(True)  # Hacer que el QLineEdit sea de solo lectura
                self.layout.addWidget(line_edit)

            else:
                # Manejar otros tipos de entrada si es necesario
                pass

        # Crear un widget que contendrá el diseño principal
        contenido_widget = QWidget()
        contenido_widget.setLayout(self.layout)
        self.scroll_area.setWidget(contenido_widget)

        # Establecer el diseño principal de la ventana como el área de desplazamiento
        self.main_layout = QVBoxLayout()
        self.main_layout.addWidget(self.scroll_area)
        self.setLayout(self.main_layout)

    def on_selection_change(self, index):
        combo_box = self.sender()
        selected_value = combo_box.currentText()
        print(f"Opción seleccionada: {selected_value}")

    def on_text_change(self, text):
        line_edit = self.sender()
        print(f"Texto escrito: {text}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    ventana = Ventana()
    ventana.show()
    sys.exit(app.exec_())
