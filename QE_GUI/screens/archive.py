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
from PySide2.QtWidgets import QApplication, QMainWindow, QPushButton, QPlainTextEdit, QVBoxLayout, QWidget, QLabel, QLineEdit, QGridLayout, QFileDialog, QDialog, QTextEdit, QSizePolicy, QMessageBox, QHBoxLayout
from PySide2.QtGui import QPixmap, QColor, QFont, QIcon
from PySide2.QtCore import Qt, Signal, QSize
import os
import subprocess
import matplotlib.pyplot as plt
import file_operations.quantum_espresso_io as qe_io
from file_operations.run_quantum_espresso import RunQuantumEspresso

class ArchiveWindow(QMainWindow):
    message_signal = Signal(str)

    def __init__(self):
        super().__init__()

        # Configurar la ventana principal
        self.setGeometry(100, 100, 1200, 800)
        self.setWindowTitle('Iteration Files')

        # Crear un layout vertical
        layout = QVBoxLayout()

        # Agregar una etiqueta con una imagen en la esquina
        imagen_path = './screens/images/GUI.jpeg'
        abs_imagen_path = os.path.abspath(imagen_path)
        print(f"Ruta absoluta de la imagen: {abs_imagen_path}")  # Imprimir la ruta absoluta
        imagen_label = QLabel(self)

        # Verificar que la imagen se cargue correctamente
        pixmap = QPixmap(abs_imagen_path)
        if not pixmap.isNull():
            # Ajustar el tamaño de la imagen a más pequeña
            pixmap = pixmap.scaled(100, 100, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            imagen_label.setPixmap(pixmap)
        else:
            print(f"Error: No se pudo cargar la imagen desde {abs_imagen_path}.")

        imagen_label.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        layout.addWidget(imagen_label)

        # Crear un widget para los controles
        controls_widget = QWidget()
        controls_layout = QVBoxLayout(controls_widget)

        # Crear la fila para los puntos K
        self.k_points_label = QLabel("Puntos K:")
        controls_layout.addWidget(self.k_points_label)

        self.grid_layout = QGridLayout()
        self.k_points_editors = [QLineEdit() for _ in range(6)]  # Fila de 6 QLineEdit
        for i in range(6):
            self.k_points_editors[i].setStyleSheet("background-color: white")  # Establecer fondo blanco
            self.grid_layout.addWidget(self.k_points_editors[i], 0, i)
        controls_layout.addLayout(self.grid_layout)

        # Crear el campo para la energía de corte
        self.energy_label = QLabel("Energía de Corte:")
        controls_layout.addWidget(self.energy_label)
        self.energy_editor = QLineEdit()
        self.energy_editor.setStyleSheet("background-color: white")  # Establecer fondo blanco
        controls_layout.addWidget(self.energy_editor)

        # Crear el campo para el número de iteraciones
        self.iterations_label = QLabel("Número de Iteraciones:")
        controls_layout.addWidget(self.iterations_label)
        self.iterations_editor = QLineEdit()
        self.iterations_editor.setStyleSheet("background-color: white")  # Establecer fondo blanco
        controls_layout.addWidget(self.iterations_editor)

        # Agregar botones para cargar archivos .in con imágenes y las opciones adicionales
        file_buttons_layout = QHBoxLayout()

        self.create_example_button(file_buttons_layout, "./screens/images/grafeno_1.png", "./examples/example_1.in")
        self.create_example_button(file_buttons_layout, "./screens/images/grafeno_2.png", "./examples/example_2.in")
        self.create_example_button(file_buttons_layout, "./screens/images/grafeno_3.png", "./examples/example_3.in")

        controls_layout.addLayout(file_buttons_layout)

        layout.addWidget(controls_widget)

        self.text = QPlainTextEdit()
        self.text.setReadOnly(True)
        layout.addWidget(self.text)

        # Crear un widget central y establecer el layout
        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

        # Establecer el color de fondo
        color_fondo = QColor(135, 158, 197)
        central_widget.setStyleSheet(f'background-color: {color_fondo.name()};')

    def message(self, s):
        self.text.appendPlainText(s)

    def create_example_button(self, layout, image_path, input_file):
        # Botón de imagen para ver el archivo .in
        btn_image = QPushButton()
        btn_image.setIcon(QIcon(image_path))
        btn_image.setIconSize(QSize(200, 200))
        btn_image.setStyleSheet("border: none;")
        btn_image.clicked.connect(lambda: self.show_file_content(input_file))
        layout.addWidget(btn_image)

        # Botón para ejecutar Quantum ESPRESSO una sola vez
        btn_run_qe = QPushButton("Run Quantum ESPRESSO")
        btn_run_qe.clicked.connect(lambda: self.run_single_qe(input_file))
        layout.addWidget(btn_run_qe)

        # Botón para iterar Quantum ESPRESSO
        btn_iterate_qe = QPushButton("Iterate Quantum ESPRESSO")
        btn_iterate_qe.clicked.connect(lambda: self.iterate_qe(input_file))
        layout.addWidget(btn_iterate_qe)

    def show_file_content(self, input_file):
        with open(input_file, 'r') as file:
            content = file.read()
        dialog = FileContentDialog(self, content)
        dialog.exec_()

    def run_single_qe(self, input_file):
        if not self.validate_inputs():
            self.message("Por favor, rellene todos los campos de puntos K y energía de corte.")
            return
        output_file = input_file.replace('.in', '.out')
        self.run_qe(input_file, output_file)

    def iterate_qe(self, input_file):
        if not self.validate_inputs():
            self.message("Por favor, rellene todos los campos de puntos K, energía de corte y número de iteraciones.")
            return

        # Obtener los valores de los puntos K y la energía de corte
        k_points = self.get_k_points()
        energy = float(self.energy_editor.text())
        iterations = self.iterations_editor.text()

        if not iterations.isdigit():
            self.message("Por favor, ingrese un número válido para las iteraciones.")
            return

        iterations = int(iterations)

        # Modificar el archivo de entrada con los nuevos valores
        qe_io.modify_k_points(input_file, k_points)
        qe_io.modify_ecut(input_file, energy)

        # Configurar los parámetros de iteración
        energy_diff_threshold = 0.001
        ecutwfc_increment = 5.0
        max_iterations = iterations

        # Variables para el proceso iterativo
        converged = False
        previous_energy = None
        iteration = 0
        energy_list = []
        ecutwfc_list = []

        run = RunQuantumEspresso()
        output_file = input_file.replace('.in', '.out')

        while not converged and iteration < max_iterations:
            # Ejecutar Quantum ESPRESSO
            run.run_qe_process(input_file, output_file)

            # Leer la energía total del sistema
            total_energy_str = qe_io.find_total_energy(output_file)
            total_energy = float(total_energy_str.split()[0])  # Asume que el formato es algo como "total_energy Ry"
            energy_list.append(total_energy)
            ecutwfc_list.append(energy)

            # Comprobar convergencia
            if previous_energy is not None:
                energy_diff = abs(total_energy - previous_energy)
                if energy_diff < energy_diff_threshold:
                    converged = True
                    self.message(f'Convergencia alcanzada en iteración {iteration+1} con ecutwfc={energy} y energía total={total_energy}')
                else:
                    # Ajustar ecutwfc
                    qe_io.sum_ecut(input_file, ecutwfc_increment)
            else:
                # Si es la primera iteración, simplemente ajustar ecutwfc
                qe_io.sum_ecut(input_file, ecutwfc_increment)

            # Actualizar la energía anterior
            previous_energy = total_energy
            energy += ecutwfc_increment  # Incrementar la energía de corte para la siguiente iteración
            iteration += 1

        if not converged:
            self.message(f'No se alcanzó la convergencia después de {max_iterations} iteraciones. Última energía total={total_energy}')

        # Generar la gráfica al final de la iteración
        self.plot_energies(ecutwfc_list, energy_list)

    def validate_inputs(self):
        return all(editor.text().strip() for editor in self.k_points_editors) and self.energy_editor.text().strip() and self.iterations_editor.text().strip()

    def get_k_points(self):
        k_points = []
        for i in range(0, 6, 3):
            k_points.append([int(self.k_points_editors[i].text()), int(self.k_points_editors[i+1].text()), int(self.k_points_editors[i+2].text())])
        return k_points

    def run_qe(self, input_file, output_file):
        run = RunQuantumEspresso()
        run.run_qe_process(input_file, output_file)
        self.message(f"Proceso QE ejecutado. Archivo de salida: {output_file}")

    def plot_energies(self, ecutwfc_list, energy_list):
        plt.figure(figsize=(10, 6))
        plt.plot(ecutwfc_list, energy_list, marker='o')
        plt.xlabel('Cut Energy (ecutwfc)')
        plt.ylabel('Total Energy (Ry)')
        plt.title('Total Energy vs Cut Energy')
        plt.grid(True)
        plt.show()

class FileContentDialog(QDialog):
    def __init__(self, parent, content):
        super(FileContentDialog, self).__init__(parent)
        self.setWindowTitle("File Content")
        self.resize(600, 400)  # Ajusta el tamaño de la ventana

        layout = QVBoxLayout()
        text_edit = QTextEdit()
        text_edit.setReadOnly(True)
        text_edit.setPlainText(content)
        layout.addWidget(text_edit)
        self.setLayout(layout)

def run_archive():
    # Verificar si ya existe una instancia de QApplication
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    archive_window = ArchiveWindow()
    archive_window.show()
    return app

if __name__ == "__main__":
    app = run_archive()
    sys.exit(app.exec_())


 



 
