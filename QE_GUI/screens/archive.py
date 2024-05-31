###################################################
##                                               ##
##    archive.py                                 ##
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
from PySide2.QtCore import Qt, Signal, QSize, QThread, QTimer
import os
import subprocess
import threading
import traceback
import matplotlib.pyplot as plt
import file_operations.quantum_espresso_io as qe_io
from file_operations.run_quantum_espresso import RunQuantumEspresso

class LoadingDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Loading")
        self.setWindowModality(Qt.ApplicationModal)  # Lock the entire app
        # Adjust window flags to remove close button
        self.setWindowFlags(self.windowFlags() & ~Qt.WindowCloseButtonHint)

        layout = QVBoxLayout()
        layout.addWidget(QLabel("Running Quantum ESPRESSO, please wait..."))
        self.setLayout(layout)
        self.setFixedSize(400, 100)

class QEIterateThread(QThread):
    progress_signal = Signal(str)
    finished_signal = Signal(str)

    def __init__(self, parent, input_file, k_points, energy, max_iterations, energy_diff_threshold=0.001, ecutwfc_increment=5.0):
        super().__init__(parent)
        self.input_file = input_file
        self.k_points = k_points
        self.energy = energy
        self.max_iterations = max_iterations
        self.energy_diff_threshold = energy_diff_threshold
        self.ecutwfc_increment = ecutwfc_increment
        self.energy_list = []
        self.ecutwfc_list = []

    def run(self):
        import file_operations.quantum_espresso_io as qe_io
        from file_operations.run_quantum_espresso import RunQuantumEspresso

        converged = False
        previous_energy = None
        iteration = 0
        run = RunQuantumEspresso()
        output_file = self.input_file.replace('.in', '.out')

        while not converged and iteration < self.max_iterations:
            try:
                run.run_qe_process(self.input_file, output_file)
                total_energy_str = qe_io.find_total_energy(output_file)
                total_energy = float(total_energy_str.split()[0])
                self.energy_list.append(total_energy)
                self.ecutwfc_list.append(self.energy)

                self.progress_signal.emit(f'Iteration {iteration+1} with ecutwfc={self.energy} and total energy={total_energy}')

                if previous_energy is not None:
                    energy_diff = abs(total_energy - previous_energy)
                    if energy_diff < self.energy_diff_threshold:
                        converged = True
                        self.progress_signal.emit(f'Convergence achieved in iteration {iteration+1} with ecutwfc={self.energy} and total energy={total_energy}')
                    else:
                        qe_io.sum_ecut(self.input_file, self.ecutwfc_increment)
                else:
                    qe_io.sum_ecut(self.input_file, self.ecutwfc_increment)

                previous_energy = total_energy
                self.energy += self.ecutwfc_increment
                iteration += 1

            except Exception as e:
                self.progress_signal.emit(f'Error during iteration {iteration+1}: {str(e)}')
                traceback.print_exc()
                break

        if not converged:
            self.finished_signal.emit(f'Convergence was not reached after {self.max_iterations} iterations. Last total energy={total_energy}')
        else:
            self.finished_signal.emit(f'Convergence achieved in iteration {iteration+1} with ecutwfc={self.energy} and total energy={total_energy}')

class RunQESingleThread(QThread):
    progress_signal = Signal(str)
    finished_signal = Signal(str)

    def __init__(self, input_file, parent=None):
        super().__init__(parent)
        self.input_file = input_file

    def run(self):
        try:
            run = RunQuantumEspresso()
            output_file = self.input_file.replace('.in', '.out')
            run.run_qe_process(self.input_file, output_file)
            self.finished_signal.emit(f"Quantum ESPRESSO process executed successfully. Output file: {output_file}")
        except Exception as e:
            self.progress_signal.emit(f"Error: {str(e)}")
            traceback.print_exc()

class ArchiveWindow(QMainWindow):
    message_signal = Signal(str)

    def __init__(self):
        super().__init__()

        # Configurar la ventana principal
        self.setGeometry(100, 100, 1200, 800)
        self.setWindowTitle('Iteration Files')

        # Crear una instancia de LoadingDialog
        self.loading_dialog = LoadingDialog(self)

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
        self.k_points_label = QLabel("K Points:")
        controls_layout.addWidget(self.k_points_label)

        self.grid_layout = QGridLayout()
        self.k_points_editors = [QLineEdit() for _ in range(6)]  # Fila de 6 QLineEdit
        for i in range(6):
            self.k_points_editors[i].setStyleSheet("background-color: white")  # Establecer fondo blanco
            self.grid_layout.addWidget(self.k_points_editors[i], 0, i)
        controls_layout.addLayout(self.grid_layout)

        # Crear el campo para la energía de corte
        self.energy_label = QLabel("Kinetic energy cutoff (Ry):")
        controls_layout.addWidget(self.energy_label)
        self.energy_editor = QLineEdit()
        self.energy_editor.setStyleSheet("background-color: white")  # Establecer fondo blanco
        controls_layout.addWidget(self.energy_editor)

        # Crear el campo para el número de iteraciones
        self.iterations_label = QLabel("Number of Iterations:")
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
        message = self.validate_inputs()
        if message != '':
            self.message(message)
            return

        # Mostrar el diálogo de carga
        self.loading_dialog.show()

        qe_io.modify_ecut(input_file, float(self.energy_editor.text()))
        qe_io.modify_k_points(input_file, self.get_k_points())

        self.thread = RunQESingleThread(input_file)
        self.thread.progress_signal.connect(self.message)
        self.thread.finished_signal.connect(self.handle_run_single_finished)
        self.thread.start()

    def handle_run_single_finished(self, result):
        # Ocultar el diálogo de carga
        self.loading_dialog.hide()
        self.message(result)

    def iterate_qe(self, input_file):
        message = self.validate_inputs()
        if message != '':
            self.message(message)
            return

        # Mostrar el diálogo de carga
        self.loading_dialog.show()

        k_points = self.get_k_points()
        energy = float(self.energy_editor.text())
        max_iterations = int(self.iterations_editor.text())

        qe_io.modify_ecut(input_file, energy)
        qe_io.modify_k_points(input_file, k_points)

        self.thread = QEIterateThread(self, input_file, k_points, energy, max_iterations)
        self.thread.progress_signal.connect(self.message)
        self.thread.finished_signal.connect(self.handle_iterate_finished)
        self.thread.start()

    def handle_iterate_finished(self, result):
        # Ocultar el diálogo de carga
        self.loading_dialog.hide()
        self.message(result)
        self.plot_energies(self.thread.ecutwfc_list, self.thread.energy_list)

    def validate_inputs(self):
        message = ''
        if not (all(editor.text().strip() for editor in self.k_points_editors) and self.energy_editor.text().strip() and self.iterations_editor.text().strip()):
            message = 'Please fill out all fields for K points, cutting energy and number of iterations.'
            return message
        try:
            k_points = self.get_k_points()
            for k_point in k_points:
                k_point[0] = int(k_point[0])
                k_point[1] = int(k_point[1])
                k_point[2] = int(k_point[2])

            energy_cutoff = float(self.energy_editor.text())
            max_iterations = int(self.iterations_editor.text())
        except:
            message = 'K Points need to be integer number\nKinetic energy cutoff need to be float number\nIterations need to be integer number'
            return message
        
        if any(k_point >= 10 for k_point in k_points[0]) or any(k_point < 0 for k_point in k_points[0]):
            message += 'nk1, nk2, nk3 in K_POINTS they need to be more than 0 and less than 10\n'
        
        if any(k_point > 1 for k_point in k_points[1]) or any(k_point < 0 for k_point in k_points[1]):
            message += 'sk1, sk2, sk3 in K_POINTS they need to be 0 or 1\n'
        
        if energy_cutoff < 0 or energy_cutoff >= 200:
            message += 'Energy cutoff need to be more than 0 and less than 200\n'
        
        if max_iterations < 0 or max_iterations >= 30:
            message += 'Iterations need to be more than 0 and less than 30\n'

        return message

    def get_k_points(self):
        k_points = []
        for i in range(0, 6, 3):
            k_points.append([int(self.k_points_editors[i].text()), int(self.k_points_editors[i+1].text()), int(self.k_points_editors[i+2].text())])
        return k_points

    def run_qe(self, input_file, output_file):
        run = RunQuantumEspresso()
        run.run_qe_process(input_file, output_file)
        self.message(f"Quantum ESPRESSO process executed. Output file: {output_file}")

    def plot_energies(self, ecutwfc_list, energy_list):
        plt.figure(figsize=(10, 6))
        plt.plot(ecutwfc_list, energy_list, marker='o')
        plt.xlabel('Kinetic energy cutoff (Ry)')
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
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    archive_window = ArchiveWindow()
    archive_window.show()
    return app

if __name__ == "__main__":
    app = run_archive()
    sys.exit(app.exec_())

