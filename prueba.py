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

from PySide2.QtWidgets import (QApplication, QMainWindow, QPushButton, QPlainTextEdit,
                                QVBoxLayout, QWidget, QLabel, QLineEdit, QGridLayout,
                                QFileDialog, QDialog, QTextEdit, QSizePolicy, QMessageBox)
from PySide2.QtGui import QPixmap, QColor, QImageReader, QImage, QFont
from PySide2.QtCore import Qt, Signal, QSize , QObject, QThread
import sys
import subprocess
import quantum_espresso_io as qe_io

class QeWorker(QObject):
    message_signal = Signal(str)

    def run_qe_process(self):
        for i in range(5):
            self.message_signal.emit(f"Execute process {i}")
            print(f"Execute process {i}")
            subprocess.run('pw.x < ./C.in > ./C.out', shell=True)
            self.message_signal.emit("Execute Finish")

            energy = qe_io.find_total_energy("C.out")
            self.message_signal.emit(f"La energía total en la iteración {i} es: {energy}")
            print(f"La energía total en la iteración {i} es: {energy}")

            try:
                qe_io.sum_ecut("C.in", 5)
            except Exception as e:
                self.message_signal.emit(f"Something is wrong: {e}")


class MainWindow(QMainWindow):
    message_signal = Signal(str)

    def __init__(self):
        super().__init__()

        # Configurar la ventana principal
        self.setGeometry(100, 100, 1200, 800)  
        self.setWindowTitle('Mi Aplicación PySide2')

        # Crear un layout vertical
        layout = QVBoxLayout()

        # Agregar una etiqueta con una imagen en la esquina
        imagen_path = 'GUI.jpeg'  
        imagen_label = QLabel(self)

        # Verificar que la imagen se cargue correctamente
        pixmap = QPixmap(imagen_path)
        if not pixmap.isNull():
            # Ajustar el tamaño de la imagen
            pixmap = pixmap.scaled(300, 300, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            imagen_label.setPixmap(pixmap)
        else:
            print("Error: No se pudo cargar la imagen.")

        imagen_label.setAlignment(Qt.AlignTop | Qt.AlignLeft)  
        layout.addWidget(imagen_label)

        # Crear un widget para los controles
        controls_widget = QWidget()
        controls_layout = QVBoxLayout(controls_widget)

        # Agregar un botón y un cuadro de texto para ejecutar QE
        self.btn_qe = QPushButton("Ejecutar QE")
        self.btn_qe.pressed.connect(self.start_process)
        self.btn_qe.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.btn_qe.setStyleSheet("background-color: #a8a8a8; color: black; font-size: 18px;")
        controls_layout.addWidget(self.btn_qe)

        # Agregar un botón y un cuadro de texto para iterar QE
        self.btn_iterate_qe = QPushButton("Iterar QE")
        self.btn_iterate_qe.pressed.connect(self.iterate_process)
        self.btn_iterate_qe.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.btn_iterate_qe.setStyleSheet("background-color: #a8a8a8; color: black; font-size: 18px;")
        controls_layout.addWidget(self.btn_iterate_qe)

        # Crear la matriz para los puntos K
        self.k_points_label = QLabel("Puntos K:")
        controls_layout.addWidget(self.k_points_label)

        self.grid_layout = QGridLayout()
        self.k_points_editors = [[QLineEdit() for _ in range(3)] for _ in range(2)]  # Matriz 2x3 de QLineEdit
        for i in range(2):
            for j in range(3):
                self.k_points_editors[i][j].setStyleSheet("background-color: white")  # Establecer fondo blanco
                self.grid_layout.addWidget(self.k_points_editors[i][j], i, j)
        controls_layout.addLayout(self.grid_layout)

        # Crear el campo para la energía de corte
        self.energy_label = QLabel("Energía de Corte:")
        controls_layout.addWidget(self.energy_label)
        self.energy_editor = QLineEdit()
        self.energy_editor.setStyleSheet("background-color: white")  # Establecer fondo blanco
        controls_layout.addWidget(self.energy_editor)

        # Botón para guardar y modificar el archivo
        self.save_btn = QPushButton("Guardar y Modificar Archivo")
        self.save_btn.pressed.connect(self.save_and_modify_file)
        self.save_btn.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.save_btn.setStyleSheet("background-color: #a8a8a8; color: black; font-size: 18px;")
        controls_layout.addWidget(self.save_btn)

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

    def start_process(self):
        self.message("Executing process.")
        subprocess.run('pw.x < ./C.in > ./C.out', shell=True)
        self.message("Execute Finish")
    
    # Iteration process with thread
    def iterate_process(self):
        self.qe_worker = QeWorker()
        self.qe_worker.message_signal.connect(self.message_signal.emit)

        self.qe_thread = QThread()
        self.qe_worker.moveToThread(self.qe_thread)

        self.qe_thread.started.connect(self.qe_worker.run_qe_process)
        self.qe_thread.finished.connect(self.qe_thread.quit)

        self.qe_thread.start()

    def save_and_modify_file(self):
        # Obtener el archivo .in
        file_path, _ = QFileDialog.getOpenFileName(self, "Seleccionar archivo .in", "", "Archivo .in (*.in)")
        if not file_path:
            self.message("No se seleccionó ningún archivo .in.")
            return

        # Leer el contenido del archivo .in y mostrarlo en una ventana de texto
        with open(file_path, 'r') as file:
            content = file.read()
            dialog = FileContentDialog(self, content)
            dialog.exec_()


class FileContentDialog(QDialog):
    def __init__(self, parent, content):
        super(FileContentDialog, self).__init__(parent)
        self.setWindowTitle("Contenido del Archivo .in")
        self.setGeometry(100, 100, 800, 600)

        layout = QVBoxLayout()

        self.text_edit = QTextEdit()
        self.text_edit.setPlainText(content)
        layout.addWidget(self.text_edit)

        self.setLayout(layout)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
