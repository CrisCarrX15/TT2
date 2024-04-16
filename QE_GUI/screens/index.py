import sys
from PySide2.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QWidget, QLabel, QPushButton, 
                                QGroupBox, QAction, QMenu, QMenuBar, QGridLayout, QPlainTextEdit, QLineEdit,
                                QTextEdit, QFileDialog, QDialog, QHBoxLayout, QTabWidget)
from PySide2.QtGui import QPixmap, QColor
from PySide2.QtCore import Qt, Signal, QObject, QThread

import subprocess
import file_operations.quantum_espresso_io as qe_io

class QeWorker(QObject):
    message_signal = Signal(str)

    def run_qe_process(self):
        for i in range(5):
            self.message_signal.emit(f"Execute process {i}")
            print(f"Execute process {i}")
            subprocess.run('pw.x < ./test_files/C.in > ./test_files/C.out', shell=True)
            self.message_signal.emit("Execute Finish")

            energy = qe_io.find_total_energy("./test_files/C.out")
            self.message_signal.emit(f"La energía total en la iteración {i} es: {energy}")
            print(f"La energía total en la iteración {i} es: {energy}")

            try:
                qe_io.sum_ecut("./test_files/C.in", 5)
            except Exception as e:
                self.message_signal.emit(f"Something is wrong: {e}")

class AppStyle:
    @staticmethod
    def apply(window):
        # Establecer el tamaño de la ventana
        window.setGeometry(100, 100, 1200, 800)  
        window.setWindowTitle('Mi Aplicación PySide2')

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

class WelcomeWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # Aplicar el estilo
        AppStyle.apply(self)

        # Crear un layout vertical
        layout = QVBoxLayout()

        # Agregar una etiqueta con el texto "Bienvenido"
        welcome_label = QLabel("Bienvenido", self)
        welcome_label.setAlignment(Qt.AlignCenter)
        welcome_label.setStyleSheet("font-size: 24px; font-weight: bold;")
        layout.addWidget(welcome_label)

        # Agregar una etiqueta con la imagen GUI
        imagen_path = './images/GUI.jpeg'  
        imagen_label = QLabel(self)
        pixmap = QPixmap(imagen_path)
        if not pixmap.isNull():
            pixmap = pixmap.scaled(300, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            imagen_label.setPixmap(pixmap)
        layout.addWidget(imagen_label, alignment=Qt.AlignCenter)

        # Crear un grupo para contener el menú y el botón
        group_box = QGroupBox()
        group_layout = QVBoxLayout(group_box)

        # Agregar un menú en el grupo
        menu_bar = QMenuBar(self)
        proyecto_menu = menu_bar.addMenu("Proyecto")

        tutorial_action = QAction("Tutorial", self)
        subir_action = QAction("Subir Archivo", self)
        crear_action = QAction("Crear", self)
        abrir_action = QAction("Abrir", self)

        proyecto_menu.addAction(crear_action)
        proyecto_menu.addAction(abrir_action)

        menu_bar.addAction(subir_action)
        menu_bar.addAction(tutorial_action)

        group_layout.addWidget(menu_bar)

        # Agregar un botón de "Aceptar" en el grupo
        accept_button = QPushButton("Aceptar", self)
        group_layout.addWidget(accept_button, alignment=Qt.AlignCenter)

        layout.addWidget(group_box)

        # Crear un widget central y establecer el layout
        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

class ProyectoWindow(QWidget):
    def __init__(self):
        super().__init__()

        # Aplicar el estilo
        AppStyle.apply(self)

        # Layout principal
        layout = QVBoxLayout()

        # Imagen de GUI en la esquina superior derecha
        image_label = QLabel()
        pixmap = QPixmap("./images/GUI.jpeg")
        pixmap = pixmap.scaledToHeight(100)
        image_label.setPixmap(pixmap)
        layout.addWidget(image_label, alignment=Qt.AlignRight)

        # Etiqueta del título en la esquina superior izquierda
        title_label = QLabel("Proyecto")
        title_label.setStyleSheet("font-size: 24px; font-weight: bold;")
        layout.addWidget(title_label, alignment=Qt.AlignLeft)

        # Menú de opciones
        menu_bar = self.create_menu_bar()
        layout.addWidget(menu_bar)

        # Botón "Seleccionar Proyecto"
        select_button = QPushButton("Seleccionar Proyecto")
        select_button.setStyleSheet("font-size: 18px; padding: 10px 20px;")
        layout.addWidget(select_button, alignment=Qt.AlignCenter)

        # Establecer el layout principal
        self.setLayout(layout)

    def create_menu_bar(self):
        menu_bar = QMenuBar()

        # Opciones del menú
        crear_action = QAction("Crear", self)
        abrir_action = QAction("Abrir", self)
        borrar_action = QAction("Borrar", self)

        # Menú "Proyecto" con las opciones
        proyecto_menu = menu_bar.addMenu("Proyecto")
        proyecto_menu.addAction(crear_action)
        proyecto_menu.addAction(abrir_action)
        proyecto_menu.addAction(borrar_action)

        return menu_bar

class MainWindow(QMainWindow):
    message_signal = Signal(str)

    def __init__(self):
        super().__init__()

        # Aplicar el estilo
        AppStyle.apply(self)

        # Crear un layout vertical
        layout = QVBoxLayout()

        # Agregar una etiqueta con una imagen en la esquina
        imagen_path = './images/GUI.jpeg'  
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
        controls_layout.addWidget(self.btn_qe)

        # Agregar un botón y un cuadro de texto para iterar QE
        self.btn_iterate_qe = QPushButton("Iterar QE")
        self.btn_iterate_qe.pressed.connect(self.iterate_process)
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
        controls_layout.addWidget(self.save_btn)

        layout.addWidget(controls_widget)

        self.text = QPlainTextEdit()
        self.text.setReadOnly(True)
        layout.addWidget(self.text)

        # Crear un widget central y establecer el layout
        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

    def message(self, s):
        self.text.appendPlainText(s)

    def start_process(self):
        self.message("Executing process.")
        subprocess.run('pw.x < ./test_files/C.in > ./test_files/C.out', shell=True)
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

        # Recuperar los valores de los puntos K
        k_points = []
        for row in self.k_points_editors:
            k_point = []
            for editor in row:
                text = editor.text()
                if text.strip():
                    try:
                        value = int(text)
                        k_point.append(value)
                    except ValueError:
                        pass  # Ignorar valores no enteros
            if len(k_point) == 3:
                k_points.append(k_point)

        # Recuperar el valor de la energía de corte
        energy_text = self.energy_editor.text()
        try:
            energy = float(energy_text)

            # Modificar el archivo con los nuevos valores
            qe_io.modify_k_points(file_path, k_points)
            qe_io.modify_ecut(file_path, energy)

            self.message("Archivo modificado y guardado.")
        except Exception as e:
            print("Error al modificar el archivo:", e)
        
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

    tab_widget = QTabWidget()
    tab_widget.addTab(WelcomeWindow(), "Bienvenida")
    tab_widget.addTab(ProyectoWindow(), "Proyecto")
    tab_widget.addTab(MainWindow(), "Principal")

    tab_widget.show()

    sys.exit(app.exec_())
