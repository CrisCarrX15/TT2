from PySide2.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QLabel, QPushButton, QGroupBox, QVBoxLayout, QAction, QMenu, QMenuBar, QGridLayout
from PySide2.QtGui import QPixmap, QColor
from PySide2.QtCore import Qt

class WelcomeWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # Configurar la ventana de bienvenida
        self.setGeometry(200, 200, 600, 400)
        self.setWindowTitle('')

        # Crear un layout vertical
        layout = QVBoxLayout()

        # Agregar una etiqueta con el texto "Bienvenido"
        welcome_label = QLabel("Bienvenido", self)
        welcome_label.setAlignment(Qt.AlignCenter)
        welcome_label.setStyleSheet("font-size: 24px; font-weight: bold;")
        layout.addWidget(welcome_label)

        # Agregar una etiqueta con la imagen GUI
        imagen_path = 'GUI.jpeg'
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

        # Establecer el color de fondo
        color_fondo = QColor(135, 158, 197)
        self.setStyleSheet(f'background-color: {color_fondo.name()};')

        # Crear un widget central y establecer el layout
        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

if __name__ == '__main__':
    app = QApplication([])
    window = WelcomeWindow()
    window.show()
    app.exec_()
