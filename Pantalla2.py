import sys
from PySide2.QtWidgets import QApplication, QMainWindow, QLabel, QVBoxLayout, QAction, QPushButton, QWidget, QMenuBar
from PySide2.QtGui import QPixmap
from PySide2.QtCore import Qt

class ProyectoWindow(QWidget):
    def __init__(self):
        super().__init__()

        # Configurar la ventana principal
        self.setWindowTitle("Proyecto")
        self.setGeometry(100, 100, 800, 600)
        self.setStyleSheet("background-color: #879ec5;")

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

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ProyectoWindow()
    window.show()
    sys.exit(app.exec_())
