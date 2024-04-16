import sys

from PySide2.QtWidgets import QApplication, QMainWindow, QLabel, QVBoxLayout, QAction, QPushButton, QWidget, QMenuBar, QDialog, QLineEdit
from PySide2.QtGui import QPixmap,QColor
from PySide2.QtCore import Qt
from screens.parameters_screen import run_parameters

class AppStyle:
    @staticmethod
    def apply(window):
        # Configurar el tamaño de la ventana
        window.setGeometry(100, 100, 1200, 800)  
        color_fondo = QColor(135, 158, 197)
        window.setStyleSheet(f'background-color: {color_fondo.name()};')

class ProjectWIndow(QWidget):
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
        title_label = QLabel("Project")
        title_label.setStyleSheet("font-size: 24px; font-weight: bold;")
        layout.addWidget(title_label, alignment=Qt.AlignLeft)

        # Menú de opciones
        menu_bar = self.create_menu_bar()
        layout.addWidget(menu_bar)

        # Botón "Seleccionar Proyecto"
        select_button = QPushButton("Seleccionar Proyecto")
        select_button.setStyleSheet("font-size: 18px; padding: 10px 20px;")
        layout.addWidget(select_button, alignment=Qt.AlignCenter)

        # Botón "CrearProyecto"
        create_button = QPushButton("Create Proyect")
        create_button.setStyleSheet("font-size: 18px; padding: 10px 20px;")
        create_button.clicked.connect(self.show_popup)
        layout.addWidget(create_button, alignment=Qt.AlignCenter)

        # Establecer el layout principal
        self.setLayout(layout)

    def create_menu_bar(self):
        menu_bar = QMenuBar()

        # Opciones del menú
        crear_action = QAction("Create", self)
        abrir_action = QAction("Open", self)
        borrar_action = QAction("Delete", self)

        # Menú "Proyecto" con las opciones
        proyecto_menu = menu_bar.addMenu("Project")
        proyecto_menu.addAction(crear_action)
        proyecto_menu.addAction(abrir_action)
        proyecto_menu.addAction(borrar_action)

        return menu_bar
    
    def show_popup(self):
        # Create and show the popup window
        dialog = PopupWindow(self)
        if dialog.exec_():
            entered_text = dialog.get_text()
            # If the user enters the name of the project, the parameters screen is executed
            if entered_text:
                run_parameters(entered_text)
                


class PopupWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Project Name")
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()

        # Text field to enter text
        self.text_input = QLineEdit()
        layout.addWidget(self.text_input)

        # "Create" Button
        create_button = QPushButton("Create")
        create_button.clicked.connect(self.accept)
        layout.addWidget(create_button, alignment=Qt.AlignCenter)

        self.setLayout(layout)

    def get_text(self):
        # Return the text entered by the user in the QLineEdit
        return self.text_input.text()

def run_project():
    app = QApplication(sys.argv)
    window = ProjectWIndow()
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    run_project()
