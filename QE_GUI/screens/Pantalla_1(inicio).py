import os
from PySide2.QtWidgets import QApplication, QMainWindow, QHBoxLayout, QWidget, QLabel, QVBoxLayout, QListWidget, QListWidgetItem, QFileDialog
from PySide2.QtGui import QPixmap, QColor, QIcon
from PySide2.QtCore import Qt

class AppStyle:
    @staticmethod
    def apply(window):
        window.setGeometry(100, 100, 1200, 800)  
        window.setWindowTitle('')
        color_fondo = QColor(135, 158, 197)
        window.setStyleSheet(f'background-color: {color_fondo.name()};')

class CreateOpenWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setGeometry(400, 400, 400, 200)
        self.setWindowTitle('Create/Open Window')
        color_fondo = QColor(135, 158, 197)
        self.setStyleSheet(f'background-color: {color_fondo.name()};')

        layout = QVBoxLayout()

        button_layout = QHBoxLayout()

        create_button = QLabel(self)
        create_button.setPixmap(QPixmap('./images/create.png').scaled(150, 150, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        create_button.mousePressEvent = self.create_clicked  
        button_layout.addWidget(create_button)

        button_layout.addStretch(1)

        open_button = QLabel(self)
        open_button.setPixmap(QPixmap('./images/open.png').scaled(150, 150, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        open_button.mousePressEvent = self.open_clicked
        button_layout.addWidget(open_button)

        layout.addLayout(button_layout)

        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

    def create_clicked(self, event):  
        print("Create button clicked")

    def open_clicked(self, event):
        print("Open button clicked")



class WelcomeWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        AppStyle.apply(self)
        layout = QVBoxLayout()

        imagen_label = QLabel(self)
        pixmap = QPixmap('./images/GUI.jpeg')
        if not pixmap.isNull():
            pixmap = pixmap.scaled(200, 150, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            imagen_label.setPixmap(pixmap)
        layout.addWidget(imagen_label, alignment=Qt.AlignTop | Qt.AlignRight)

        welcome_label = QLabel("Bienvenido", self)
        welcome_label.setAlignment(Qt.AlignCenter)
        welcome_label.setStyleSheet("font-size: 24px; font-weight: bold;")
        layout.addWidget(welcome_label)

        self.recent_documents_list = QListWidget(self)
        layout.addWidget(self.recent_documents_list)

        self.show_in_files()

        buttons_layout = QHBoxLayout()

        tutorial_button = QLabel(self)
        tutorial_button.setPixmap(QPixmap('./images/tutorial.png').scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        tutorial_button.mousePressEvent = self.tutorial_clicked
        buttons_layout.addWidget(tutorial_button)

        subir_button = QLabel(self)
        subir_button.setPixmap(QPixmap('./images/upload.png').scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        subir_button.mousePressEvent = self.subir_archivo
        buttons_layout.addWidget(subir_button)

        nuevo_calculo_button = QLabel(self)
        nuevo_calculo_button.setPixmap(QPixmap('./images/new_calculation.png').scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        nuevo_calculo_button.mousePressEvent = self.nuevo_calculo
        buttons_layout.addWidget(nuevo_calculo_button)

        proyect_button = QLabel(self)
        proyect_button.setPixmap(QPixmap('./images/proyect.png').scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        proyect_button.mousePressEvent = self.show_create_open_window
        buttons_layout.addWidget(proyect_button)

        layout.addLayout(buttons_layout)

        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

    def show_in_files(self):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        for filename in os.listdir(current_dir):
            if filename.endswith(".in"):
                item = QListWidgetItem(QIcon('./images/file_icon.png'), filename)
                self.recent_documents_list.addItem(item)

    def tutorial_clicked(self, event):
        print("Tutorial button clicked")

    def subir_archivo(self, event):
        file_dialog = QFileDialog(self)
        file_dialog.setNameFilter("Archivos (*.in)")
        if file_dialog.exec_():
            file_path = file_dialog.selectedFiles()[0]
            print("Archivo seleccionado:", file_path)

    def nuevo_calculo(self, event):
        print("Nuevo cálculo button clicked")

    def show_create_open_window(self, event):
        self.create_open_window = CreateOpenWindow()
        self.create_open_window.show()

if __name__ == '__main__':
    app = QApplication([])
    window = WelcomeWindow()
    window.show()
    app.exec_()

