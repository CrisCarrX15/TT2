import os
from PySide2.QtWidgets import QApplication, QMainWindow, QHBoxLayout, QWidget, QLabel, QVBoxLayout, QListWidget, QListWidgetItem, QFileDialog, QLineEdit, QPushButton, QDialog
from PySide2.QtGui import QPixmap, QColor, QIcon
from PySide2.QtCore import Qt
from screens.parameters_screen import run_parameters

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
        self.setWindowTitle('Create/Open')
        color_fondo = QColor(135, 158, 197)
        self.setStyleSheet(f'background-color: {color_fondo.name()};')

        layout = QVBoxLayout()

        button_layout = QHBoxLayout()

        create_button = QLabel(self)
        create_button.setPixmap(QPixmap('./screens/images/create.png').scaled(150, 150, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        create_button.mousePressEvent = self.create_clicked
        #create_button.clicked.connect(self.create_clicked) 
        button_layout.addWidget(create_button)

        button_layout.addStretch(1)

        open_button = QLabel(self)
        open_button.setPixmap(QPixmap('./screens/images/open.png').scaled(150, 150, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        open_button.mousePressEvent = self.open_clicked
        button_layout.addWidget(open_button)

        layout.addLayout(button_layout)

        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

    def create_clicked(self, event=None):  
        # Create and show the popup window
        dialog = PopupWindow(self)
        if dialog.exec_():
            entered_text = dialog.get_text()
            # If the user enters the name of the project, the parameters screen is executed
            if entered_text:
                directory = self.select_folder()
                if directory:
                    run_parameters(entered_text, directory, False)
                
    def open_clicked(self, event=None):
        file_path = self.select_file()
        if file_path:
            file_name = os.path.splitext(os.path.basename(file_path))[0]
            directory = os.path.dirname(file_path)
            run_parameters(file_name, directory, True)
    
    def select_file(self):
        file_dialog = QFileDialog(self)
        file_dialog.setNameFilter("Files (*.qg)")
        filename = None
        if file_dialog.exec_():
            file_path = file_dialog.selectedFiles()[0]
            print("Archivo seleccionado:", file_path)
            return file_path
        return filename
    
    def select_folder(self):
        folder_dialog = QFileDialog(self)
        folder_dialog.setFileMode(QFileDialog.Directory)
        folder_dialog.setOption(QFileDialog.ShowDirsOnly, True)
        folder_dialog.setOption(QFileDialog.DontResolveSymlinks, True)
        directory = None
        if folder_dialog.exec_():
            directory = folder_dialog.selectedFiles()[0]
            print("Directorio seleccionado:", directory)
            return directory
        return directory


# ========== POPUP TO CREATE NEW PROJECT ==========
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


# ========== MAIN WINDOW ========== 
class WelcomeWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        AppStyle.apply(self)
        layout = QVBoxLayout()

        imagen_label = QLabel(self)
        pixmap = QPixmap('./screens/images/GUI.jpeg')
        if not pixmap.isNull():
            pixmap = pixmap.scaled(200, 150, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            imagen_label.setPixmap(pixmap)
        layout.addWidget(imagen_label, alignment=Qt.AlignTop | Qt.AlignRight)

        welcome_label = QLabel("Welcome", self)
        welcome_label.setAlignment(Qt.AlignCenter)
        welcome_label.setStyleSheet("font-size: 24px; font-weight: bold;")
        layout.addWidget(welcome_label)

        self.recent_documents_list = QListWidget(self)
        layout.addWidget(self.recent_documents_list)

        self.show_in_files()

        buttons_layout = QHBoxLayout()

        tutorial_button = QLabel(self)
        tutorial_button.setPixmap(QPixmap('./screens/images/tutorial.png').scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        tutorial_button.mousePressEvent = self.tutorial_clicked
        buttons_layout.addWidget(tutorial_button)

        subir_button = QLabel(self)
        subir_button.setPixmap(QPixmap('./screens/images/upload.png').scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        subir_button.mousePressEvent = self.subir_archivo
        buttons_layout.addWidget(subir_button)

        nuevo_calculo_button = QLabel(self)
        nuevo_calculo_button.setPixmap(QPixmap('./screens/images/new_calculation.png').scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        nuevo_calculo_button.mousePressEvent = self.nuevo_calculo
        buttons_layout.addWidget(nuevo_calculo_button)

        proyect_button = QLabel(self)
        proyect_button.setPixmap(QPixmap('./screens/images/project.png').scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation))
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
                item = QListWidgetItem(QIcon('./screens/images/file_icon.png'), filename)
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

def run_main():
    app = QApplication([])
    window = WelcomeWindow()
    window.show()
    app.exec_()

if __name__ == '__main__':
    run_main()

