###################################################
##                                               ##
##    main_screen.py                             ##
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


import os
from PySide2.QtWidgets import QApplication, QMainWindow, QHBoxLayout, QWidget, QLabel, QVBoxLayout, QListWidget, QListWidgetItem, QFileDialog, QLineEdit, QPushButton, QDialog
from PySide2.QtGui import QPixmap, QColor, QIcon
from PySide2.QtCore import Qt
from screens.parameters_screen import run_parameters
from screens.archive import run_archive

RECENT_FILES_PATH = './file_operations/recent_files.txt'

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
        dialog = PopupWindow(self)
        if dialog.exec_():
            entered_text = dialog.get_text()
            if entered_text:
                directory = self.select_folder()
                if directory:
                    run_parameters(entered_text, directory, False)
                    self.add_recent_file(os.path.join(directory, entered_text + ".qg"))

    def open_clicked(self, event=None):
        file_path = self.select_file()
        if file_path:
            file_name = os.path.splitext(os.path.basename(file_path))[0]
            directory = os.path.dirname(file_path)
            run_parameters(file_name, directory, True)
            self.add_recent_file(file_path)

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

    def add_recent_file(self, file_path):
        with open(RECENT_FILES_PATH, 'a') as f:
            f.write(file_path + '\n')


class PopupWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Project Name")
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()

        self.text_input = QLineEdit()
        layout.addWidget(self.text_input)

        create_button = QPushButton("Create")
        create_button.clicked.connect(self.accept)
        layout.addWidget(create_button, alignment=Qt.AlignCenter)

        self.setLayout(layout)

    def get_text(self):
        return self.text_input.text()


class WelcomeWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        AppStyle.apply(self)
        layout = QVBoxLayout()

        header_layout = QHBoxLayout()

        top_left_image_label = QLabel(self)
        pixmap_top_left = QPixmap('./screens/images/ESCOM.png')
        if not pixmap_top_left.isNull():
            pixmap_top_left = pixmap_top_left.scaled(250, 250, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            top_left_image_label.setPixmap(pixmap_top_left)
        header_layout.addWidget(top_left_image_label, alignment=Qt.AlignLeft)

        top_right_image_label = QLabel(self)
        pixmap_top_right = QPixmap('./screens/images/GUI.jpeg')
        if not pixmap_top_right.isNull():
            pixmap_top_right = pixmap_top_right.scaled(100, 100, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            top_right_image_label.setPixmap(pixmap_top_right)
        header_layout.addWidget(top_right_image_label, alignment=Qt.AlignRight)

        layout.addLayout(header_layout)

        welcome_label = QLabel("Welcome", self)
        welcome_label.setAlignment(Qt.AlignCenter)
        welcome_label.setStyleSheet("font-size: 24px; font-weight: bold;")
        layout.addWidget(welcome_label)

        self.recent_documents_list = QListWidget(self)
        self.recent_documents_list.itemClicked.connect(self.handle_file_click)
        layout.addWidget(self.recent_documents_list)

        self.show_in_files()

        buttons_layout = QHBoxLayout()

        tutorial_button = QLabel(self)
        tutorial_button.setPixmap(QPixmap('./screens/images/tutorial.png').scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        tutorial_button.mousePressEvent = self.tutorial_clicked
        buttons_layout.addWidget(tutorial_button)

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
        if os.path.exists(RECENT_FILES_PATH):
            with open(RECENT_FILES_PATH, 'r') as f:
                recent_files = f.readlines()

            for file_path in recent_files:
                file_path = file_path.strip()
                if os.path.exists(file_path):
                    filename = os.path.basename(file_path)
                    item = QListWidgetItem(QIcon('./screens/images/file_icon.png'), filename)
                    item.setToolTip(file_path)
                    self.recent_documents_list.addItem(item)

    def handle_file_click(self, item):
        file_path = item.toolTip()
        if file_path.endswith('.qg'):
            file_name = os.path.splitext(os.path.basename(file_path))[0]
            directory = os.path.dirname(file_path)
            run_parameters(file_name, directory, True)
        else:
            with open(file_path, 'r') as file:
                content = file.read()
                print(f"Contents of {item.text()}:\n{content}")

    def tutorial_clicked(self, event):
        print("Tutorial button clicked")

    def nuevo_calculo(self, event):
        self.archive = run_archive()

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



