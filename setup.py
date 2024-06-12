from setuptools import setup, find_packages

def main():
    import QE_GUI.main
    QE_GUI.main.main()

setup(
    name='QE_GUI',
    version='1.0.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'importlib-metadata==4.6.4',
        'numpy==1.26.3',
        'packaging==23.2',
        'paramiko==2.9.3',
        'pexpect==4.8.0',
        'Pillow==9.0.1',
        'psutil==5.9.8',
        'ptyprocess==0.7.0',
        'PySide2==5.15.2.1',
        'python-dateutil==2.8.1',
        'python-time==0.3.0'
    ],
    entry_points={
        'console_scripts': [
            'qg_gui = QE_GUI.main:main',  # Ajusta según tu punto de entrada
        ],
    },
)
