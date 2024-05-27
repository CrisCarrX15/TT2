import psutil
import subprocess
from PySide2.QtCore import QThread, Signal

# Run Quantum ESPRESSO with another thread
class RunQuantumEspressoThread(QThread):
    finished = Signal()
    def __init__(self, file_path_in, file_path_out):
        super().__init__()
        self.file_path_in = file_path_in
        self.file_path_out = file_path_out

    def run(self):
        subprocess.run(f'pw.x < {self.file_path_in} > {self.file_path_out}', shell=True)
        self.finished.emit()


class RunQuantumEspresso():
    def run_qe_process(self, file_path_in, file_path_out):
            #print(f"Execute process")
            subprocess.run(f'pw.x < {file_path_in} > {file_path_out}', shell=True)
    
    def run_qe_process_thread(self, file_path_in, file_path_out, callback):
        self.thread = RunQuantumEspressoThread(file_path_in, file_path_out)
        self.thread.finished.connect(callback)
        self.thread.start()
    

    def run_qe_cores(self, file_path_in, file_path_out):
        # Obtener información del sistema
        num_cores, total_memory = self.get_system_info()
        
        # Ajustar los parámetros según las características del sistema
        if total_memory >= 16:
            mpi_processes = min(num_cores, 4)  # Utilizar hasta 4 procesos MPI si hay más de 16 GB de RAM
        else:
            mpi_processes = min(num_cores, 2)  # Utilizar hasta 2 procesos MPI si hay menos de 16 GB de RAM
        
        #print(f'mpi: {mpi_processes}')

        subprocess.run(f'mpirun -np {mpi_processes} pw.x < {file_path_in} > {file_path_out}', shell=True)
        
        #print("Cores identificados correctamente.")
    
    def get_system_info(self):
        # Obtener el número de núcleos de la CPU, logical=False -> núcleos físicos
        num_cores = psutil.cpu_count(logical=False)
        
        # Obtener la cantidad de memoria RAM disponible en GB
        total_memory = psutil.virtual_memory().total / (1024**3)
        
        return num_cores, total_memory