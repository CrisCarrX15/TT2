import psutil
import subprocess

class RunQuantumEspresso():
    def run_qe_process(self, file_path_in, file_path_out):
            subprocess.run(f'pw.x < {file_path_in} > {file_path_out}', shell=True)

    def run_qe_cores(self, file_path_in, file_path_out):
        # Get system information
        num_cores, total_memory = self.get_system_info()
        
        if total_memory >= 16:
            mpi_processes = min(num_cores, 4)  # Use up to 4 MPI processes if more than 16 GB of RAM
        else:
            mpi_processes = min(num_cores, 2)  # Use up to 2 MPI processes if there is less than 16 GB of RAM

        subprocess.run(f'mpirun -np {mpi_processes} pw.x < {file_path_in} > {file_path_out}', shell=True)
    
    def get_system_info(self):
        # Get the number of CPU cores, logical=False -> physical cores
        num_cores = psutil.cpu_count(logical=False)
        
        # Get the amount of RAM available in GB
        total_memory = psutil.virtual_memory().total / (1024**3)
        
        return num_cores, total_memory