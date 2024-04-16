import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


color_of_molecules = {'H': 'red', # hydrogen
                      'C': 'black', # carbon
                      'O': 'blue', # oxygen
                      'N': 'green' # nitrogen
                      } 

def graph_in_file(xyz_file):
    # Load data from file
    data = np.loadtxt(xyz_file, skiprows=2, dtype=str)  # Assuming the first two lines contain metadata

    # Extract x, y, z coordinates and atom type
    x = data[:, 1].astype(float)
    y = data[:, 2].astype(float)
    z = data[:, 3].astype(float)
    atom_type = data[:, 0]

    # Create a dictionary to map atom types to colors
    atom_color_mapping = {'H': 'red', 'C': 'black', 'O': 'blue', 'N': 'green'}

    # Get colors based on atom type
    colors = [atom_color_mapping[atom] for atom in atom_type]

    # Create a 3D figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the points in 3D with colors according to the type of atom
    scatter = ax.scatter(x, y, z, c=colors, s=50, label='Átomos')

    # Additional settings as needed
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Visualization of Atoms')

    # Show the figure
    plt.show()

graph_in_file('./atom3.xyz')