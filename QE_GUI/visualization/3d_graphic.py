import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

color_of_atoms = {'H': 'red',  # hydrogen
                  'C': 'black',  # carbon
                  'O': 'blue',  # oxygen
                  'N': 'green'  # nitrogen
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

    atoms_info = {}

    # Function to draw bonds between atoms
    def draw_bonds(ax, x, y, z, atom_type):
        #if atoms_info.get(atom_type, False) != '':
        #    atoms_info[atom_type] = atom_color_mapping.get(color_of_atoms, "gray")

        for i in range(len(x)):
            for j in range(i + 1, len(x)):
                # Calculate distance between atoms
                dist = np.sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2 + (z[i] - z[j]) ** 2)
                # Draw a line if the distance is less than a threshold
                if dist < 2.0:  # Adjust this threshold as needed
                    ax.plot([x[i], x[j]], [y[i], y[j]], [z[i], z[j]], color='gray', alpha=0.5)

    # Call the function to draw bonds
    draw_bonds(ax, x, y, z, atom_type)

    # Create custom legends for atom type and color
    custom_legends = []
    for atom, color in color_of_atoms.items():
        custom_legends.append(plt.Line2D([0], [0], marker='o', color='w', label=atom, markersize=10, markerfacecolor=color))

    # Additional settings as needed
    #ax.set_xlabel('X')
    #ax.set_ylabel('Y')
    #ax.set_zlabel('Z')
    # Get rid of colored axes planes
    # First remove fill
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    # Now set color to white (or whatever is "invisible")
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')
    # Bonus: To get rid of the grid as well:
    ax.grid(False)
    ax.set_title('3D Visualization of Atoms')

    # Show legends
    ax.legend(handles=custom_legends)

    # Show the figure
    plt.show()


graph_in_file('./atom3.xyz')
