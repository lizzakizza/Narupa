# Timing routine to test the speed of custom cython code generating dynamic bonds
import timeit

cython_100 = """
import numpy as np
import accelerated_bonds
import ctypes

numatoms = 100
positions = ((ctypes.c_double * (3 * numatoms))(*[float(x) for x in range(3*numatoms)]))
starting_bonds = np.random.randint(low=0, high=numatoms, size=(numatoms*3,2))
b = accelerated_bonds.fast_bonds(starting_bonds,  starting_bonds.shape[0], 10.0)
"""

cython_1000 = """
import numpy as np
import accelerated_bonds
import ctypes

numatoms = 1000
positions = ((ctypes.c_double * (3 * numatoms))(*[float(x) for x in range(3*numatoms)]))
starting_bonds = np.random.randint(low=0, high=numatoms, size=(numatoms*3,2))
b = accelerated_bonds.fast_bonds(starting_bonds,  starting_bonds.shape[0], 10.0)
"""

cython_10000 = """
import numpy as np
import accelerated_bonds
import ctypes

numatoms = 10000
positions = ((ctypes.c_double * (3 * numatoms))(*[float(x) for x in range(3*numatoms)]))
starting_bonds = np.random.randint(low=0, high=numatoms, size=(numatoms*3,2))
b = accelerated_bonds.fast_bonds(starting_bonds,  starting_bonds.shape[0], 10.0)
"""

cython_100000 = """
import numpy as np
import accelerated_bonds
import ctypes

numatoms = 100000
positions = ((ctypes.c_double * (3 * numatoms))(*[float(x) for x in range(3*numatoms)]))
starting_bonds = np.random.randint(low=0, high=numatoms, size=(numatoms*3,2))
b = accelerated_bonds.fast_bonds(starting_bonds,  starting_bonds.shape[0], 10.0)
"""

narupa_100 = """
import ctypes
import numpy as np

numatoms = 100
positions = ((ctypes.c_double * (3 * numatoms))(*[float(x) for x in range(3*numatoms)]))
positions_3n = np.ctypeslib.as_array(positions, shape=(numatoms * 3)).reshape(numatoms, 3)
bonds_list = np.random.randint(low=0, high=numatoms, size=(numatoms*3,2))
"""

narupa_1000 = """
import ctypes
import numpy as np

numatoms = 1000
positions = ((ctypes.c_double * (3 * numatoms))(*[float(x) for x in range(3*numatoms)]))
positions_3n = np.ctypeslib.as_array(positions, shape=(numatoms * 3)).reshape(numatoms, 3)
bonds_list = np.random.randint(low=0, high=numatoms, size=(numatoms*3,2))
"""

narupa_10000 = """
import ctypes
import numpy as np

numatoms = 10000
positions = ((ctypes.c_double * (3 * numatoms))(*[float(x) for x in range(3*numatoms)]))
positions_3n = np.ctypeslib.as_array(positions, shape=(numatoms * 3)).reshape(numatoms, 3)
bonds_list = np.random.randint(low=0, high=numatoms, size=(numatoms*3,2))
"""

narupa_100000 = """
import ctypes
import numpy as np

numatoms = 100000
positions = ((ctypes.c_double * (3 * numatoms))(*[float(x) for x in range(3*numatoms)]))
positions_3n = np.ctypeslib.as_array(positions, shape=(numatoms * 3)).reshape(numatoms, 3)
bonds_list = np.random.randint(low=0, high=numatoms, size=(numatoms*3,2))
"""

narupa_run = """

def reduce_logic(booly_boy):
    # Compare each row against the last and evaluate if the elements are true, this is v.fast
    output = booly_boy[0]  # separate out first element
    output = np.logical_and(output, booly_boy[1])
    output = np.logical_and(output, booly_boy[2])
    return output

atom1 = np.empty([len(bonds_list), 3])
atom2 = np.empty([len(bonds_list), 3])
# Extract the atoms index by bond pairs into separate lists
positions_3n.take(bonds_list[:, 0], axis=0, out=atom1)
positions_3n.take(bonds_list[:, 1], axis=0, out=atom2)
# Get the absolute distance between atoms along all axis
distances = np.abs((atom1 - atom2))

# Compare get a bool list of all distances for each axis
# The transpose makes the vector operation more efficient.
booly_boy = np.less(distances.T, 10.0)
# Reduce logic is faster than other approaches as it accumulates the bool along each axis
booly_boy2 = reduce_logic(booly_boy)
# Generate the new bond list based on the bool evaluate of the distances
updated_bond_list = np.compress(booly_boy2, bonds_list, axis=0)


"""

cython_run = """
adjusted_bonds = b.extract_valid_bonds_c_no_gil(positions)
"""

narupa1_100 = timeit.timeit(stmt=narupa_run, setup=narupa_100, number=100)
narupa1_1000 = timeit.timeit(stmt=narupa_run, setup=narupa_1000, number=100)
narupa1_10000 = timeit.timeit(stmt=narupa_run, setup=narupa_10000, number=100)
narupa1_100000 = timeit.timeit(stmt=narupa_run, setup=narupa_100000, number=100)

cython3_100 = timeit.timeit(stmt=cython_run, setup=cython_100, number=100)
cython3_1000 = timeit.timeit(stmt=cython_run, setup=cython_1000, number=100)
cython3_10000 = timeit.timeit(stmt=cython_run, setup=cython_10000, number=100)
cython3_100000 = timeit.timeit(stmt=cython_run, setup=cython_100000, number=100)


print("\n")
print("The following results are for a given number of atom and 3x the number of bonds")
print("\n")

print("Num Atoms :    {:d}      {:d}      {:d}   {:d}".format(100, 1000, 10000, 100000))
print("Numpy        :      {:.6f} {:.6f} {:.6f} {:.6f}".format(narupa1_100, narupa1_1000, narupa1_10000, narupa1_100000))
print("cython float :      {:.6f} {:.6f} {:.6f} {:.6f}".format(cython3_100, cython3_1000, cython3_10000, cython3_100000))

