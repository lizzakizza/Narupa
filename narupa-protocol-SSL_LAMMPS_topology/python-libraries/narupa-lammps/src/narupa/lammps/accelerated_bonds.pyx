# cython: boundscheck=False
# cython: cdivision=True
# cython: initializedcheck=False
# cython: cdivision_warnings=False
# cython: wraparound=False
# cython: binding=False
# cython: initializedcheck=False
# cython: nonecheck=False
# cython: overflowcheck=False
# cython: profile=False
# cython: linetrace=False
# distutils: define_macros=CYTHON_TRACE_NOGIL=0

from __future__ import print_function
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
import cython
cimport cython

from libc.math cimport fabs
from libc.stdlib cimport abs as c_abs
from cython.parallel import prange

cdef class fast_bonds:
    # Assign variable types so the ar available to the rest of the space.
    cdef long num_bonds
    cdef double cutoff
    cdef long[:,::1] adjusted_view
    cdef long[:,::1] starting_bonds
    cdef double[::1] *positions_ptr
    cdef long [::1] starting_bonds1_view
    cdef long [::1] starting_bonds2_view


    def __cinit__(self, long[:,::1] starting_bonds, long num_bonds, double cutoff ) : #, adjusted_view_ptr):
        self.num_bonds = num_bonds
        # Create contiguous memoryview for the new bond_list
        adjusted_bonds = np.empty([num_bonds, 2], dtype=int)
        cdef long [:,::1] adjusted_view = adjusted_bonds
        # Access adjusted bond list to flat pointer for speed

        # Convert access the pointer directly for the positions.
        cdef double [:,::1]  distances = np.empty([num_bonds, 3], dtype=np.double)

        # Split the topology into two contiguous arrays, this is clear during the main loop and avoids
        # two indexing additions per loop
        cdef long [::1] starting_bonds1_view =  np.ascontiguousarray(starting_bonds[:,0], dtype=int)
        cdef long [::1] starting_bonds2_view =  np.ascontiguousarray(starting_bonds[:,1], dtype=int)

        # Add to namespace for use in workhorse loop
        self.adjusted_view = adjusted_view
        self.starting_bonds = starting_bonds
        self.starting_bonds1_view = starting_bonds1_view
        self.starting_bonds2_view = starting_bonds2_view
        self.cutoff = cutoff

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    @cython.cdivision(True)
    cpdef extract_valid_bonds_c_no_gil(self, positions):
        cdef long x, y, z = 0
        cdef long new_number_of_bonds = 0
        cdef bint keep_pair
        cdef double distance #  = self.distances_view
        cdef long *starting_bonds = &(self.starting_bonds[0][0])
        cdef long num_bonds = self.num_bonds
        # Because the pointer address in lammps changes between MD steps we have to update it constantly.
        # Unfortunately this causes a factor of 10x slowdown compared to a fixed pointer.

        cdef double *positions_ptr = self.get_double_array_ptr(positions)
        cdef long *new_bond_list = &(self.adjusted_view[0][0])
        cdef long *starting_bonds1 = &(self.starting_bonds1_view[0])
        cdef long *starting_bonds2 = &(self.starting_bonds2_view[0])

        # Switching this to the nogil here could result in significant speedups with openMP, however
        # it requires reworking the code to avoid the z iterator and for the moment is not an issue.
        for x in range(num_bonds):
            # If any of the three distances evaluates to false eliminate it
            keep_pair = True
            for y in range(3):
               #take the absolute distance between each atom pair, this is slightly faster than taking the square
               distance  = fabs(positions_ptr[starting_bonds1[x]*3+y] - positions_ptr[starting_bonds2[x]*3+y])

               # if any bonds are longer than the required distances then skip
               # this is the part that slows down the code.
               if distance > self.cutoff:
                    keep_pair = False
                    break # Adds a tiny bit of speedup
            if keep_pair :
                new_bond_list[z]   = starting_bonds1[x]
                new_bond_list[z+1] = starting_bonds2[x]
                z += 2

        new_number_of_bonds = z/2
        return self.adjusted_view[0:new_number_of_bonds,:]

    cdef double * get_double_array_ptr(self, double_array):
        """Return the raw C pointer to the passed double array which was
           created using create_double_array
        """
        cdef double [::1] a = double_array
        return &(a[0])
