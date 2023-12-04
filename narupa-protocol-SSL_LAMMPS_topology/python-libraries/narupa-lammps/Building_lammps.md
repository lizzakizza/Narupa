# Installing Lammps and linking with narupa
The following guide shows how to install lammps on an ubuntu 20.04 cloud instance and make sure that narupa works with it.

This part covers setting up the system with compilers and anaconda so that it can build lammps and narupa
    
    1  sudo apt install build-essential
    2  wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
    3  sh Anaconda3-2020.02-Linux-x86_64.sh 
    4  conda activate base
    5  source .bashrc    
    
Installing narupa

    1  git clone https://gitlab.com/intangiblerealities/narupa-protocol.git
    2  cd narupa-protocol/
    3  git checkout feature/LAMMPS_topology
    4  conda install mpi4py
    5  ./compile.sh 

Optional: Check that narupa is working
   
    1  cd python-libraries/narupa-lammps/tests
    2  pytest test_lammps_converter.py 
    
Building lammps     
    
    1  git clone -b unstable https://github.com/lammps/lammps.git mylammps
    2  mkdir build
    3  cd build/
    4  cmake -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++  -D PKG_PYTHON=on -D BUILD_LIB=on -D BUILD_SHARED_LIBS=on -D LAMMPS_EXCEPTIONS=on -D BUILD_MPI=on -D PKG_USER-OMP=on -D PKG_MOLECULE=on -D PKG_CLASS2=on -D PKG_KSPACE=on -D PKG_USER-REAXC=yes -D PKG_USER-CGDNA=yes -D PYTHON_EXECUTABLE=~/anaconda3/envs/bin/python3.7 ../cmake
    5  make -j
    6  make install-python
    7  echo 'export PYTHONPATH=${PYTHONPATH}:${HOME}/mylammps/python/' >> ~/.bashrc 
    8  echo 'export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/mylammps/build' >> ~/.bashrc
    9  source .bashrc    

    
At this point you should see a lmp executable in this folder, you could add this executable to your path or move it to
the folder of an input file.   
   
The final stage is to open the firewall so that connections can be made. You may wish to limit the adresses that
are allowed though the firewall. 
   
    1 sudo iptables -I INPUT 1 -p tcp --dport 38801 -j ACCEPT; sudo iptables -I OUTPUT 1 -p tcp --dport 38801 -j ACCEPT 

How to modify lammps files so that they run.
Add the following line to your input somewhere before the lammps command.

        python         post_force_callback here """
        from narupa.lammps import LammpsImd
        from lammps import lammps
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        h = LammpsImd()
        
        def post_force_callback(lmp,v):
          h.lammps_hook(lmp,comm)
        """

and add this after any other fixes in the the input file

        fix     3 all python/invoke 1 post_force post_force_callback

this means that after every times step narupa can gather and scatter the positions and forces.


Running lammps with MPI over 32 cores:
     
     1  mpirun -n 32 ./lmp < FAU-acetone-suprcell.in 
     
Topologies:

Topologies are auto-loaded from any .data file in the same folder as the input, or a user can specify the topology 
by adding the following line in the python input (customised ti the filename of your data file):

    h.data_file_for_bonds="./data.rdx"      

Which in the context of the full file

    python         post_force_callback here """
    from narupa.lammps import LammpsImd
    from lammps import lammps
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    h = LammpsImd()
    h.data_file_for_bonds="./data.rdx"
    
    def post_force_callback(lmp,v):
      h.lammps_hook(lmp,comm)
    """

