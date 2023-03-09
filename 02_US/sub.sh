
#To run the simulation
gmx_mpi grompp -f nvt.mdp -c struc.gro -r struc.gro -n index.ndx -p topol.top  -o  us.tpr 
mpirun -np 1 gmx_mpi mdrun -v -deffnm us -nb gpu -pme gpu  -bonded gpu  -plumed plumed.dat




