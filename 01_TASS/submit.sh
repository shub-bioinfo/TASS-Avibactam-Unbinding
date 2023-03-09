cd $PBS_O_WORKDIR
module load codes/intel/gmx-2018.8
module load codes/intel/plumed2-intel
export PLUMED_USE_LEPTON=yes

gmx_mpi grompp -f nvt-1.mdp  -c  struc-1.gro -p topol.top -n index.ndx -o tass.tpr
mpirun -np 1 gmx_mpi mdrun -v  -deffnm tass  -plumed plumed.dat -nb gpu   -pme gpu 

