# HPC_Final-Project
This repository is for Fnial Project of HPC course

This repository has three version code, explicit.c and implicit.c is the final version with HDF5 I/O to achieve program restart.
To run the code
### First
use "make implicit" or "make explicit" to compile the c++ code

### Second 
use "mpirun -np 'n' ./explicit" to run the code . you can change time step or mesh size with " -n 'n' " or " -dt 'n' ",where 'n' is integer number

### Third
If you run the code with cluster, use "bsub < FinalPro_jobscript" to submit job.
