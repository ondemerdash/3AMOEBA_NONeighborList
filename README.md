3AMOEBA Src, Parallelized with MPI, Polarization Implemented Without Neighbor Lists.
Below is an explanation of the important modifications of and additions to the standard Tinker code:

1. analyze_parallel_bcast_evenload.f :  Performs multiple iterations of an energy/gradient/virial calculation
   with the polarization calculation polarized with MPI.  This is essentially a means of testing the parallel 
   implementation, as ultimately this needs to be implemented in Tinker's standard molecular dynamics routine,
   dynamic.f (Work in progress).  

2. innerloop_empole1a_3b.f :  Contains the subroutines Innerloop2 and Innerloop3 which respectively calculate 
   the 2- and 3-body polarization energy, gradient, and virial for a chunk of work. (Work is trivially divided 
   among the MPI tasks in analyze_parallel_bcast_evenload.f)

3. empole1c_1a_3bi_reform.f : Contains the subroutine empole1a_3b_Polar, which calculates the polarization energy,
   gradient, and virial using standard Coulombic (as opposed to Ewald) electrostatics with an infinite cutoff.   
