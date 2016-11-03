# WeakLMoka
until now the parameter file weaklensingMOKA.ini
has 17 parameters .. each of them is followed in order by a whole string line

!...
 1. omega matter
!...
 2. omega lambda 
!...
 3. hubble constant in unit of 100 
!...
 4. dark energy equation of state parameter
!...
 5. source redshift
!...
 6. radius at which cut the density profile in 1/Rvir
!...
 7. path of files
!...
 8. total number of planes in which the catalogues are divided no only up to the source redshift
!...
 9. simulation type name
!...
10. planelist file produced by MapSim (ask cgiocoli@gmail.com for it eventually!)
!...
11. number of pixels in x
!...
12. number of pixels in y
!...
13. size of the field of view in x (degree)
!...
14. size of the field of view in y (degree)
!...
15. file with lm s relation to be adopted (if not uses Neto et al. evolved as Bullock et al.)
!...
16. log-normal scatter in the c-M relation if a c-M model is used
!...
17. if equal to vir use FOF haloes and masses if different 200rhoc haloes and masses

  - - - ---  WARNING --- - - - 
  (1)
  /**
   *  check on cutR (parameter 6.):  
   *  this value could be positive or negative ... not null
   *  if positive the integral in z of kappa is extended as the profile
   *  on the plane of the sky
   *  if negative on the plane of the sky is threated as it was positive
   *  while alone of the line-of-sight the integral of the density is 
   *  extended up to infinity
   **/
   (2)
   /**
    *   
    *
    *
    *
    *
    *
    **/
  