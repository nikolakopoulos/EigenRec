EIGENREC Computes the Recommendation Matrix
   for the Cosine Inter-Item Proximity matrix
   INPUTS:
   R: Is the user x item Rating matrix
   f: Is the number of latent factors
   d: Is the scaling parameter
   OUTPUTS:
   Pi: The recommendation matrix (users x items). The rows of this matrix
       contain the recommendation vectors for every user in the system.
   dt: Time needed for Building the Latent Model.
   
  EXAMPLE:
  [ Pi, dt ] = EigenRec( R, 150, 0.2 );
    
  If you use our code please cite our paper

  "EIGENREC: Generalizing PureSVD for Effective and Efficient Top-N Recommendations" 

 -------------------------------------------
 Primary developers (in alphabetical order):
 -------------------------------------------
 1) Vassilis Kalantzis,    kalantzi@cs.umn.edu
 2) Athanasios N. Nikolakopoulos, anikolak@umn.edu
 For any questions or bug reports, please contact
 one of the developers listed above.
 
 Copyright (C) 2017 Vassilis Kalantzis, Athanasios N. Nikolakopoulos
 Release 1.0
