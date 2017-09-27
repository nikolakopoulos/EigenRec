function [ Pi, dt ] = EigenRec( R, f, d )
%EIGENREC Computes the Recommendation Matrix
%   for the Cosine Inter-Item Proximity matrix
%   INPUTS:
%   R: Is the user x item Rating matrix
%   f: Is the number of latent factors
%   d: Is the scaling parameter
%   OUTPUTS:
%   Pi: The recommendation matrix (users x items). The rows of this matrix
%       contain the recommendation vectors for every user in the system.
%   dt: Time needed for Building the Latent Model.
%   
%  EXAMPLE:
%  [ Pi, dt ] = EigenRec( R, 150, 0.2 );
%    
%  If you use our code please cite our paper
%
%  "EIGENREC: Generalizing PureSVD for Effective and Efficient Top-N Recommendations" 
%
% -------------------------------------------
% Primary developers (in alphabetical order):
% -------------------------------------------
% 1) Vassilis Kalantzis,    kalantzi@cs.umn.edu
% 2) Athanasios N. Nikolakopoulos, anikolak@umn.edu
% For any questions or bug reports, please contact
% one of the developers listed above.
% 
% Copyright (C) 2017 Vassilis Kalantzis, Athanasios N. Nikolakopoulos
% Release 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Example using the Cosine Similarity Matrix (see our paper for details)
W = R*diag(diag(sparse(diag(sqrt(diag(R'*R)))).^(d-1))); % The inter-item similarity matrix A_cos = W'*W
[n,m] = size(W);
OPTS.issym = 1; OPTS.tol = 1e-8;
tic
[V,~] = eigs(@(x)Afunc(x,W),m,f,'LM',OPTS); % V is the reduced latent-model
dt = toc;
Pi = R*V*V';

% Similar constructions can be used for the other similarity matrices
% proposed in our paper. When we have the matrix expressed in terms of sparse components
% it suffices to provide a function that computes the MV product, similar
% to the Afunc we provide for the Cosine matrix.
% When the final inter-item proximity matrix is formed explicitly we run 
% [Q,~] = eigs(A,m,f,'LM',OPTS);
% instread. 
%
% Notice that the general diagonal scaling matrix is given by
% S = diag(diag(sparse(diag(sqrt(diag(R'*R)))).^d)); 

end

