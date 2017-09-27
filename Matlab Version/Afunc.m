function [ Y ] = Afunc(x,W)
%AFUNC This Function Computes the MatrixVector product with the Inter-Item
%Proximity matrix that can be expressed in terms of sparse and/or low-rank
%components. Here as an example we give the simple Cosine Inter-Proximity
%Matrix. vector x is fed by the eigs function. (see EIGENREC.m) 

Y = (W'*(W*x));  % General Cosine Inter-Item Proximity Matrix. 
                 % For Simple PureSVD choose (W==R)
end

