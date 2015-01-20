function [M] = computeMatrix (A,B,Bc,C,Omega,delta)
% COMPUTEMATRIX Compute the Liouville operator governing the density matrix
% evolution.
    M = A+C*delta;
    Omegac = conj(Omega);
    for i=1:3
        M(:,:) = M(:,:) + Omega(i) * B(:,:,i) + Omegac(i) * Bc(:,:,i);
    end
end
