function [Energy,Pg,Pe] = atomicHamiltonian (Atom,Hhfs,Iz,Jz,ng,ne,Bz)
% ATOMICHAMILTONIAN Compute the hamiltonien in the presence of an external
% magnetic field BZ.
%
%   See also HYPERFINEHAMILTONIAN
muB = 1.399624604;  %MHz/Gauss
ntot = ng + ne;
Energy = zeros (ntot, 1);
gJpIz = Atom.gI*Iz +...
blkdiag(Atom.ground.gJ*eye(ng),Atom.excited.gJ*eye(ne)).*Jz;
Htot = Hhfs + muB * gJpIz * Bz;
[Pg, energy] = eig(Htot(1:ng,1:ng));
Energy(1:ng) = diag (energy);
[Pe, energy] = eig(Htot(ng+1:end,ng+1:end));
Energy(ng+1:ntot) = diag (energy);
end
