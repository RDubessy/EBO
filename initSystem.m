function [A,B,Bc,C,rho0,mask] = initSystem (Atom, Bz)
% INITSYSTEM Initialize matrices used to describe the Liouvilian evolution 
% of an ATOM density matrix in an external magnetic field BZ (measured in 
% Gauss).
%   A contains the decay rates (populations and coherences), B and BC the 
%   coupling terms due to the laser and C the free evolution terms related 
%   to the laser detuning. RHO0 is an initial atomic density matrix with 
%   population equidistributed over the fundamental states and MASK is a 
%   convenient way to extract the populations: RHO0(MASK) contains only the
%   term describing the populations (diagonal of the density matrix).
%   Calling this function is only an intermediate step before calling the
%   COMPUTEMATRIX function.
%   Note that the density matrix elements are defined in the exact system 
%   basis imposed by the external bias field.
%
%   See also ATOMICHAMILTONIAN, COMPUTEDIPOLARCOUPLING
%muB = 1.399624604;  %MHz/Gauss
[Hhfs,Iz,Jz,ng,ne] = hyperfineHamiltonian (Atom);
[Energy,Pg,Pe] = atomicHamiltonian (Atom,Hhfs,Iz,Jz,ng,ne,Bz);
ntot = ng+ne;
A = zeros (ntot*ntot, ntot*ntot);
B = zeros (ntot*ntot, ntot*ntot, 3);
Bc = zeros (ntot*ntot, ntot*ntot, 3);
C = zeros (ntot*ntot, ntot*ntot);
rho0 = zeros (ntot*ntot, 1);
mask = zeros (1, ntot);
for i = 1:ntot
    mask(1,i) = i + (i-1)*ntot;
end
%Energy = zeros (ntot, 1);
%gJpIz = Atom.gI*Iz +...
%    blkdiag(Atom.ground.gJ*eye(ng),Atom.excited.gJ*eye(ne)).*Jz;
%Htot = Hhfs + muB * gJpIz * Bz;
%[Pg, energy] = eig(Htot(1:ng,1:ng));
%Energy(1:ng) = diag (energy);
%[Pe, energy] = eig(Htot(ng+1:end,ng+1:end));
%Energy(ng+1:ntot) = diag (energy);
dab = computeDipolarCoupling (Atom, ng, ne, Iz, Jz, Pg, Pe);
%Populations (fundamental) {{{
for i=1:ng
    %Source term
    rho0((i-1)*ntot+i,1) = 1.0 / ng;
    %Laser coupling to excited states
    for j=ng+1:ntot
        dip = dab(i,j-ng,:);
        Bc((i-1)*ntot+i,(j-1)*ntot+i,:) = 1i*dip/2;
        B((i-1)*ntot+i,(i-1)*ntot+j,:) = -1i*dip/2;
    end
    %Free evolution
    A((i-1)*ntot+i,(i-1)*ntot+i) = 0;
    for j=ng+1:ntot
        dip = reshape(dab(i,j-ng,:),1,3);
        ratio = (2*Atom.excited.J+1)/(2*Atom.ground.J+1)*(dip*dip');
        A((i-1)*ntot+i,(j-1)*ntot+j) = Atom.Gamma * ratio;
    end
end
% }}}
%Coherences (fundamental) {{{
for i=1:ng
    for ip=i+1:ng
        %Free evolution
        A((i-1)*ntot+ip,(i-1)*ntot+ip) = -1i*(Energy(i)-Energy(ip));
        A((ip-1)*ntot+i,(ip-1)*ntot+i) = -1i*(Energy(ip)-Energy(i));
        %Laser coupling to excited state
        for j=ng+1:ntot
            dip = dab(i,j-ng,:);
            Bc((i-1)*ntot+ip,(j-1)*ntot+ip,:) = 1i*dip/2;
            B((ip-1)*ntot+i,(ip-1)*ntot+j,:) = -1i*dip/2;
            dip = dab(ip,j-ng,:);
            B((i-1)*ntot+ip,(i-1)*ntot+j,:) = -1i*dip/2;
            Bc((ip-1)*ntot+i,(j-1)*ntot+i,:) = 1i*dip/2;
        end
    end
end
% }}} 
%Populations (excited states) {{{
for j=ng+1:ntot
    %Free evolution
    A((j-1)*ntot+j,(j-1)*ntot+j) = -Atom.Gamma;
    %Laser coupling to fundamental
    for i=1:ng
        dip = dab(i,j-ng,:);
        B((j-1)*ntot+j,(i-1)*ntot+j,:) = 1i*dip/2;
        Bc((j-1)*ntot+j,(j-1)*ntot+i,:) = -1i*dip/2;
    end
end
% }}}
%Coherences (excited states) {{{
for j=ng+1:ntot
    for jp=j+1:ntot
        %Free evolution
        A((j-1)*ntot+jp,(j-1)*ntot+jp) = -1i*(Energy(j)-Energy(jp))-Atom.Gamma;
        A((jp-1)*ntot+j,(jp-1)*ntot+j) = -1i*(Energy(jp)-Energy(j))-Atom.Gamma;
        %Laser coupling to fundamental
        for i=1:ng
            dip = dab(i,j-ng,:);
            B((j-1)*ntot+jp,(i-1)*ntot+jp,:) = 1i*dip/2;
            Bc((jp-1)*ntot+j,(jp-1)*ntot+i,:) = -1i*dip/2;
            dip = dab(i,jp-ng,:);
            Bc((j-1)*ntot+jp,(j-1)*ntot+i,:) = -1i*dip/2;
            B((jp-1)*ntot+j,(i-1)*ntot+j,:) = 1i*dip/2;
        end
    end
end
% }}}
%Coherences (excited - fundamental) {{{
for j=ng+1:ntot
    for i=1:ng
        %Free evolution
        A((j-1)*ntot+i,(j-1)*ntot+i) = -1i*(Energy(j)-Energy(i))-Atom.Gamma/2;
        A((i-1)*ntot+j,(i-1)*ntot+j) = -1i*(Energy(i)-Energy(j))-Atom.Gamma/2;
        %Detunings (one photon)
        C((j-1)*ntot+i,(j-1)*ntot+i) = 1i;
        C((i-1)*ntot+j,(i-1)*ntot+j) = -1i;
        %Laser coupling
        for ip=1:ng
            dip = dab(ip,j-ng,:);
            B((j-1)*ntot+i,(ip-1)*ntot+i,:) = 1i*dip/2;
            Bc((i-1)*ntot+j,(i-1)*ntot+ip,:) = -1i*dip/2;
        end
        for jp=ng+1:ntot
            dip = dab(i,jp-ng,:);
            B((j-1)*ntot+i,(j-1)*ntot+jp,:) = -1i*dip/2;
            Bc((i-1)*ntot+j,(jp-1)*ntot+j,:) = 1i*dip/2;
        end
    end
end
% }}}
end
