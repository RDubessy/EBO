function [A,B,Bc,C,rho0,pmask,cmask,lmask,dab] = initSystem (Atom, Bz)
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
%   LMASK is an auxiliary logical mask.
%
%   See also ATOMICHAMILTONIAN, COMPUTEDIPOLARCOUPLING
if (strcmp(Atom.name,'2Levels'))
    ng = 1;
    ne = 1;
    Energy = [Atom.ground.Energy(1),Atom.excited.Energy(1)];
    dab = zeros (1,1,3);
    dab(1,1,1) = 1/sqrt(2);
elseif (strcmp(Atom.name,'3Levels_LType'))
    ng = 2;
    ne = 1;
    Energy = [Atom.ground.Energy(1),Atom.ground.Energy(2),Atom.excited.Energy(1)];
    dab = zeros (2,1,3);
    dab(:,1,1) = [1;1]/sqrt(2);
else
    [Hhfs,Iz,Jz,ng,ne] = hyperfineHamiltonian (Atom);
    [Energy,Pg,Pe] = atomicHamiltonian (Atom,Hhfs,Iz,Jz,ng,ne,Bz);
    dab = computeDipolarCoupling (Atom, ng, ne, Iz, Jz, Pg, Pe);
end
%Initialize return values
ntot = ng+ne;
A = zeros (ntot*ntot, ntot*ntot);
B = zeros (ntot*ntot, ntot*ntot, 3);
Bc = zeros (ntot*ntot, ntot*ntot, 3);
C = zeros (ntot*ntot, ntot*ntot);
rho0 = zeros (ntot*ntot, 1);
pmask = false(1, ntot*ntot);
cmask = false(1, ntot*ntot);
lmask = false(1, ntot*ntot);
for i = 1:ntot
    for j = 1:ntot
        k = j + (i-1) * ntot;
        if (i == j)
            pmask(k) = true;
        elseif ((i<=ng && j<=ng)||(i>ng && j>ng))
            lmask(k) = true;
        else
            cmask(k) = true;
        end
    end
end
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
