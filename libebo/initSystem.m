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

function dab = computeDipolarCoupling (Atom, ng, ne, Iz, Jz, Pg, Pe)
% COMPUTEDIPOLARCOUPLING Compute the reduced dipolar electric coupling
% coefficients of a fine structure line.
%   This gives access to the so called Clebsch-Gordan coefficient between
%   particular sub-levels.
%   The return value DAB is the NG times NE coupling matrix between the 
%   ground and excited states of the ATOM. IZ and JZ are the projection of 
%   the nuclear spin and the total angular momentum operators and PG and PE
%   are the actual spin eigenstates of the ATOM.
%   Note that this function is only called in the INITSYSTEM function and 
%   should not be directly called.
%
%   See also INITSYSTEM
dab = zeros (ng, ne,3);
for i=1:ng
    for j=1:ne
        tmp = [0,0,0];
        for ii=1:ng
            for jj=1:ne
                mI = Iz(ii,ii);
                mIp = Iz(jj+ng,jj+ng);
                if (mI == mIp)
                    mJ = Jz(ii,ii);
                    mJp = Jz(jj+ng,jj+ng);
                    scale = Pg(ii,i) * Pe(jj,j) *...
                        (-1)^(Atom.excited.J-1+mJ)*...
                        sqrt(2*Atom.ground.J+1);
                    tmp(1) = tmp(1) + scale*w3j(Atom.excited.J,mJp,1,-1,...
                                Atom.ground.J,-mJ);
                    tmp(2) = tmp(2) + scale*w3j(Atom.excited.J,mJp,1,0,...
                                Atom.ground.J,-mJ);
                    tmp(3) = tmp(3) + scale*w3j(Atom.excited.J,mJp,1,1,...
                                Atom.ground.J,-mJ);
                end
            end
        end
        dab(i,j,:) = tmp;
    end
end
end

function w=w3j(j1,m1,j2,m2,j3,m3)
% W = W3J(j1,m1,j2,m2,j3,m3)
%
% j1,j2,j3 must satisfy triangle condition |j2 - j3|<=j1<=j2 + j3
%
% Wigner 3-j symbol is evaluated using equation found in 'Angular 
% Momentum: An Illustrated guide to Rotational Symmetries for
% Physical Systems', W. J. Thompson
%
%J. Pritchard Durham University 2009

%Check Triangular relation
if(((j3<abs(j1-j2))||(j3>(j1+j2))))
    error('Addition of angular momentum requires triangle relation\n\t|j1-j2|<=j3<=j1+j2');
%Evaluate w3j
else
    if(m1+m2+m3~=0)
        w=0;
    elseif(j2==0)
        w=1;
    else
        w=(-1)^(j1-j2-m3)*...
        exp(0.5*(lgf(j3+j1-j2)+lgf(j3-j1+j2)+lgf(j1+j2-j3)...
        +lgf(j3-m3)+lgf(j3+m3)- lgf(j1+j2+j3+1)...
        -lgf(j1-m1)-lgf(j1+m1)-lgf(j2-m2)-lgf(j2+m2)))...
        *ksum(j1,m1,j2,m2,j3,m3);
    end
end
end

%Summation performed for all values of k which give non-negative fs
function s=ksum(j1,m1,j2,m2,j3,m3)
s=0;
kmin=max([0,-j1+j2-m3]);
kmax=min([j3-j1+j2,j3-m3]);
for k=kmin:kmax
    s = s +(-1)^(k+j2+m2)*exp(lgf(j2+j3+m1-k)+lgf(j1-m1+k)...
        -lgf(k)-lgf(j3-j1+j2-k)-lgf(j3-m3-k)-lgf(k+j1-j2+m3));
end
end

%Stirlings approximation ln(n!)=n*ln(n)-n+0.5*ln(2*pi*n)
function y=lgf(x)
if(x<170)
    y=log(factorial(x));
else
    y=x*log(x)-x+0.5*log(2*pi*x)+1/(12*x)-1/(360*x^3)+1/(1260*x^5)...
        -1/(1680*x^7)+1/(1188*x^9);
end
end