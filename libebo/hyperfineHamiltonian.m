function [Hhfs, Iz, Jz, ng, ne] = hyperfineHamiltonian (Atom)
% HYPERFINEHAMILTONIAN Initialize matrices used to describe the internal 
% state of an ATOM, in the |I,m_I,J,m_J> representation.
%   HHFS is the hyperfine hamiltonien, IZ and JZ are the matrices of the 
%   operators I and J along the quantization axis "z". NG is the number of 
%   ground states, NE is the number of excited states.
%
%   See also ANGULARMOMENTUM
I = Atom.I;
[Ip,Im,Iz] = angularMomentum (I);
%Ground state
J = Atom.ground.J;
[Jp,Jm,Jz] = angularMomentum (J);
Ahfs = Atom.ground.dip(1);
%Bhfs = Atom.ground.dip(2);         %Bhfs and Chfs are zero for the ground
%Chfs = Atom.ground.dip(3);         %state
IsJ = (kron(Ip,Jm)+kron(Im,Jp))/2+kron(Iz,Jz);
Hghfs = Ahfs * IsJ;
gIz = kron (Iz,eye(size(Jz)));
gJz = kron (eye(size(Iz)),Jz);
%Excited state
J = Atom.excited.J;
[Jp,Jm,Jz] = angularMomentum (J);
Ahfs = Atom.excited.dip(1);
Bhfs = Atom.excited.dip(2);
Chfs = Atom.excited.dip(3);
IsJ = (kron(Ip,Jm)+kron(Im,Jp))/2+kron(Iz,Jz);
Id = eye (size (IsJ));
Hehfs = Ahfs * IsJ +...
    Bhfs * (3*IsJ^2+3./2*IsJ-Id*I*(I+1)*J*(J+1))/(2*I*(2*I-1)*J*(2*J-1))+...
    Chfs * (10*IsJ^3+20*IsJ^2+2*IsJ*(I*(I+1)+J*(J+1)+3)-Id*(...
    3*I*(I+1)*J*(J+1)+5*I*(I+1)*J*J+1))/(I*(I-1)*(2*I-1)*J*(J-1)*(2*J-1));
eIz = kron (Iz,eye(size(Jz)));
eJz = kron (eye(size(Iz)),Jz);
Hhfs = blkdiag (Hghfs, Hehfs);
Iz = blkdiag (gIz, eIz);
Jz = blkdiag (gJz, eJz);
ng = length (gIz);
ne = length (eIz);
end

function [Jp,Jm,Jz] = angularMomentum (J)
% ANGULARMOMENTUM Compute the three matrices representing an angular
% momentum operator J.
%   [JP,JM,JZ] = ANGULARMOMENTUM(J)
%   JZ is the angular momentum projection onto the quantization axis.
%   JM decreases by one the angular momentum projection.
%   JP increases by one the angular momentum projection.
n = 2*J+1;
Jz = zeros (n,n);
Jp = Jz;
Jm = Jz;
for i=1:n
    Jz(i,i) = -J+i-1;
end
for i=1:n-1
    m = Jz(i,i);
    Jp(i+1,i) = sqrt(J*(J+1)-m*(m+1));
    m = Jz(i+1,i+1);
    Jm(i,i+1) = sqrt(J*(J+1)-m*(m-1));
end
end