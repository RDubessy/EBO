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
