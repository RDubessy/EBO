function [rho, Mp] = stationaryState (M, mask)
% STATIONARYSTATE Compute the system stationary state for the given
% Liouville operator M.
%   The stationary state is definied by the matrix equation M*rho=0. As M 
%   is not invertible, there is a non trivial solution. Basically this is
%   due to the fact that the sum of the density matrix elements is constant
%   (equal to one). We use this to construct a modified Liouville operator
%   and the associated source term.
rho = zeros (length(M),1);
Mp = M(2:end,2:end);
Sp = zeros (length(Mp),1);
maskp = mask(2:end);
for i = 2:length(M)
    if ((M(i,1) ~= 0))
        Sp(i-1,1) = M(i,1);
        for j=maskp
            Mp(i-1,j-1) = Mp(i-1,j-1) - M(i,1);
        end
    end
end
rho(2:end,1) = -Mp\Sp;
rho(1,1) = 1-sum(rho(mask));
end
