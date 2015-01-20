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
