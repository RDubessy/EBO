function [M,Tp1,Tm1] = rateMatrix1Beam (A,B,omega,pmask,cmask,lmask)
Id = eye(size(A));
Tp1 = -(A(cmask,cmask)-1i*omega*Id(cmask,cmask))\B(cmask,pmask);
Tm1 = -(A(cmask,cmask)+1i*omega*Id(cmask,cmask))\B(cmask,pmask);
M = A(pmask,pmask)+B(pmask,cmask)*(Tp1+Tm1);
end
