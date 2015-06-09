function [M,T_p1_0,T_0_p1] = rateMatrix2Beam (A,B,omega,Bp,omegap,pmask,cmask,lmask)
Id = eye(size(A));
% coh <- pop
Gppp = (A(cmask,cmask)+1i*omegap*Id(cmask,cmask))\Bp(cmask,pmask);
Gm00 = (A(cmask,cmask)-1i*omega*Id(cmask,cmask))\B(cmask,pmask);
Gmpp = (A(cmask,cmask)-1i*omegap*Id(cmask,cmask))\Bp(cmask,pmask);
Gp00 = (A(cmask,cmask)+1i*omega*Id(cmask,cmask))\B(cmask,pmask);
% coh <- intra
Gpp0 = (A(cmask,cmask)+1i*omegap*Id(cmask,cmask))\B(cmask,lmask);
Gm0p = (A(cmask,cmask)-1i*omega*Id(cmask,cmask))\Bp(cmask,lmask);
Gmp0 = (A(cmask,cmask)-1i*omegap*Id(cmask,cmask))\B(cmask,lmask);
Gp0p = (A(cmask,cmask)+1i*omega*Id(cmask,cmask))\Bp(cmask,lmask);
% intra <- pop
T_p1_m1 = (A(lmask,lmask)-1i*(omega-omegap)*Id(lmask,lmask)-B(lmask,cmask)*Gpp0-Bp(lmask,cmask)*Gm0p)\(B(lmask,cmask)*Gppp+Bp(lmask,cmask)*Gm00);
T_m1_p1 = (A(lmask,lmask)+1i*(omega-omegap)*Id(lmask,lmask)-B(lmask,cmask)*Gmp0-Bp(lmask,cmask)*Gp0p)\(B(lmask,cmask)*Gmpp+Bp(lmask,cmask)*Gp00);
% coh <- pop
T_p1_0 = -Gm00-Gm0p*T_p1_m1;
T_m1_0 = -Gp00-Gp0p*T_m1_p1;
T_0_p1 = -Gmpp-Gmp0*T_m1_p1;
T_0_m1 = -Gppp-Gpp0*T_p1_m1;
% pop <- pop
M = A(pmask,pmask)+B(pmask,cmask)*(T_p1_0+T_m1_0)+Bp(pmask,cmask)*(T_0_p1+T_0_m1);
end
