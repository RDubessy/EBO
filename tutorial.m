%% Optical Bloch Equation solver tutorial
% This short tutorial shows how to setup and use the OBE solver to model
% the interaction of atoms with a laser field.
% The aim of this tutorial is to show how the OBE toolbox can be used to go
% beyond simple two level approximations.
%
% First we clear the prompt, display and variables
clear all;
close all;
clc;
addpath ('./libebo');   %add the library path to matlab path
%% System definition
% First, we need to define the atomic parameters, here for a two level atom,
% and the magnetic field value $B_z$. This magnetic field defines the
% quantization axis.
Atom = getAtom('2Levels');
Bz = 1;
%%
% Then we call the initSystem() function that effectively initialize the
% matrices used to describe the system evolution.
[A,B,Bc,C,rho0,pmask,cmask,lmask,dab] = initSystem (Atom,Bz);
%% 
% The return values are:
%
% * A : a matrix describing the system free evolution,
% * B and Bc : two matrices describing the coupling with laser fields,
% * C : a matrix accounting for the laser detunings,
% * rho0 : a vector that contains the equilibrium system state with 
% population distributed over the ground states,
% * pmask, cmask and lmask : some logical masks that we will use later on,
% * dab : the matrix of Clebsch Gordan coefficients coupling the different
% states.
%
%%
% Then we define the laser parameters: namely its intensity and
% polarization. Using our notations, the laser intensity is defined in
% terms of Rabi pulsation $\Omega$, with the convention that:
%
% $$ \left(\frac{\Omega}{\Gamma}\right)^2 = \frac{I}{I_{sat}} $$
%
% where $I$ is the laser intensity on the atoms, $\Gamma$ is the linewidth
% of the transition and $I_{sat}$ is the saturation intensity of the closed
% cycle transition $\left|F=m_F=2\right>$ to
% $\left|F^\prime=m_F^\prime=3\right>$ (for Sodium atoms).
% The polarization is defined with respect to the quantization axis by a
% unit three component vector
% $\overrightarrow{\epsilon}=(\epsilon_{1},\epsilon_0,\epsilon_{-1})$. Each
% component refers to a polarization with respect to the quantization axis:
%
% * $\epsilon_{\pm1}$ are the $\sigma^{\pm}$ components,
% * $\epsilon_0$ is the $\pi$ component.
%
% With our notation [0,0,1] is the $\sigma^-$ polarization, [1,0,0] the
% $\sigma^+$ polarization and [0,1,0] the $\pi$ polarization.
%
% Here for a two level atom the polarization does not matter but for
% pratical purposes the only allowed polarization is $\sigma^+$.
s0 = 10;
OmegaPump = sqrt(s0) * Atom.Gamma * [1,0,0];
%%
% Finally we define the reference frequency of the lasers, here we assume
% that it is given by the $\left|F=2\right>$ to $\left|F^\prime=3\right>$
% transition.
deltaPump = Atom.excited.Energy(1) - Atom.ground.Energy(1);
%%
% In order to "see" the optical resonance we will span the laser frequency
% around this reference frequency. As $\Gamma$ is the "natural" frequency
% scale, we specify the detunings in units of $\Gamma$:
detunings = Atom.Gamma * linspace (-5,5,51);
%%
% We then compute the stationary state spectrum
[myA,myB,~] = rateMatrices (Atom,A,B,Bc,C,OmegaPump);
pop = zeros (sum (pmask),length (detunings));
for i=1:length(detunings)
    omega = Atom.omega0 + deltaPump + detunings(i);
    [M,Tp1,Tm1] = rateMatrix1Beam (myA,myB,omega,pmask,cmask,lmask);
    sol = null (M);
    pop(:,i) = sol / sum(sol);
end
figure(1);clf;
model = @(x,p) 0.5*p(1)./(1+p(1)+4*(x/p(2)).^2);
box on;
hold on;
plot (detunings,pop(1,:),'b+');
plot (detunings,pop(2,:),'r+');
plot (detunings,model(detunings,[s0,Atom.Gamma]),'k-');
hold off;
xlabel ('Detunings [MHz]');
ylabel ('Population');
legend ('Fundamental','Excited','Analytical two level model');
%%
% We notice an excellent agreement between the numerically computed
% stationary state and the exact theoretical solution.
% 
%% Dynamics
M = A + C * deltaPump + myB;
times = linspace (0,0.3,301);
dt = times(2)-times(1);
U = expm(2*pi*dt*M);
rho = zeros(length(rho0),length(times));
rho(:,1) = rho0;
omega = Atom.omega0 + deltaPump;
[M,Tp1,Tm1] = rateMatrix1Beam (myA,myB,omega,pmask,cmask,lmask);
Upop = expm(2*pi*dt*M);
pop = zeros(sum(pmask),length(times));
pop(:,1) = rho0(pmask);
for i=2:length(times)
    rho(:,i) = U * rho(:,i-1);
    pop(:,i) = Upop * pop(:,i-1);
end
model = @(x,p) 0.5*p(1)/(1+p(1))*(1-exp(-3*x*p(2)/4).*(...
    cos(x*p(2)/4*sqrt(8*p(1)-1))+3/sqrt(8*p(1)-1)*sin(x*p(2)/4*sqrt(8*p(1)-1))...
    ));
figure(2);clf;
data = real(rho(pmask,:));
hold on;
plot(times,data(2,:),'r+');
plot(times,pop(2,:),'m-');
plot([min(times) max(times)],0.5*s0/(1+s0)*[1 1],'k-');
plot(times,model(times,[s0 2*pi*Atom.Gamma]));
hold off;