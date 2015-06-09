function Atom = getAtom (name)
% GETATOM Return an ATOM structure containing the relevant atomic
% parameters for a particular fine structure line.
%   NAME is a string containing a valid atom name. For now the only 
%   authorized NAME value is 'Sodium'.
if (nargin ~= 1)
    error ('Invalid call to getAtom (name)');
end
Atom.name = name;
%Sodium properties (COMMON)
Atom.m = 0.381754035e-25;               %kg
Atom.k = 1./0.589;                      %MHz / [m/s]
Atom.Gamma = 9.7946;                    %MHz
Atom.omega0 = 508.848716e6;             %MHz
if (strcmp(name,'2Levels'))
    %Groundstate
    Atom.ground.Energy = [0];     %MHz
    Atom.ground.J = 0;
    %Excited state
    Atom.excited.Energy = [0]; %MHz
    Atom.excited.J = 1./2;
elseif (strcmp(name,'3Levels_LType'))
    %Groundstate
    Atom.ground.Energy = [-850,850];     %MHz
    Atom.ground.J = 0;
    %Excited state
    Atom.excited.Energy = [0]; %MHz
    Atom.excited.J = 0;
elseif (strcmp(name,'Sodium'))
    %Sodium D_2 line
    Atom.I = 3./2;
    Atom.gI = -0.00080461080;
    %Groundstate
    Atom.ground.J = 1./2;
    Atom.ground.F = [1,2];
    Atom.ground.Energy = [-1107.266, 664.3598];     %MHz
    Atom.ground.dip = [885.81306440,0,0];           %MHz
    Atom.ground.gJ = 2.00229600;
    %Excited state
    Atom.excited.J = 3./2;
    Atom.excited.F = [0,1,2,3];
    Atom.excited.Energy = [-66.097,-50.288,-15.944,42.382]; %MHz
    Atom.excited.dip = [18.534,2.724,0];                    %MHz
    Atom.excited.gJ = 1.33420;
else
    error('Unknown atom!');
end
end
