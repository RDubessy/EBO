function [myA,myB,myBp] = rateMatrices (varargin)
if (nargin == 6 || nargin == 7)
    Atom = varargin{1};
    A = varargin{2};
    B = varargin{3};
    Bc = varargin{4};
    C = varargin{5};
    Omega = varargin{6};
    myA = A - Atom.omega0 * C;
    myB = zeros(size(A));
    myBp = zeros(size(A));
    for i = 1:3
        myB = myB + Omega(i) * B(:,:,i) + conj(Omega(i)) * Bc(:,:,i);
    end
    if (nargin == 7)
        Omegap = varargin{7};
        for i = 1:3
            myBp = myBp + Omegap(i) * B(:,:,i) + conj(Omegap(i)) * Bc(:,:,i);
        end
    end
end
end
