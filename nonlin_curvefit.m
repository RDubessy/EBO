function p = nonlin_curvefit (f, p0, x, y)
% NONLIN_CURVEFIT Redefinition for compatibility with Octave.
p = lsqcurvefit (f,p0,x,y);
end