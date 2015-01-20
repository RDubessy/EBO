function p = nonlin_curvefit (f, p0, x, y)
% NONLIN_CURVEFIT Redefinition for compatibility with Octave.
options = optimset('Display','off');
p = lsqcurvefit (f,p0,x,y,[],[],options);
end