function [velocities] = velocityClasses (Atom,detuning)
    velocities = (-abs(detuning)-2500):100:(abs(detuning)+2500);
    span = 100 + floor(max(Atom.excited.Energy)-min(Atom.excited.Energy));
    step = 5;
    for Eg = Atom.ground.Energy
        zoom = -span:step:span;
        velocities = [velocities,zoom-detuning-Eg,zoom+detuning+Eg];
    end
    velocities = sort(unique(velocities))/Atom.k;
end