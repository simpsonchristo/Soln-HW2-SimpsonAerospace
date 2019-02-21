function [E] = keplerseqn(e, n, t, M0)
%KEPLERSEQN Solve eccentric anomaly for given mean anomaly and eccentricity
%Mean anomaly is determined by the mean motion, time, and time of perifocus
%passage. Employing a root-finding method, one can determine the eccentric
%anomaly when provided eccentricity. ToDo: Verify Newton-Raphson method is working properly
M = n*t + M0;
g = 1.0e+10;
maxiter = 1000;
i = 0;
%initial guess
E = M;
while abs(g)>1.0e-6
    gprime = 1 - e*cos(E);
    g = E - (e*sin(E)) - M;
    if abs(g)>1e-6
        E = E - (g/gprime);
    end
    i= i+1;
    if i>maxiter
        break
    end
end
end