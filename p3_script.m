clear;
clc;
%p3
%J2 propagation (perturbations) for GLONASS
mu = 3.9860044e+14; %m^3/s^2, Earth gravitational parameter
we = (2*pi()/86164);%rad/sec, Earth avg rotational rate
re = 6378137;%meters, spherical Earth radius
ag = 0;%Greenwich angle
J2 = 0.001082636;
C2_0 = -0.0004841695;%normalized value of harmonic coefficient

%solve for initial state
t = 0:60:(48*60*60);%seconds
oe0 = [25500000, 0.0015, 63.0, 300, 0.01, 0.0];%m, [], deg, deg, deg, deg
f0 = 0;
[X0] = oe2rv(oe0, f0);%provide intial state

%provide ephemeris
[Xt, oe, tn] = J2perturbedephemeris(X0, t);
%oe = [a e i raan w M0]

%ecf and subsatellite points
for i=1:length(tn)
    Txyz = ecf2ecisimple(tn(i),ag);
    xt(:,i) = transpose(Txyz)*transpose(Xt(i,1:3));
    latlonalt(i,:) = ecf2spherical(transpose(xt(:,i)));
    latlonalt(i,1) = rad2deg(latlonalt(i,1));
    latlonalt(i,2) = rad2deg(latlonalt(i,2));
end

%Jacobi constant
for i=1:length(tn)
    R = Xt(i,1:3);
    V = Xt(i,4:6);
    h(i,:) = cross(R,V);%specific angular momentum
    R = norm(R);
    V = norm(V);
    Uprime = (-mu/R)*((re/R)^2)*J2*(((3*sind(latlonalt(i,1))^2) - 1)/2);
    K(i) = V^2/2 - (dot([0 0 we],h(i,:)) + (mu/R) + Uprime);
end

%periodic and secular variations of orbital elements
navg = sqrt(mu/oe(1,1)^3);%avg mean motion
aavg = oe(1,1);%avg semimajor axis
for i=1:length(oe)
    n = sqrt(mu/oe(i,1)^3);%mean motion
    %1st order estimate of secular d(raan)/dt
    dANdt(i) = (-3/2)*J2*(n/((1-oe(i,2)^2)^2))*((re/oe(i,1))^2)*cosd(oe(i,3));
    dAOPdt(i)= (3/4)*J2*(n/((1-oe(i,2)^2)^2))*((re/oe(i,1))^2)*(5*cosd(oe(i,3)^2 - 1));
    dMdt(i)  = navg + (3/4)*J2*(n/((1-oe(i,2)^2)^2))*((re/oe(i,1))^2)*(5*cosd(oe(i,3)^2 - 1));
end
%compare node rate
for i=2:length(oe)
    %taylor series estimate from numerical integration of d(raan)/dt
    estdANdt(i) = (oe(i,4)-oe(i-1,4))/(tn(i)-tn(i-1));
    estdAOPdt(i) = (oe(i,5)-oe(i-1,5))/(tn(i)-tn(i-1));
    estdMdt(i) = (oe(i,6)-oe(i-1,6))/(tn(i)-tn(i-1));
    estdadt(i) = (oe(i,1)-oe(i-1,1))/(tn(i)-tn(i-1));
    dadt(i) = ((aavg + 3*navg*aavg*J2*((re/oe(i,1))^2)*(sin(oe(i,3))^2)...
        *((cosd(2*oe(i,5) + 2*oe(i,6)))/(2*dAOPdt(i) + 2*dMdt(i))))...
        -(aavg + 3*navg*aavg*J2*((re/oe(i-1,1))^2)*(sin(oe(i-1,3))^2)...
        *((cosd(2*oe(i-1,5) + 2*oe(i-1,6)))/(2*dAOPdt(i-1) + 2*dMdt(i-1)))));
    dadt(i) = dadt(i)/(tn(i)-tn(i-1));
end



