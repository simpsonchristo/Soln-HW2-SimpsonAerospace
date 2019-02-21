clear;
clc;
%example 2.3.1 (p. 61)
%J2 propagation (perturbations)
t = 0:60:11000;%seconds
X0 = [5492000.34 3984001.40 2955.81 -3931.046491 5498.676921 3665.980697];
[Xt, oe, tn] = J2perturbedephemeris(X0, t);
%oe = [a e i raan w M0]
h=figure('Name','Orbital Elements');
subplot(3,2,1)
plot(tn,oe(:,1))
xlabel('Time (sec)')
ylabel('Semimajor axis (m)')

subplot(3,2,2)
plot(tn,oe(:,2))
xlabel('Time (sec)')
ylabel('Eccentricity')
saveas(h,'ecc.png')

subplot(3,2,3)
plot(tn,oe(:,3))
xlabel('Time (sec)')
ylabel('Inclination (deg)')

subplot(3,2,4)
plot(tn,oe(:,4))
xlabel('Time (sec)')
ylabel('Ascending node (deg)')

subplot(3,2,5)
plot(tn,oe(:,5))
xlabel('Time (sec)')
ylabel('Arg. of Perigee (deg)')

saveas(h,'ex231_oe.png')