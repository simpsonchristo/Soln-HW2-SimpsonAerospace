%p3_script plot script
%% plot ground track
hfig=figure('Renderer','opengl');
worldmap('world')
geoshow('landareas.shp')
plotm(latlonalt(:,1), latlonalt(:,2),'r--')
saveas(hfig,'p3_groundtrack.png');

%% plot Jacobi constant and hz
hfig=figure('Name','Jacobi constant');
subplot(1,2,1)
plot(tn, K, 'k--', 'DisplayName', 'Jacobi')
xlabel('Time (sec)')
ylabel('Jacobi Constant')
subplot(1,2,2)
plot(tn, h(:,3), 'k-', 'DisplayName', 'h_z')
xlabel('Time (sec)')
ylabel('h_z (m^2/sec)')
saveas(hfig,'JacobiAndhz.png')

%% plot rates and estimated rates
hfig=figure('Name','Ascending Node Rate Estimated');
subplot(1,2,1)
plot(tn,estdANdt, 'k--', 'DisplayName', 'estimated')
xlabel('Time (sec)')
ylabel('dAN/dt (deg/sec)')
subplot(1,2,2)
plot(tn, dANdt, 'k-', 'DisplayName', 'predicted')
xlabel('Time (sec)')
ylabel('Eq. 2.3.31 (deg/sec)')
saveas(hfig,'compare_dANdt.png')

hfig=figure('Name','Arg. of Perigee Rate Estimated');
subplot(1,2,1)
plot(tn,estdAOPdt, 'k--', 'DisplayName', 'estimated')
xlabel('Time (sec)')
ylabel('dAOP/dt (deg/sec)')
subplot(1,2,2)
plot(tn, dAOPdt, 'k-', 'DisplayName', 'predicted')
xlabel('Time (sec)')
ylabel('Eq. 2.3.32 (deg/sec)')
saveas(hfig,'compare_dAOPdt.png')

hfig=figure('Name','Mean Anomaly Rate Estimated');
subplot(1,2,1)
plot(tn,estdMdt, 'k--', 'DisplayName', 'estimated')
xlabel('Time (sec)')
ylabel('dM/dt (deg/sec)')
subplot(1,2,2)
plot(tn, dMdt, 'k-', 'DisplayName', 'predicted')
xlabel('Time (sec)')
ylabel('Eq. 2.3.33 (deg/sec)')
saveas(hfig,'compare_dMdt.png')

hfig=figure('Name','Semimajor axis Rate Estimated');
subplot(1,2,1)
plot(tn,estdadt, 'k--', 'DisplayName', 'estimated')
xlabel('Time (sec)')
ylabel('da/dt (m/sec)')
subplot(1,2,2)
plot(tn, dadt, 'k-', 'DisplayName', 'predicted')
xlabel('Time (sec)')
ylabel('Taylor-series expansion of Eq. 2.3.34 (m/sec)')
saveas(hfig,'compare_dadt.png')

%% plot orbital elements over time
hfig=figure('Name','Orbital Elements');
subplot(3,2,1)
plot(tn,oe(:,1))
xlabel('Time (sec)')
ylabel('Semimajor axis (m)')

subplot(3,2,2)
plot(tn,oe(:,2))
xlabel('Time (sec)')
ylabel('Eccentricity')
saveas(hfig,'ecc.png')

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
saveas(hfig,'p3_oe.png')
