c = physconst('lightspeed');
freq = 300e6;
lambda = c/freq;

L1 = lambda/2; % length of the dipole
dz = lambda/20; % length of a discreet segment

z = [-L1/2:dz:L1/2]';

a = 0.01; % radius of the wires

R = sqrt((z-z').^2 + a^2);

k = 2*pi/lambda; % wavenumber

G1 = (- 1 - 1i*k*R + k^2*R.^2) / (R.^3);
G2 = (  3 + 3i*k*R - k^2*R.^2) / (R.^5);

A = (G1 + (z-z').^2.*G2).*exp(-1i*k*R);

