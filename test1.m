c = physconst('lightspeed');
freq = 300e6;
lambda = c/freq;

L1 = lambda/2; % length of the dipole
dz = lambda/40; % length of a discreet segment

z = [-L1/2:dz:L1/2]';
N = length(z);

a = 0.0001; % radius of the wires

R = sqrt((z-z').^2 + a^2);

k = 2*pi/lambda; % wavenumber

G1 = (- 1 - 1i*k*R + k^2*R.^2) / (R.^3);
G2 = (  3 + 3i*k*R - k^2*R.^2) / (R.^5);

A = (G1 + (z-z').^2.*G2).*exp(-1i*k*R);

% According to the slides
%   Enforcing this boundary condition for the E component will yield
%   zero at every observation point, except at the feed point to the dipole
%   where the field is assumed to be unity (scaled).
Ez = [zeros((N-1)/2,1); 1; zeros((N-1)/2,1)];

Jz = A\Ez


figure(1)
stem(z,abs(Jz))
xlabel('z')
ylabel('|J|')