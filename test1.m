c = physconst('lightspeed');
freq = 300e6;
lambda = c/freq;

L1 = lambda/2; % length of the dipole
L2 = lambda/4;
dz = lambda/40; % length of a discrete segment

z1 = [-L1/2:dz:L1/2]';
z2 = [-L2/2:dz:L2/2]';
z = [z1; z2];
N = length(z);
x1 = ones(length(z1),1)*0;
x2 = ones(length(z2),1)*lambda/8;
x = [x1; x2];

a = 0.0001; % radius of the wires

dx = x-x';
dx = dx + (dx==0)*a;

R = sqrt((z-z').^2 + dx.^2);

k = 2*pi/lambda; % wavenumber

G1 = (- 1 - 1i*k*R + k^2*R.^2) / (R.^3);
G2 = (  3 + 3i*k*R - k^2*R.^2) / (R.^5);

A = (G1 + (z-z').^2.*G2).*exp(-1i*k*R);

% According to the slides
%   Enforcing this boundary condition for the E component will yield
%   zero at every observation point, except at the feed point to the dipole
%   where the field is assumed to be unity (scaled).
Ez = [zeros((length(z1)-1)/2,1); 1; zeros(N-(length(z1)+1)/2,1)];

Jz = A\Ez


figure(1)
subplot(211)
stem(z1,abs(Jz(1:length(z1))))
subplot(212)
stem(z2,abs(Jz(length(z1)+1:end)))
xlim([min(z1) max(z1)])
xlabel('z')
ylabel('|J|')