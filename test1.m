c = physconst('lightspeed');
freq = 300e6;
lambda = c/freq;

L = [lambda/2,...   % length of the dipole
    lambda/4,...    % length of parasite 1
    lambda/8];      % length of parasite 2
x = [0,...
    lambda/8,...
    lambda/4];  % location of wires on the x axis
dz = lambda/40;     % length of a discrete segment
a = 0.0001;         % radius of the wires

[R, z, Ez, N] = computeR(L,x,dz,a);

k = 2*pi/lambda; % wavenumber

G1 = (- 1 - 1i*k*R + k^2*R.^2) / (R.^3);
G2 = (  3 + 3i*k*R - k^2*R.^2) / (R.^5);

A = (G1 + (z-z').^2.*G2).*exp(-1i*k*R);

% According to the slides
%   Enforcing this boundary condition for the E component will yield
%   zero at every observation point, except at the feed point to the dipole
%   where the field is assumed to be unity (scaled).

Jz = A\Ez

dummyN = cumsum([0 N]);
nfigures = length(L);


figure(1)
for i=2:nfigures+1
    subplot(nfigures,1,i-1)
    stem(z(dummyN(i-1) + 1 : dummyN(i)), abs(Jz(dummyN(i-1) + 1 : dummyN(i))))
	xlim([min(z) max(z)])
end