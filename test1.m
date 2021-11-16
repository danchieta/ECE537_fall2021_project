c = physconst('lightspeed');
freq = 300e6;
lambda = c/freq;

L = [lambda/2,...   % length of the dipole
    lambda/4,...    % length of parasite 1
    lambda/5,...    % length of parasite 2
    lambda/5];      % length of parasite 3
x = [0,...
    lambda/8,...
    lambda/4,...
    lambda/2];      % location of wires on the x axis
dz = lambda/40;     % length of a discrete segment
a = 0.005;          % radius of the wires

[R, z, Ez, N] = computeR(L,x,dz,a);

k = 2*pi/lambda; % wavenumber

G1 = (- 1 - 1i*k*R + k^2*R.^2) / (R.^3);
G2 = (  3 + 3i*k*R - k^2*R.^2) / (R.^5);

A = (G1 + (z-z').^2.*G2).*exp(-1i*k*R);

Jz = A\Ez

dummyN = cumsum([0 N]);
nfigures = length(L);

figure(1)
clf
for i=2:nfigures+1
    subplot(nfigures,1,i-1)
    stem(z(dummyN(i-1) + 1 : dummyN(i)), abs(Jz(dummyN(i-1) + 1 : dummyN(i))))
	xlim([min(z) max(z)])
    ylim([min(abs(Jz)) max(abs(Jz))])
    xlabel('z')
    ylabel('|J|')
end