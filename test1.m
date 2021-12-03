c = physconst('lightspeed');
freq = 300e6;
lambda = c/freq;

L = [2*lambda, lambda, lambda];      % length of parasite 2
x = [0, lambda/4, lambda/2];      % location of wires on the x axis
dz = lambda/45;    % length of a discrete segment
a = 0.005;          % radius of the wires

t0 = atan(dz/2/a);

[R, z, Ez, N] = computeR(L,x,dz,a);

k = 2*pi/lambda; % wavenumber

G1 = (- 1 - 1i*k*R + k^2*R.^2) ./ (R.^3);
G2 = (  3 + 3i*k*R - k^2*R.^2) ./ (R.^5);

A = (G1 + (z-z').^2.*G2).*exp(-1i*k*R);
% A = A - diag(diag(A)) + A(1,2)*1.5;
% A = A - diag(real(diag(A)));

[M,~] = size(A);
for m=1:M
    A(m,m) = (2/(a^2*dz))*(sin(t0*((a*k)^2)-1) + (1j*a*k*t0)/2 - ...
        (3/4)*(1j*a*k*sin(2*t0)) + (sin(t0)^3));
end

Jz = A\Ez

dummyN = cumsum([0 N]);
nfigures = length(L);

figure(1)
clf
% Jz = abs(real(Jz));
for i=2:nfigures+1
    subplot(nfigures,1,i-1)
    plot(z(dummyN(i-1) + 1 : dummyN(i)), db(abs(real(Jz(dummyN(i-1) + 1 : dummyN(i))))))
    hold on
    plot(z(dummyN(i-1) + 1 : dummyN(i)), db(abs(imag(Jz(dummyN(i-1) + 1 : dummyN(i))))))
	xlim([min(z) max(z)]*1.05)
%     ylim(db([min(Jz) max(Jz)]))
    grid on
    xlabel('z')
    ylabel('J (dB)')
end
legend('Real','Imag')

figure(2)
clf
for i=2:nfigures+1
    subplot(nfigures,1,i-1)
    plot(z(dummyN(i-1) + 1 : dummyN(i)), db(abs(Jz(dummyN(i-1) + 1 : dummyN(i)))))
	xlim([min(z) max(z)]*1.05)
    ylim(db([min(abs(Jz)) max(abs(Jz))]))
    xlabel('z')
    ylabel('|J| (dB)')
end

figure(3)
imagesc(db(abs(A)))
axis equal
colorbar
title('A')

for i=1:3
    figure(i)
    exportgraphics(gcf, sprintf('fig_%d.png', i), 'resolution', 400)
end