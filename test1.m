c = physconst('lightspeed');
freq = 300e6;
lambda = c/freq;

L = [4*lambda];...   % length of parasite 2
x = [0];      % location of wires on the x axis
dz = lambda/100;    % length of a discrete segment
a = 0.005;          % radius of the wires

[R, z, Ez, N] = computeR(L,x,dz,a);

k = 2*pi/lambda; % wavenumber

G1 = (- 1 - 1i*k*R + k^2*R.^2) ./ (R.^3);
G2 = (  3 + 3i*k*R - k^2*R.^2) ./ (R.^5);

A = (G1 + (z-z').^2.*G2).*exp(-1i*k*R);

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
    legend('Real','Imag')
end

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
imagesc(abs(A))
colorbar
title('10log10(A)')