function [R, z, Ez, N] = computeR(L,xloc,dz,a)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = floor(L/dz)+1;
ntotalsegments = sum(N);
sumN = cumsum([0 N]);
nelements = length(L);

xv = zeros(ntotalsegments,1)*nan;
z = zeros(ntotalsegments,1)*nan;
Ez = zeros(ntotalsegments,1);
Ez((sumN(2)-1)/2+1) = 1;

for i = 2:nelements+1
    z( sumN(i-1) + 1 : sumN(i)) = -L(i-1)/2:dz:L(i-1)/2;
    xv(sumN(i-1) + 1 : sumN(i)) = xloc(i-1)*ones(N(i-1),1);
end

dx = xv-xv';
dx = dx + (dx==0)*a;

R = sqrt((z-z').^2 + dx.^2);
end

