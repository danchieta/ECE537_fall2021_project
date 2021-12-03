function [R, z, Ez, N] = computeR(L,xloc,dz,a)
%computeR Summary of this function goes here
%   L - length of the elements
%   xloc - location of the elements in the x axis
%   dz - size of each segment in the z axis
%   a - radius of the elements

if length(L)~=length(xloc)
    error('L and xloc must have the same length.')
end

N = floor(L/dz)+1;          % Number of segments at each element

ntotalsegments = sum(N);    % Total number of segments
sumN = cumsum([0 N]);       % will use this to calculate z and x coordinates for each segment
nelements = length(L);      % total number of elements

% Initializing vectors of coordinates
xv = zeros(ntotalsegments,1)*nan;
z = zeros(ntotalsegments,1)*nan;

% According to the slides
%   Enforcing this boundary condition for the E component will yield
%   zero at every observation point, except at the feed point to the dipole
%   where the field is assumed to be unity (scaled).
Ez = zeros(ntotalsegments,1);
Ez(floor(sumN(2)/2)+1) = 1;

% for each element calculate the coordinates of the segments
for i = 2:nelements+1
    %z( sumN(i-1) + 1 : sumN(i)) = -L(i-1)/2 : dz: L(i-1)/2;
    z( sumN(i-1) + 1 : sumN(i)) = linspace(-L(i-1)/2,L(i-1)/2,N(i-1));
    xv(sumN(i-1) + 1 : sumN(i)) = xloc(i-1)*ones(N(i-1),1);
end

% x-axis distance between segments
dx = xv-xv';
dx = dx + (dx==0)*a; % add a whenever distance is zero

R = sqrt((z-z').^2 + dx.^2);
end

