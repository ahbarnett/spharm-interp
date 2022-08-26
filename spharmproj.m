function cnm = spharmproj(u, z, w, P)
% SPHARMPROJ  projects function on tensor grid on S2 into degree<=P sph harms
%
% cnm = spharmproj(u, z, w, P)
% Inputs:
%  u = complex M*N array of values on grid, with z being fast (1st index), phi
%        slow (2nd index). The phi grid values are assumed to be 2pi*(0:N-1)/N,
%        ie, zero-indexed.
%  z,w = length-M lists of quadrature nodes and weights in z=cos(theta)
%  P = maximum degree of harmonic projection.
% Outputs:
%  cnm = complex coefficients returned in rectangular (P+1)-by-(2P+1) array;
%        only around half the entries are used, due to the triangular shape.
%        Eg c(0,0) is cnm(1,P+1), while c(1,-1) is at cnm(2,P), etc.
%
% See also: FLATTENCNM

[nz nph] = size(u);
if numel(z)~=nz, error('number of z nodes must match # rows of u!'); end
if numel(w)~=nz, error('number of z weights must match # rows of u!'); end
if isreal(u), u = complex(u); end
u = u.';               % fix to Leslie's phi-fast z-slow ordering.
P1 = P+1; P2 = 2*P+1;  % size of local coeffs output array

mex_id_ = 'projloc3d(i int, i int, i int, i int, i double[x], i double[x], i dcomplex[xx], o dcomplex[xx])';
[cnm] = gateway(mex_id_, P, P, nz, nph, z, w, u, nz, nz, nph, nz, P1, P2);

% tidy up cnm: just zero the bad values for now incase are unassigned
for n=0:P, for m=[-P:-n-1, n+1:P], cnm(n+1,P+1-m) = 0; end, end


% =============================================================================
