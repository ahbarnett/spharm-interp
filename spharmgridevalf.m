function u = spharmgridevalf(cnm,z,phi)
% SPHARMGRIDEVALF  evaluate spherical harmonic expansion on tensor grid on S2.
%
% u = spharmgridevalf(cnm,z,phi) evaluates u a M*N complex array given size-M
%  grid of z, and size-N grid of phi, defining a tensor-product grid on the
%  sphere, and the length-(P+1)^2 list of coefficients c_nm.
% 
%     u =    sum      sum  c_nm Y_nm(theta,phi)   where z = cos(theta)
%           n=0..P  |m|<=n
%
%  The fast (row) index of u is z, the slow (col) is phi, and normalization
%  matches projloc3d.
%
% Notes: MEX interface to Leslie's shevalsphere. Run spharmgrideval for test.
%
% See also: SPHARMGRIDEVAL (all-Matlab version)

local = stackcnm(cnm);
[P1 P2] = size(local);
P = P1-1;
nquad = numel(z);
nquadm = numel(phi);
mex_id_ = 'shevalsphere(i dcomplex[xx], o dcomplex[xx], i int, i int, i int, i int, i double[x])';
[u] = gateway(mex_id_, local, P, P, nquad, nquadm, z, P1, P2, nquadm, nquad, nquad);
u = u.';  % make z fast, phi slow

