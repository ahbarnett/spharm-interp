function c = fastapplyAT(y,L,P,zg,wg)
% c = fastapplyAT(y,L,P,zg,wg) applies Hermitian transpose of A to a data
%  vector y, returning truncated Fourier coeffs in -P:P. The sparse
%  interpolation matrix L (returned by fillsparseL) is needed, along with the
%  intermediate projection grid z-nodes zg, and weights wg.
%
% Barnett 8/21/15
n = size(L,2);
Fg = L'*y;
M = numel(zg); N = round(n/M);  % projection grid size
% note ahb debugged 2*M to 2*N in the following, so correct adjoint (6/24/16):
Fg = reshape(Fg,[M N]) .* repmat((2*N)./wg(:), [1 N]);  % undo proj weights
% ...in order to transpose the sph harm grid eval op (prefac here is guessed).
c = flattencnm(spharmproj(Fg, zg, wg, P));  % project & order the coeffs
