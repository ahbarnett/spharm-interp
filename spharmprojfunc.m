function [cnm err] = spharmprojfunc(f,P,task)
% SPHARMPROJFUNC  Get spherical harmonic expansion coeffs of a function.
%
% cnm = spharmprojfunc(f,P) returns cnm coefficients of function f(z,phi) up to
%  maximum degree P. f is a handle to a function of z and phi which must be
%  able to handle vector or array inputs (z and phi arrays of same size).
%
% cnm = spharmprojfunc(f,P,'show') also plots f on the sphere.
%
% [cnm err] = spharmprojfunc(...) also returns the root mean square error of the
%  truncated expansion on the projection grid; this indicates how accurate the
%  expansion is for the given function (since ppw>>2, aliasing is small).
%
% Without arguments, a self-test is done.
%
% See also: LSQSOLVESPHARM for definition of expansion
%           SPHARMPROJ

% Barnett 8/21/15
% made nz!=nph 9/1/15
if nargin==0, test_spharmprojfunc; return; end
PP = (P+1)^2;
ppw = 5.0;          % eval on overresolved grid to get coeffs
N = round(ppw*P);   % # phi nodes
M = round(ppw*P/2); % # z nodes, scaled to # phi nodes by Nyquist
[zg wg] = gauss(M);
phig =  2*pi*(0:N-1)/N;
[pphi zz] = meshgrid(phig,zg);
ug = f(zz,pphi);
if nargin>2 & strcmp(task,'show'), showsphgrid(ug,zg); end
cnm = flattencnm(spharmproj(reshape(ug, [M N]), zg, wg, P));
if nargout>1
  u = spharmgrideval(cnm,zg,phig);      % check coeffs go back to func on grid
  err = norm(u(:)-ug(:)) / sqrt(M*N);
end

%%%%
function test_spharmprojfunc
P = 30;      % >= n
n = 27; m = 5;                                    % a single spherical harmonic
f = @(z,phi) legendrenm(n,m,z) .* exp(1i*m*phi);
PP=(P+1)^2; ce=zeros(P+1,2*P+1); ce(n+1,P+1+m) = 1; ce = flattencnm(ce); % true
ce = -ce/sqrt(2);                      % empirical prefactor!
[cnm err] = spharmprojfunc(f,P,'show');
fprintf('cnm err vs true coeffs = %.3g,  claimed err = %.3g\n',norm(cnm-ce),err)
