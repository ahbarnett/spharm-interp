function showspharmexp(cnm, o)
% SHOWSPHARMEXP  show a spherical harmonic expansion as color on sphere
%
% showspharmexp(cnm) treats cnm as an unrolled list of (P+1)^2 spherical
%  harmonic coefficients, and evaluates and plots on a grid suitable to capture
%  the finest features. Real part only is shown.
%
% showspharmexp(cnm,opts) controls options:
%  opts.ppw = points per shortest wavelength (controls resolution)
%  opts.sc = color scale, passed to showsphgrid.

% Barnett 8/27/22
if nargin==0, test_showspharmexp; return; end
if nargin<2, o = []; end
if ~isfield(o,'ppw'), o.ppw = 5.0; end   % points per shortest wavelength

P = sqrt(length(cnm))-1;
M = round(o.ppw*P);   % # phi nodes
N = round(o.ppw*P/2);   % # z nodes
zg = fejer(N);
phig = 2*pi*(0:M-1)/M;
ug = spharmgrideval(cnm,zg,phig);
if ~isfield(o,'sc')          % don't or do pass option through
  showsphgrid(ug,zg);
else
  showsphgrid(ug,zg,o.sc);
end

%%%%%%%%%%
function test_showspharmexp
P=20;
cnm = randn((P+1)^2,1) + 1i*randn((P+1)^2,1);
showspharmexp(cnm)
