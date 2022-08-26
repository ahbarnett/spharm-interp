function [x,w] = fejer(N)
% FEJER  nodes x (1/2-offset Chebyshev points) and weights w
%        for N-node Fejer's "first-rule" quadrature on [-1,1]
%
% [x,w] = fejer(N)
%
% nodes x are in increasing order, and are "half-way between" usual Chebyshev
% nodes.
%
% Called without arguments does self-test.
%
% Modified from clencurt.m by a 1/2-node phase shift and not halving the endpts.
% See Waldvogel 2006 BIT paper.
% Alex Barnett 8/20/15
if nargin==0, test_fejer; return; end
theta = pi*(N-.5:-1:.5)'/N; x = cos(theta);
k = 1:floor(N/2);                             % frequency indices
W = kron(exp(1i*pi*k/N)./(1/4-k.^2), [0 1]);
if mod(N,2)==1          % pi-freq term is zero
  w = ifft([4 W 0 conj(W(end:-1:1))]);  % 4 is the zero-freq term  
else, W(end) = real(W(end));   % pi-freq term should be made real
  w = ifft([4 W conj(W(end-1:-1:1))]);
end
w = w(1:N);

%%%%
function test_fejer
f = @(x) sin(1+x); fint = cos(0)-cos(2);
[x w] = fejer(20); w*f(x)-fint   % even
[x w] = fejer(21); w*f(x)-fint   % odd
