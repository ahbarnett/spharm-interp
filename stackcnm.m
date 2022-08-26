function c = stackcnm(a)
% STACKCNM  convert single vector storage of c_nm into (P+1)-by-(2P+1) array
%
% No arguments, runs a self-test

% Barnett 9/1/15
if nargin==0, test_stackcnm; return; end

P = round(sqrt(numel(a))-1);
if numel(a)~=(P+1)^2, error('a is not of length (P+1)^2 for any integer P'); end
c = zeros(P+1,2*P+1); % output
count = 1;
for n=0:P, for m=-n:n, c(n+1,m+P+1) = a(count); count=count+1; end, end

%%%%%
function test_stackcnm
fprintf('testing stackcnm, should give zero:\n')
P = 10; PP = (P+1)^2; c = randn(PP,1)+1i*randn(PP,1);
norm(c - flattencnm(stackcnm(c)))
