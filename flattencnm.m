function a = flattencnm(c)
% FLATTENCNM  convert a (P+1)-by-(2P+1) storage of C_nm into linear column vec
%
% No arguments, runs a self-test

% Barnett 8/1/15
if nargin==0, test_flattencnm; return; end

P = size(c,1)-1;
if size(c,2)~=2*P+1, error('c is not correct size (P+1)-by-(2P+1)'); end
a = zeros((P+1)^2,1);  % output
count = 1;
for n=0:P, for m=-n:n, a(count) = c(n+1,m+P+1); count=count+1; end, end

%%%%%
function test_flattencnm
fprintf('testing flattencnm, should give zero:\n')
P = 10; c = nan(P+1,2*P+1);
count=1; for n=0:P, for m=-n:n, c(n+1,m+P+1) = count; count = count+1; end, end
a = flattencnm(c);
norm(a - (1:(P+1)^2)')
