function newpsi = nystrom_extension(psi, X, Y, sig, Lambda)

d = pdist2(X', Y')';
A = exp(-d.^2/(2*sig^2));
M = bsxfun(@rdivide, A, sum(A,2));
Mpsi = M*psi;
newpsi = bsxfun(@times, Mpsi', 1./Lambda);

