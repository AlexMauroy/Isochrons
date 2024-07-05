function dX=van_der_pol(t,X,param)

mu=param.mu;

if size( X, 2) >1
    X=X(:);
end

nPoints = length( X) / 2;
ix = (0 * nPoints + 1) : (1 * nPoints);
iy = (1 * nPoints + 1) : (2 * nPoints);

dX=[X(iy);mu*(1-X(ix).^2).*X(iy)-X(ix)];
