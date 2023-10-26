restart
m = 3; n = 3; -- m x n matrices
KK = ZZ/101;
R = KK[x_(1,1)..x_(m,n)]**KK[y_(1,1)..y_(m,n)]**KK[a_1..a_m,b_1..b_n];

matX = matrix for i in 1..m list for j in 1..n list x_(i,j);
matY = matrix for i in 1..m list for j in 1..n list y_(i,j);

xx = toList(x_(1,1)..x_(m,n));

-- ideal of the determinantal variety of matrices of rank at most r-1
r = 3; IX = minors(r, matX);

-- ideal of the singular locus of X
jacX = diff(matrix{{x_(1,1)..x_(m,n)}},transpose gens IX);
IXsing = minors(r-1, matX);

-- ideal of the conormal variety of X (see Ottaviani, Spaelenhauer, Sturmfels)
IconX = IX+minors(m+2-r,matY)+(ideal flatten entries(transpose(matX)*matY))+(ideal flatten entries(matX*transpose(matY)));
multidegree IconX

suby = flatten apply(m, i-> apply(n, j-> y_(i+1,j+1)=>x_(i+1,j+1)));

-- ideal of the dual variety of X.
-- It coincides with the determinantal variety of corank r-1 matrices
dualX = sub(eliminate(xx, IconX), suby);

-- (1) intersection of X with a generic complete intersection variety Y
cY = 2; -- codimension of Y
DY = {1,1}; -- degrees of the generators of the ideal of Y
LY =  apply(cY, i-> (symmetricPower(DY#i,matrix{{x_(1,1)..x_(m,n)}})*random(KK^(binomial(m*n+DY#i-1,DY#i)),KK^1))_(0,0));
IY = trim ideal LY;
IZ = IY + IX;

-- (2) subvariety of matrices with the first row equal to zero
IZ = ideal(x_(1,1),x_(1,2),x_(1,3));

-- (3) another subvariety Z
IZ = ideal(x_(1,1),x_(1,2),x_(2,1)*x_(3,2)-x_(2,2)*x_(3,1));

-- ideal of the conormal variety of X (in PP^{n-1}xPP^{n-1})
time IconXZ = saturate(IZ+minors(m+2-r,matY)+(ideal flatten entries(transpose(matX)*matY))+(ideal flatten entries(matX*transpose(matY))), IXsing);
multidegree IconXZ

-- ideal of the dual variety of X relative to Z
dualXZ = sub(eliminate(xx, IconXZ), suby);

(codim dualX, degree dualX)
(codim dualXZ, degree dualXZ)
