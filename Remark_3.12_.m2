restart
m = 2; n = 3; -- m x n matrices
KK = ZZ/101;
R = KK[x_(1,1)..x_(m,n)]**KK[y_(1,1)..y_(m,n)]**KK[a_1..a_m,b_1..b_n];

xx = toList(x_(1,1)..x_(m,n));
yy = toList(y_(1,1)..y_(m,n));

matX = matrix for i in 1..m list for j in 1..n list x_(i,j);
matY = matrix for i in 1..m list for j in 1..n list y_(i,j);

-- ideal of the determinantal variety X of m x n matrices of rank at most r-1
r = 2; IX = minors(r, matX);

jacX = diff(matrix{xx},transpose gens IX);
IXsing = minors(r-1, matX);

-- ideal of the conormal variety of X (see Ottaviani, Spaelenhauer, Sturmfels)
IconX = IX+minors(m+2-r,matY)+(ideal flatten entries(transpose(matX)*matY))+(ideal flatten entries(matX*transpose(matY)));
multidegree IconX

suby = flatten apply(m, i-> apply(n, j-> y_(i+1,j+1)=>x_(i+1,j+1)));
dualX = sub(eliminate(xx, IconX), suby);

-- consider a subvariety of the type L1xL2, L1, L2 subspaces
-- in the remark, we consider first (c1,c2)=(1,0), then (c1,c2)=(0,2)
Isegre = ideal flatten apply(m, i-> apply(n, j-> x_(i+1,j+1)-a_(i+1)*b_(j+1)));
c1 = 1
c2 = 0
IH1 = ideal flatten entries (random(KK^c1,KK^m)*transpose(matrix{{a_1..a_m}}));
IH2 = ideal flatten entries (random(KK^c2,KK^n)*transpose(matrix{{b_1..b_n}}));
IZ = eliminate(toList(a_1..a_m)|toList(b_1..b_n),IH1+IH2+Isegre);

-- ideal of the projective conormal variety of X relative to Z
time IconXZ = saturate(IZ+minors(m+2-r,matY)+(ideal flatten entries(transpose(matX)*matY))+(ideal flatten entries(matX*transpose(matY))), IXsing);
multidegree IconXZ

-- ideal of the dual variety of X relative to Z
dualXZ = sub(eliminate(xx, IconXZ), suby);

(codim dualX, degree dualX)
(codim dualXZ, degree dualXZ)
