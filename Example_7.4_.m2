restart
m = 2; n = 3; -- m x n matrices
KK = ZZ/101;
R = KK[x_(1,1)..x_(m,n)]**KK[y_(1,1)..y_(m,n)]**KK[a_1..a_m,b_1..b_n];

xx = toList(x_(1,1)..x_(m,n));
matX = matrix for i in 1..m list for j in 1..n list x_(i,j);
matY = matrix for i in 1..m list for j in 1..n list y_(i,j);

-- ideal of the determinantal variety X of m x n matrices of rank at most one
IX = minors(2, matX);

-- ideal of the singular locus of X
IXsing = minors(1, matX);

-- consider a subvariety of the type L1xL2, L1, L2 subspaces
-- in the example, we consider (c1,c2)=(0,1)
Isegre = ideal flatten apply(m, i-> apply(n, j-> x_(i+1,j+1)-a_(i+1)*b_(j+1)));
c1 = 0;
c2 = 1;
IH1 = ideal flatten entries (random(KK^c1,KK^m)*transpose(matrix{{a_1..a_m}}));
IH2 = ideal flatten entries (random(KK^c2,KK^n)*transpose(matrix{{b_1..b_n}}));
IZ = eliminate(toList(a_1..a_m)|toList(b_1..b_n),IH1+IH2+Isegre);

-- ideal of the conormal variety of X relative to Z
time IconXZ = saturate(IZ+minors(m,matY)+(ideal flatten entries(transpose(matX)*matY))+(ideal flatten entries(matX*transpose(matY))), IXsing);
multidegree IconXZ

xy = matrix{apply(flatten entries(matX-matY), i-> random(KK)*i)}; -- (run this line to use a generic inner product)
xy = matrix{flatten entries(matX-matY)}; -- (run this line to use the Frobenius inner product)

-- ideal of the conditional ED correspondence variety of X given Z
jacX = diff(matrix{{x_(1,1)..x_(m,n)}},transpose gens IX);
IEDcorrXZ = saturate(IZ + minors((codim IX)+1,xy||jacX), IXsing);

-- ideal of the ED data locus of X given Z
DLXZ = eliminate(xx, IEDcorrXZ);
