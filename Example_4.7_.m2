restart
m = 3; n = 3; -- m x n matrices
KK = ZZ/101;
R = KK[x_(1,1)..x_(m,n)]**KK[y_(1,1)..y_(m,n)]**KK[a_1..a_m,b_1..b_n];

xx = toList(x_(1,1)..x_(m,n));
matX = matrix for i in 1..m list for j in 1..n list x_(i,j);
matY = matrix for i in 1..m list for j in 1..n list y_(i,j);

compoundMatrix = (M,s) -> (
    nrow := numRows(M);
    ncol := numColumns(M);
    if s> min(nrow,ncol) then return "error: s larger than min(nrow,ncol)";
    srow := sort subsets(toList(0..nrow-1),s);
    scol := sort subsets(toList(0..ncol-1),s);
    return matrix for i in srow list for j in scol list det((M_i)^j)
    )

s = 2;
cM2 = compoundMatrix(matX,s);
J = trim ideal flatten entries(cM2-transpose(cM2));

decJ = decompose J;
-- for (m,n)=(3,3) and s=2, decJ has two irreducible components.
-- one irreducible component is the subspace of symmetric matrices
-- we check that the second irreducible component is the dual variety of PP^2xPP^2 relative to the subspace of rank one symmetric matrices

r = 2;
IX = minors(r, matX);

jacX = diff(matrix{xx},transpose gens IX);
IXsing = minors(r-1, matX);

-- ideal of the conormal variety of X (see Ottaviani, Spaelenhauer, Sturmfels)
IconX = IX+minors(m+2-r,matY)+(ideal flatten entries(transpose(matX)*matY))+(ideal flatten entries(matX*transpose(matY)));
multidegree IconX

suby = flatten apply(m, i-> apply(n, j-> y_(i+1,j+1)=>x_(i+1,j+1)));

-- ideal of the dual variety of X
dualX = sub(eliminate(xx, IconX), suby);

IY = trim ideal(flatten entries(matX-transpose(matX))); IZ = IX + IY; -- Y subspace of symmetric matrices
IY = trim ideal(flatten entries(matX+transpose(matX))); IZ = IX + IY; -- Y subspace of antisymmetric matrices

-- ideal of the conormal variety of X relative to Z
time IconXZ = saturate(IZ+minors(m+2-r,matY)+(ideal flatten entries(transpose(matX)*matY))+(ideal flatten entries(matX*transpose(matY))), IXsing);
multidegree IconXZ

-- ideal of the dual variety of X relative to Z
dualXZ = sub(eliminate(xx, IconXZ), suby);

(codim dualX, degree dualX)
(codim dualXZ, degree dualXZ)

dualXZ == decJ#1
