restart
n = 4; -- the associated projective space has dimension n-1
KK = QQ;
R = KK[x_1..x_n]**KK[y_1..y_n];
xx = toList(x_1..x_n); yy = toList(y_1..y_n);

-- ideal of the projective closure of Cayley's cubic surface X
IX = ideal(16*x_1*x_2*x_3+4*x_4*(x_1^2+x_2^2+x_3^2)-x_4^3);

-- ideal of the singular locus of X
jacX = diff(matrix{xx},transpose gens IX);
IXsing = trim(IX+minors(codim(IX), jacX));

-- ideal of the projective closure of a line contained in X and joining two of its singular points
IZ = ideal(x_1-x_2,2*x_3+x_4);

-- ideal of the conormal variety of X (in PP^{n-1}xPP^{n-1})
IconX = saturate(IX + minors((codim IX)+1,matrix{yy}||jacX), IXsing);
multidegree IconX

-- ideal of the conormal variety of X relative to Z (in PP^{n-1}xPP^{n-1})
IconXZ = saturate(IZ + minors((codim IX)+1,matrix{yy}||jacX), IXsing);
multidegree IconXZ

suby = apply(n, i-> y_(i+1)=>x_(i+1));

-- ideal of the dual variety of X (seen in the same space of X)
dualX = sub(eliminate(xx, IconX), suby);

-- ideal of the dual variety of X relative to Z (seen in the same space of X)
dualXZ = sub(eliminate(xx, IconXZ), suby);

-- the following lines are meant for checking if the relative dual variety of X given Z
-- lies in the singular locus of the dual variety of X
jacdualX = diff(matrix{yy},transpose gens dualX);
dualXsing = trim(dualX+minors(codim(dualX), jacdualX));

dualXsing+dualXZ == dualXZ -- true
