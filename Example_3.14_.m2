restart
n = 4; -- the associated projective space has dimension n-1
KK = QQ;
R = KK[x_1..x_n]**KK[y_1..y_n];

-- ideal of the discriminant of a binary cubic form
IX = ideal(x_2^2*x_3^2-4*x_1*x_3^3-4*x_2^3*x_4+18*x_1*x_2*x_3*x_4-27*x_1^2*x_4^2)

jacX = diff(matrix{{x_1..x_n}},transpose gens IX);
IXsing = trim(IX+minors(codim(IX), jacX));

-- (1) intersection of X with a generic complete intersection variety Y
cY = 2 -- codimension of Y
DY = {3,1} -- degrees of the generators of the ideal of Y
LY =  apply(cY, i-> (symmetricPower(DY#i,matrix{{x_1..x_n}})*random(KK^(binomial(n+DY#i-1,DY#i)),KK^1))_(0,0))
IY = trim ideal LY;
IZ = IY + IX;

-- (2) line Z contained in X
IZ = ideal(x_1,x_2);

-- (3) point Z contained in X
IZ = ideal(x_1,x_2,random(KK)*x_3+random(KK)*x_4);

-- ideal of the conormal variety of X (in PP^{n-1}xPP^{n-1})
IconX = saturate(IX + minors((codim IX)+1,matrix{{y_1..y_n}}||jacX), IXsing);
multidegree IconX

-- the following is the projective conormal variety of X relative to Z
IconXZ = saturate(IZ + minors((codim IX)+1,matrix{{y_1..y_n}}||jacX), IXsing);
multidegree IconXZ

suby = apply(n, i-> y_(i+1)=>x_(i+1));

dualX = sub(eliminate(toList(x_1..x_n), IconX), suby);
dualXZ = sub(eliminate(toList(x_1..x_n), IconXZ), suby);

(codim dualX, degree dualX)
(codim dualXZ, degree dualXZ)
