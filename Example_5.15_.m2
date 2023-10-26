restart
n = 3; -- dimension of the affine ambient space
KK = QQ;
R = KK[x_1..x_n]**KK[u_1..u_n];
xx = matrix{{x_1..x_n}}; uu = matrix{{u_1..u_n}};

-- ideal of the smooth quadric X
IX = ideal(x_1-x_2*x_3);

-- computing the singular locus of X
jacX = diff(xx,transpose gens IX);
IXsing = trim(IX+minors(codim(IX), jacX));

-- ideal of the line Z contained in X
IZ = ideal(x_1,x_2);

-- ideal of the conditional ED correspondence variety of X given Z
IEDcorrXZ = saturate(IZ + minors((codim IX)+1,ux||jacX), IXsing);

-- ideal of the conditional data locus of X given Z
DLXZ = eliminate(toList(x_1..x_n), IEDcorrXZ)
