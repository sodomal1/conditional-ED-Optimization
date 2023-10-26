restart
n = 3; -- dimension of the affine ambient space
KK = QQ;
R = KK[x_1..x_n]**KK[u_1..u_n]
xx = matrix{{x_1..x_n}}; uu = matrix{{u_1..u_n}};

-- ideal of Cayley's cubic surface X
IX = ideal(4*(x_1^2+x_2^2+x_3^2)+16*x_1*x_2*x_3-1)

-- ideal of the singular locus of X
jacX = diff(xx,transpose gens IX);
IXsing = trim(IX+minors(codim(IX), jacX));

-- (1) Z is a circle in X
IZ = IX + ideal(x_1)

-- (2) Z is a line in X joining two of its singular points
IZ = ideal(x_1-x_2,x_3+1/2)

-- ideal of the ED correspondence variety of X
ux = uu-xx
IEDcorrX = saturate(IX + minors((codim IX)+1,ux||jacX), IXsing);

-- ideal of the ED correspondence variety of X relative to Z
IEDcorrXZ = saturate(IZ + minors((codim IX)+1,ux||jacX), IXsing);

-- data Locus
DLXZ = eliminate(toList(x_1..x_n), IEDcorrXZ)

-- the variety of the ideal IEDcorrX+IZ might have several irreducible components
decompose eliminate(toList(x_1..x_n), IEDcorrX+IZ)
