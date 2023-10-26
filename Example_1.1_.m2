restart
n = 3; -- dimension of affine space (hence, we work in PP^(n-1))
d = 2; -- degree of the Veronese embedding
N = binomial(n-1+d,d); -- the d-th Veronese embedding of PP^(n-1) lives in PP(Sym^d(C^n)) = PP^(N-1) 
KK = ZZ/101;

-- we order the variables of PP(Sym^d(C^n)) nicely
S = KK[b_1..b_n];
symb = first entries symmetricPower(d,matrix{{b_1..b_n}});
E = apply(symb, i-> toSequence((exponents(i))#0));
--

R = KK[toList(a_1..a_n)|apply(E, e-> x_e)|apply(E, e-> y_e)|{t}];
syma = symmetricPower(d,matrix{{a_1..a_n}});
matx = matrix{apply(E, e-> x_e)};
maty = matrix{apply(E, e-> y_e)};

-- ideal of the d-th Veronese embedding of PP^(n-1), called X
IX = eliminate(toList(a_1..a_n),ideal(first entries(syma-matx)));
jacX = diff(matx,transpose gens IX);
-- ideal of the conormal variety of X
IconX = saturate(IX + minors((codim IX)+1,maty||jacX), ideal(matx));
-- ideal of the dual variety of X
IdualX = eliminate(first entries matx, IconX);

-- (for n = 3) consider a smooth conic C in PP^2
IC = ideal(a_1^2-a_2*a_3)
-- ideal of the d-th Veronese embedding of C, called Z
IZ = eliminate(toList(a_1..a_n),IC+ideal(first entries(syma-matx)));
jacZ = diff(matx,transpose gens IZ);
-- ideal of the conormal variety of Z
IconZ = saturate(IZ + minors((codim IZ)+1,maty||jacZ), ideal(matx));
-- ideal of the dual variety of Z
IdualZ = eliminate(first entries matx, IconZ);

-- ideal of the conormal variety of X relative to Z
IconXZ = saturate(IZ + minors((codim IX)+1,maty||jacX), ideal(matx));
-- ideal of the dual variety of X relative to Z, called X_Z^\vee
IdualXZ = eliminate(first entries matx, IconXZ);


(codim IdualX, degree IdualX) -- (1,3)
(codim IdualZ, degree IdualZ) -- (1,6)
(codim IdualXZ, degree IdualXZ) -- (2,6)

dec = decompose(IdualX+IdualZ);
#dec -- the intersection X^\vee \cap Z^\vee has two irreducible components
dec#0 == IdualXZ -- one of the irreducible components is X_Z^\vee
(codim(dec#1), degree(dec#1)) -- also the second irreducible component has (codimension, degree) = (2,6)

-- geometrically, the intersection X^\vee \cap Z^\vee parametrizes singular plane conics tangent to C.
-- It has 2 irreducible components, coming from the two scenarios:
-- 1) one of the lines is tangent to C, and the other line is arbitrary
-- 2) the two lines meet at a point of C. This case corresponds to the relative dual variety X_Z^\vee
