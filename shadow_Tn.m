// Allows us to scale any divisor by an integer t
function DivisorScaling(divisor,t)
    output := [];
    for D in divisor do
        Append(~output,<D[1],t*D[2]>);
    end for;
    return output;
end function;

// Checks Heegner hypothesis and ideal condition on a list of integers
function HeegnerCheck(N,D)
    return (IsSquare(ResidueClassRing(4*N)!D)) and (#HeegnerForms(N,D) ne 0);
end function;

// This program will allow us to compute the scaled shadow of T2 at any level N
// Recall that 6*Sh(Tn) = (2g-2)6*F_Tn - deg(F_Tn)6*K_X0(N) - 12*Tn(K_X0(N)) + 2*(n+1)*6*K_X0(N).
function ShadowTn(N,n)
	G := Gamma0(N);
	g := Genus(G);
	
	cusps := Cusps(G);
	elliptic_points := EllipticPoints(G);
	
	// We first compute 6*K = D' - 4*D3 - 3*D4
	// NOTE: To concatenate D_prime, D3, D4 we needed the first item in each tuple to be the same type, which is why Sprint() was used
	D_prime := [];
	for cusp in cusps do
		width := CuspWidth(G,cusp);
		Append(~D_prime,<Sprint(cusp),width-6>);
	end for;
	
	D3 := [];
	D4 := [];
	for e in elliptic_points do
		z := ExactValue(e); // Returns the value of e in its corresponding quadratic imaginary field to find its min poly over Z
		f := MinimalPolynomial(z,Integers());
		D := Discriminant(f);
		if D eq -3 then
			Append(~D3,<Sprint(e),-4>);
		elif D eq -4 then
			Append(~D4,<Sprint(e),-3>);
		else
			error "Unexpected discriminant", D, "for", e;
		end if;
	end for;
	
	sixK := D_prime cat D3 cat D4;
	
	// We now compute FTn

	Doo := [];
	for cusp in cusps do
		Append(~Doo,<Sprint(cusp),2>);
	end for;
	deg_Doo := 0;
	if (#Doo ne 0) then
		deg_Doo := &+[d[2] : d in Doo];
	end if;
	
	// Find CM elliptic curves with endomorphism ring O containing an element of norm ell
	// and an ideal I such that O/I is cyclic of order N

	// Possible discriminants of orders we could have that contain an element of norm n
	candidates := [-D : D in [1..4*n] | NormEquation(D,4*n) and (-D mod 4 in {0,1})];

	FTn_CM := [];
	for D in candidates do
		if (HeegnerCheck(N,D)) then
			mult := 1 + KroneckerSymbol(D,n);
			Append(~FTn_CM,<Sprint(D),mult>);
		end if;
	end for;
	deg_FTn_CM := 0;
	if (#FTn_CM ne 0) then
		deg_FTn_CM := &+[d[2] : d in FTn_CM];
	end if;

	deg_FTn := deg_Doo + deg_FTn_CM;
	
	// NEED: How do we compute Tn(K_X0(N))?
	// For now, we will ignore this term. Cusps map to cusps, so we will only look at when K_X0(N) has no non-cuspidal points in its support.

	non_cusps := DivisorScaling(FTn_CM,6*(2*g-2));
	
	if (#D3 eq 0) and (#D4 eq 0) then
		apply := true;
	else
		apply := false;
	end if;

    return apply, non_cusps;
end function;
