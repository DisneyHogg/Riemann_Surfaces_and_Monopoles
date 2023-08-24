from sys import displayhook
import warnings


def point_divide(A, B):
    r"""
    Given matrices A, B, if A = k*B for some k, return k. Otherwise raise an error.
    """
    if B.is_zero():
        raise ZeroDivisionError()
    
    for ai, bi in zip(A, B):
        if not bi:
            continue
        k = ai/bi
        break
        
    if not all(ai==k*bi for ai, bi in zip(A, B)):
        # return None
        raise ValueError("matrices aren't proportional to each other")
    
    return k 


def polarised_to_matrices(XY, qs, exponent, sign=1):
    r"""
    Run steps 3 and 4 of the upper-invariant algorithm.
    
    TODO: Write this out vectorially. 
    """
    ###
    # Step 3
    # Set up the objects needed for the calculation
    X, Y = XY
    nr = X.nrows()
    SM = X**ZZ(exponent)
    adY = lambda M: sign*Y.commutator(M)
    Ssp = []
    # Loop over the polynomials
    for qi in qs:
        # If qi is the zero polynomial, know the corresponding matrix is 0
        if not qi:
            Ssp.append(matrix.zero(nr))
            continue 
        # If not, iteratively apply the commutator, and sum.
        # Note this iterative approach is why we don't do the sum in a 
        # sum() function. 
        M = SM
        # ql = [co for co, _ in list(qi)]
        ql = qi.list()
        S = ql[0]*M
        for co in ql[1:]:
            M = adY(M)
            S += co*M
        Ssp.append(S)
    ###
    
    ###
    # Step 4
    K = qs[0].base_ring()
    j = (polygen(K)**2 + 1).roots(multiplicities=False)[0]
    Ss = [Ssp[0]/2 + Ssp[2], -j*Ssp[0]/2 + j*Ssp[2], sign*j*Ssp[1]]
    
    return [Ssi.change_ring(K) for Ssi in Ss]


def invariant_algorithm(XY, Q, level, sign=1):
    r"""
    Return the invariant vector corresponding to a polynomial at level $\in \lbrace -1, 0, 1 \rbrace$.
    """
    # We will assume this algorithm has only ever been called internally, and so not apply sense checks
    # to the inputs XY, Q, or level. 
    
    # Record some key variables we will reuse.
    deg = Q.degree()
    K = Q.base_ring()
    d = deg + 2*level
    
    # Polarise Q
    z0, z1 = Q.parent().gens()
    Q00, Q01, Q11 = Q.derivative(z0, 2), Q.derivative(z0, z1), Q.derivative(z1, 2)
    Polar = [Q11, Q01, Q00/2]
    
    if level == -1:
        # In the level = -1 case, we don't need to use the highest weight vector trick,
        # and can immediately convert the S^d component into component-wise polynomials
        # of adY at the polynomial level. 
        
        # Maybe eventually rewrite this for conceptual simplicity. 
        Kx = PolynomialRing(K, name='x')
        x = Kx.gen()
        def convert_zw_to_x(p):   
            cd = p.dict()
            return Kx(sum(cd[k]*x**k[0]*factorial(k[1]) for k in cd.keys()))/factorial(d)
        # Note the implicit assumption made about our ring in the above calculation, that factorial(d)
        # is a unit.
        Qts_11_01_00 = tuple(convert_zw_to_x(qi) for qi in Polar)
        
    else:
        
        ###
        # CONSTRUCT POLYNOMIAL QT MATCHING POLAR TO HIGHEST WEIGHT VECTOR
        ###
        
        # For the first S2 = 3 factor of the tensor product we will use the basis
        # z1^2, 2*z0*z1, 2*z0^2, to make computations simpler. Note choosing this 
        # basis is at the expense of symmetry. 
        adY2 = matrix([[0, 0, 0], [1, 0, 0], [0, 1, 0]])

        # This is d at level -1
        dm1 = deg - 2
        dim = dm1 + 1
        
        # In the Sd = d+1 factor we are working with the obvious basis of e_a = z0^a * z1^(d-a)
        adYc = matrix.zero(dim)
        for i in range(dm1):
            adYc[i+1, i] = dm1-i
        # adY is the tensor product of the reps S2 and S(deg-2)
        adY = adY2.tensor_product(matrix.identity(dim)) + matrix.identity(3).tensor_product(adYc)

        # Create the highest weight vector in the tensor product space
        td = matrix.zero(dim, 1)
        td[0] = 1 
        hwv = vector(matrix(3, 1, [1, 0, 0]).tensor_product(td))

        # Use this vector to generate a basis of the S(d+2)=S(deg) subspace which Q will lie in. 
        # Because this is a highest-weight subspace, we can do so with repeated applications of 
        # adY.
        basis = [hwv]
        bk = hwv
        for i in range(deg):
            bk = adY*bk
            basis.append(bk)
            
        V = VectorSpace(K, 3*dim)
        W = V.subspace_with_basis(basis)

        # Find the vector in S2 \otimes Sd corresponding to Q
        def pol2vec(p):
            # converts a polynomial in Sd to a vector (really a (d+1)x1 matrix).
            cd = p.dict()
            q = matrix.zero(K, dim, 1)
            for k in cd.keys():
                q[k[0]] = cd[k]
                nd = k[0]
            return q
        Q_vec = V(0)
        for i, qi in enumerate(Polar):
            comp = matrix.zero(K, 3, 1)
            comp[i] = 1
            Q_vec += vector(comp.tensor_product(pol2vec(qi)))

        # Rewrite the vector corresponding to Q in the basis of W = S(deg)
        Q_vec = W(Q_vec)
        c = W.coordinate_vector(Q_vec)

        # Use the coordinate vector to get the corresponding tilde polynomial
        Kx = PolynomialRing(K, name='x')
        x = Kx.gen()
        Qt = sum(co*x**i for i, co in enumerate(c))
        
        ###
        # ACT WITH QT ON HIGHEST WEIGHT VECTOR IN S2 \otimes Sd
        ###
        dim = d + 1
        adYc = matrix.zero(dim)
        for i in range(d):
            adYc[i+1, i] = d-i
        adY = adY2.tensor_product(matrix.identity(dim)) + matrix.identity(3).tensor_product(adYc)
        adY = adY.change_ring(K)

        # Build the new heighest weight vector
        # Perhaps in the future this can be made more streamlined using the intersections-of-kernels
        # definition. See cells below. 
        if level == 0:
            td = matrix.zero(dim, 1)
            td[0] = 1 
            hwv = vector(matrix(3, 1, [0, 1, 0]).tensor_product(td))
            td = matrix.zero(dim, 1)
            td[1] = 1 
            hwv += vector(matrix(3, 1, [-2, 0, 0]).tensor_product(td))
        elif level == 1:
            td = matrix.zero(dim, 1)
            td[0] = 1 
            hwv = vector(matrix(3, 1, [0, 0, 1]).tensor_product(td))
            td = matrix.zero(dim, 1)
            td[1] = 1 
            hwv += vector(matrix(3, 1, [0, -2, 0]).tensor_product(td))
            td = matrix.zero(dim, 1)
            td[2] = 1 
            hwv += vector(matrix(3, 1, [2, 0, 0]).tensor_product(td))
        else:
            raise ValueError("Invalid level of {} provided".format(level))

        # Apply the polynomial of adY to the highest weight vector
        Q_vec = Qt(adY)*hwv

        # Unpack the vector back into it's polarised components
        # This relies on our knowledge of how the tensor product of matrices orders the basis 
        # vectors of the corresponding tensor of vector spaces. This is in 
        # https://en.wikipedia.org/wiki/Kronecker_product. 
        # We will do the process of calculating the multivariate polynomial component and 
        # converting it into a univariate polynomial in a derivative operator all in one. 
        Qts_11_01_00 = []
        Kx = PolynomialRing(K, name='x')
        x = Kx.gen()
        for i in range(3):
            qi = Kx(sum(Q_vec[j + (d+1)*i]*(x**j)*factorial(d-j) for j in range(d+1)))/factorial(d)
            Qts_11_01_00.append(qi)
        
    return polarised_to_matrices(XY, Qts_11_01_00, d/2, sign=sign)


def find_invariant_vectors(Ms, Qs, sign=1, exclude_level_zero=False):
    r"""
    Given matrices Ms which gives a basis of su(2), either of the form of r1, r2, r3, or X, Y, H,
    and a polynomial of an appropriate corresponding degree, return the invariant vectors. 
    
    Returns rs as a list and Ss as a list.
    """
    # If a single polynomial is given for Qs, instead of a list of length 1, remedy this error.
    if not isinstance(Qs, list):
        if isinstance(Qs, sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular):
            Qs = [Qs]
        elif isinstance(Qs, sage.rings.polynomial.multi_polynomial_element.MPolynomial_polydict):
            Qs = [Qs]
        else:
            raise ValueError("Class of Qs not recognised as standard input")
    
    # Make sure the Qs are defined over a complex field where sqrt(-1) exists, and set this root. 
    pol = polygen(Qs[0].base_ring())**2 + 1
    if pol.is_irreducible():
        raise ValueError("Base ring must contain sqrt(-1)")
    j = pol.roots(multiplicities=False)[0]
    
    # Sort out the matrices
    if all(point_divide(Ms[i].commutator(Ms[(i + 1) % 3]), Ms[(i + 2) % 3]) == 2 
           for i in range(3)):
        r1, r2, r3 = Ms
        X = (r1 - j*r2)/2; Y = -(r1 + j*r2)/2; H = -j*r3
    elif tuple(point_divide(Ms[i].commutator(Ms[l]), Ms[k]) 
          for i, l, k in [[2, 0, 0], [2, 1, 1], [0, 1, 2]]) == (2, -2, 1):
        X, Y, H = Ms
        r1 = X - Y; r2 = j*(X + Y); r3 = j*H
    else:
        raise ValueError("Valid triple of matrices needs to be provided.")
    
    # The charge will correspond to the dimension of the matrices. 
    charge = r1.nrows()

    # Loop over the polynomials and get the corresponding invariant vectors.         
    S_invariants = []
    for Q in Qs:
        # If for any reason a zero polynomial has been provided, just skip over it
        if not Q:
            continue
            
        if not Q.is_homogeneous():
            raise ValueError("Polynomial {} is not homogeneous".format(Q))
            
        if Q.degree() % 2:
            raise ValueError("Polynomial {} is of odd degree".format(Q))
            
        if Q.degree() < 2:
            raise ValueError("Polynomial {} is of degree < 2".format(Q))
            
        top = ZZ(min((2*charge - Q.degree())/2, 2))
        for lvl in range(-1, top):
            if exclude_level_zero and lvl==0:
                continue

            S_invariants.append(invariant_algorithm([X, Y], Q, lvl, sign=sign))
        
    return [r1, r2, r3], S_invariants


def echelon_form_with_transformation(original_matrix):
    r"""
    Return a tuple (E, T) where E is the (reduced) row echelon form of M and T is the transformation s.t
    T*M = E. 
    
    This replaces Sage's '.echelon_form(transformation=True)', which sometimes quitely ignores the 
    'transformation' argument. 
    """
    # We don't use too much brain power, and use the pseudocode algorithm from wikipedia,
    # https://en.wikipedia.org/wiki/Row_echelon_form
    
    # We make two copies, one so that the argument matrix is not altered, 
    # and one for a sense check at the end. 
    M = copy(original_matrix)
    test_copy = copy(M)
    
    nr, nc = M.dimensions()
    T = matrix.identity(nr)
    
    lead = 0
    for r in range(nr):
        if nc <= lead:
            return M, T
        i = r 
        while M[i, lead] == 0:
            i += 1
            if nr == i:
                i = r
                lead += 1
                if nc == lead:
                    return M, T
        if i != r:
            M.swap_rows(i, r)
            T.swap_rows(i, r)
            
        # Warning, there was an error here in previous iterations. 
        scale = 1/M[r, lead]
        M = M.with_rescaled_row(r, scale)
        T = T.with_rescaled_row(r, scale)
        for j in range(nr):
            if j != r:
                scale = -M[j, lead]
                M = M.with_added_multiple_of_row(j, r, scale)
                T = T.with_added_multiple_of_row(j, r, scale)
        lead += 1
        
    if T*test_copy != M:
        raise ValueError("Algorithm has failed to give valid pair")

    return M, T


def solve_commutation_relations(rs, Ss, verbose=False):
    r"""
    Given invariant vectors rs and Ss (a list), look for solutions to the constants
     - aj, j = 1 .. d,
     - bjk, j,k = 1 .. d, 
     - cjk, j = 1 .. d, k = 1 .. j
     - djkl j,l = 1 .. d, k = 1 .. j
    to give the correct correlation relations. Return these as a tuple of 4 dicts, as well as a list
    of the necessary constraints on the variables for there to be a solution. 
    
    Note that solutions to Nahm data may exist even if these commutators have no solutions, but it will
    require fine-tuning of the scales of at least one invariant vector. 
    
    This method uses a different approach to previously, linearising the problem and fining the solution
    using vectorisation and linear algebra. 
    """
    # If a single invariant vector is given for Ss, instead of a list of length 1, remedy this error.
    if not isinstance(Ss[0], list):
        if isinstance(Ss[0], sage.matrix.matrix_generic_dense.Matrix_generic_dense):
            Ss = [Ss]
        else:
            raise ValueError("Class of Ss not recognised as standard input")
    
    # Define the linear algebra functions we will need.
    # These have been taken to be consistent, and s.t. 
    #
    #  vectorise(A.commutator(B)) = comm_rep(A)*vectorise(B).
    #
    vectorise = lambda M: matrix(M.base_ring(), M.nrows()*M.ncols(), 1, M.transpose().list())
    charge = rs[0].nrows()
    K = rs[0].base_ring()
    comm_rep = lambda M: matrix.identity(K, charge).tensor_product(M) - M.transpose().tensor_product(
                         matrix.identity(K, charge))
    
    
    # Initialise the rep and vectorisation of the rs and the Ss
    prs = [comm_rep(ri) for ri in rs]
    vrs = [vectorise(ri) for ri in rs]

    pSs = [[comm_rep(sij) for sij in si] for si in Ss]
    vSs = [[vectorise(sij) for sij in si] for si in Ss]
    
    # Calculate the vectorisation of the commutators that occur in the identities.
    v_rs_comms = [[prs[p]*vsi[(p+1)%3]-prs[(p+1)%3]*vsi[p] for p in range(3)] for vsi in vSs]
    
    d = len(Ss)
    v_ss_comms = [[pSs[j][p]*vSs[k][(p+1)%3]+pSs[k][p]*vSs[j][(p+1)%3] for p in range(3)] 
                  for j in range(d) for k in range(j+1)]
    
    # We concatenate together the vectors based off of their index in rs of the Ss[i].
    # This comes from the fact that we get three sets of commutation relations cycling the index.
    # These will be confusing the read because of the annoying functionality of block_matrix. 
    cat_vrs = block_matrix([[vrs[(p + 2) % 3]] for p in range(3)])
    cat_vSs = [block_matrix([[vsi[(p + 2) % 3]] for p in range(3)]) for vsi in vSs]
    cat_v_rs_comms = [block_matrix([[vrscij] for vrscij in vrsci]) for vrsci in v_rs_comms]
    cat_v_ss_comms = [block_matrix([[vsscij] for vsscij in vssci]) for vssci in v_ss_comms]
    
    # Initialise some dictionaries to store the coefficient solutions we find. 
    a_dict = dict()
    b_dict = dict()
    c_dict = dict()
    d_dict = dict()

    # Somewhat remarkably, the solutions end up coming from affine equations with 
    # only one matrix, which we calculate now. 
    M = block_matrix([cat_vrs]+[cvs for cvs in cat_vSs], ncols=d+1)

    # We now loop through the i and (i,j) pairs that occur in the equations, and 
    # look for solutions.
    # If there are no solutions, the matrix solving will raise the error we need. 
    
    # Set a flag that will capture if the commutation relations have no solutions. 
    Failed = False
    
    for j in range(d):
        # The a-b relations come from the r-s commutators
        Xabj = vector(cat_v_rs_comms[j])
        try:
            ab_sol = M.solve_right(Xabj)
            # Store the results in the corresponding dicts. 
            a_dict.update({j: ab_sol[0]})
            b_dict.update({(j,k): ab_sol[k+1] for k in range(d)})
        except ValueError:
            if verbose:
                print("No solution for {} a-b equation".format(j))
            Failed = True

        for k in range(j+1):
            # The c-d relations come from the s-s commutators
            Xcdjk = vector(cat_v_ss_comms[ZZ(j*(j+1)/2+k)])/(1 + (j==k))
            try:
                cd_sol = M.solve_right(Xcdjk)
                c_dict.update({(j,k): cd_sol[0]})
                d_dict.update({(j,k,l): cd_sol[ZZ(l+1)] for l in range(d)})
            except ValueError:
                if verbose:
                    print("no solution for ({},{}) c-d equation".format(j,k))
                Failed = True
    
    # If finding a solution to the commutations the first way failed, we take 
    # a new approach
    if Failed:
        if verbose:
            warnings.warn("Inconsistencies found, fall back to manual echelon implementation")
        # First check if M has full rank. If not, just raise an error.
        if M.rank() != d+1:
            raise NotImplementedError("Method cannot handle M not full rank.") 
            
        # Otherwise, we will find the conditions on the variables required to 
        # make a solution consistent, then solve.
        # This may eventually be the right way to go straight away. 
        K = Ss[0][0].base_ring()
        d = len(Ss)
        xys = ['x']+['y{}'.format(j) for j in range(d)]
        A = PolynomialRing(K, d+1, xys)
        gs = A.gens()
        x, ys = gs[0], gs[1:]
        V = x*sum(ys[j]*vector(cat_v_rs_comms[j]) 
                  for j in range(d)) + sum(ys[j]*ys[k]*vector(cat_v_ss_comms[ZZ(j*(j+1)/2+k)])/(1 + (j==k))
                    for j in range(d) for k in range(j+1))
        E, T = echelon_form_with_transformation(M)
        constraints = set(list(T*V)[d+1:])
        II = A.ideal([A(ci) for ci in constraints])
        constraints = II.groebner_basis()
        
        # We now solve with this restricted matrix
        Eres = E[:d+1]
        Tres = T[:d+1]
        for j in range(d):
            Xabj = vector(cat_v_rs_comms[j])
            ab_sol = Eres.solve_right(Tres*Xabj)
            a_dict.update({j: ab_sol[0]})
            b_dict.update({(j,k): ab_sol[k+1] for k in range(d)})

            for k in range(j+1):
                # The c-d relations come from the s-s commutators
                Xcdjk = vector(cat_v_ss_comms[ZZ(j*(j+1)/2+k)])/(1 + (j==k))
                cd_sol = Eres.solve_right(Tres*Xcdjk)
                c_dict.update({(j,k): cd_sol[0]})
                d_dict.update({(j,k,l): cd_sol[l+1] for l in range(d)})
        
        return (a_dict, b_dict, c_dict, d_dict), constraints

    rn = M.right_nullity()
    if rn:
        print("Nullity of {} in each equation".format(rn))
        total_nullity = rn*d*(d+3)/2 # = rn * (d + d*(d+1)/2)
        print("Total nullity is {}".format(total_nullity))
        
    return (a_dict, b_dict, c_dict, d_dict), []


def nahm_matrices(rs, Ss):
    r"""
    Calculate the Nahm matrices from the invariant vectors, including arbitrary
    parameters. 
    """
    # If a single invariant vector is given for Qs, instead of a list of 
    # length 1, remedy this error.
    if not isinstance(Ss[0], list):
        if isinstance(Ss[0], sage.matrix.matrix_generic_dense.Matrix_generic_dense):
            Ss = [Ss]
        else:
            raise ValueError("Class of Ss not recognised as standard input")
    K = Ss[0][0].base_ring()
    
    # Make sure the Ss are defined over a complex field where sqrt(-1) exists, 
    # and set this root.
    # In practice do the check at the stage where the polynomials are read in 
    # to find_invariant_vectors. 
    pol = polygen(K)**2 + 1
    j = pol.roots(multiplicities=False)[0]
    
    d = len(Ss)
    xys = ['x']+['y{}'.format(j) for j in range(d)]
    A = PolynomialRing(K, d+1, xys)
    gs = A.gens()
    x, ys = gs[0], gs[1:]
    
    Ts = [x*rs[i].change_ring(A) + sum(yj*sj[i] for yj, sj in zip(ys, Ss)) for i in range(3)]

    return Ts


def lax_pair(rs, Ss, convention="BE"):
    r"""
    Calculate the Lax pair corresponding to the Nahm data found. 
    """
    Ts = nahm_matrices(rs, Ss)
    A = Ts[0].base_ring()
    K = A.base_ring()
    # Make sure the Ts are defined over a complex field where sqrt(-1) exists,
    # and set this root. 
    pol = polygen(K)**2 + 1
    j = pol.roots(multiplicities=False)[0]
    
    Awz = PolynomialRing(A, 2, names=['w', 'z'])
    # Comment out the line below if you want the pretty_print copies of w, z to 
    # remain as w, z. 
    Awz._latex_names = [r'\eta', r'\zeta']

    _, z = Awz.gens()
    if convention == "BE":
        L = (Ts[0] + j*Ts[1]) - 2*j*Ts[2]*z + (Ts[0] - j*Ts[1])*z**2
        M = (Ts[0] - j*Ts[1])*z - j*Ts[2]
    elif convention == "HMM":
        L = -j*(Ts[0] + j*Ts[1]) + 2*j*Ts[2]*z + j*(Ts[0] - j*Ts[1])*z**2
        M = None
        warnings.warn("M matrix is not identified in the HMM convention")
    elif convention == "HS":
        L = -(Ts[0] + j*Ts[1]) + 2*j*Ts[2]*z - j*(Ts[0] - j*Ts[1])*z**2
        M = None
        warnings.warn("M matrix is not identified in the HS convention")
    else:
        raise NotImplementedError("Convention for Lax pair not recognised")  
    return L, M


def spectral_curve(rs, Ss, convention="BE"):
    r"""
    Calculate the spectral curve corresponding to the Nahm data found. 
    """
    L, _ = lax_pair(rs, Ss, convention=convention)
    charge = L.nrows()
    Awz = L.base_ring()
    w, _ = Awz.gens()
    M = w*matrix.identity(Awz, charge) - L
    return Awz(M.det())


def ode_system(rs, Ss):
    r"""
    Calculate the ODEs corresponding to Nahm's equations for given Nahm data.  
    
    The ODE's are return as a list, corresponding to the s-deriative of the  
    correspondingly indexed variable. 
    """
    
    # Get the coefficients in dictionary form.
    dicts, _ = solve_commutation_relations(rs, Ss) 
    a_dict, b_dict, c_dict, d_dict = dicts
    
    # If a single invariant vector is given for Qs, instead of a list of 
    # length 1, remedy this error.
    # We do this so we can get the base field for the variables. 
    if not isinstance(Ss[0], list):
        if isinstance(Ss[0], sage.matrix.matrix_generic_dense.Matrix_generic_dense):
            Ss = [Ss]
        else:
            raise ValueError("Class of Ss not recognised as standard input")
    K = Ss[0][0].base_ring()
    d = len(Ss)
    xys = ['x']+['y{}'.format(j) for j in range(d)]
    A = PolynomialRing(K, d+1, xys)
    gs = A.gens()
    x, ys = gs[0], gs[1:]
    
    ODEs = [2*x**2 + x*sum(a_dict[k]*ys[k] for k in range(d)) + sum(c_dict[(k,l)]*ys[k]*ys[l] 
            for k in range(d) for l in range(k+1))] + [x*sum(b_dict[(k,j)]*ys[k] 
            for k in range(d)) + sum(d_dict[(k,l,j)]*ys[k]*ys[l] for k in range(d)
            for l in range(k+1)) for j in range(d)]
    
    return ODEs


def verify_isospectrality(rs, Ss):
    r"""
    Verify that the coefficients of the spectral curve are constants of the evolution. 
    """
    _ , constraints = solve_commutation_relations(rs, Ss)
    f = spectral_curve(rs, Ss)
    ODEs = ode_system(rs, Ss)
    
    BR = f.base_ring()
    if constraints:
        message = "Coefficients are conserved in constrained ring:"
        constrained_ring = BR.quotient(BR.ideal(constraints))
    else:
        message = "Coefficients are conserved in unconstrained ring:"
        constrained_ring = BR
        
    print(message)
    print(not any(constrained_ring(sum(co.derivative(v)*dvds for v, dvds in zip(f.base_ring().gens(), ODEs)))
                  for co in f.coefficients()))
    
    return None


def is_hermitian(M):
    r"""
    Test if matrix is Hermitian. 
    """
    return not any(mi for mi in M.conjugate().transpose() - M)


def is_anti_hermitian(M):
    r"""
    Test if matrix is anti-Hermitian. 
    """
    return not any(mi for mi in M.conjugate().transpose() + M)


def make_hermitian(Ss, K=None):
    r"""
    Go through `Ss`, editing the original list, making each vector contain only
    Hermitian matrices if possible
    """
    if K is None:
        K = Ss[0][0].base_ring()
    j = (polygen(K)**2 + 1).roots(multiplicities=False)[0]
    for Ssi in Ss:
        if all(is_hermitian(Ssij) for Ssij in Ssi):
            for k in range(3):
                Ssi[k] *= j
    print("Succeeded in making all matrices anti-Hermitian:")
    print(all(is_anti_hermitian(Ssij) for Ssi in Ss for Ssij in Ssi))
    return None

