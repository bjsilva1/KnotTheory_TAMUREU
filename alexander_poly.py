import snappy
import string
import sage
from sage.rings.fraction_field import FractionField
from sage.all import Integers, vector, matrix, LaurentPolynomialRing, PolynomialRing, gcd, prod

#----------------------------------------------------------------
#
#  Abelianization of the fundamental group
#
#----------------------------------------------------------------

ZZ = Integers()

def abelianize_word(word, gens):
    return vector(ZZ, [ word.count(g) - word.count(g.swapcase()) for g in gens])

class MapToFreeAbelianization(sage.structure.sage_object.SageObject):
    def __init__(self, fund_group):
        self.domain_gens = fund_group.generators()
        R = matrix(ZZ, [abelianize_word(R, self.domain_gens) for R in fund_group.relators()]).transpose()
        D, U, V = R.smith_form()
        self.U = U
        self.elementry_divisors = [D[i,i] for i in range(D.ncols())] + [0,]*(D.nrows() - D.ncols())

    def range(self):
        return ZZ**self.elementry_divisors.count(0)

    def __call__(self, word):
        D = self.elementry_divisors
        v = self.U*abelianize_word(word, self.domain_gens)
        return vector(ZZ, [v[i] for i in range(len(D)) if D[i] == 0])


#--------------------------------------------------------------
#
#  Computing the Alexander polynomial
#
#--------------------------------------------------------------

class MapToGroupRingOfFreeAbelianization(MapToFreeAbelianization):
    def __init__(self, fund_group):
        MapToFreeAbelianization.__init__(self, fund_group)
        n = self.elementry_divisors.count(0)
        self.R = LaurentPolynomialRing(ZZ, list(string.ascii_lowercase[:n]))

    def range(self):
        return self.R

    def __call__(self, word):
        v = MapToFreeAbelianization.__call__(self, word)
        return prod([ g**v[i] for i, g in enumerate(self.R.gens())])

def fox_derivative(word, phi, var):
    R = phi.range()
    if len(word) == 0:
        return R(0)
    w = word[0]
    if w == w.lower():
        term = R(0) if w != var else R(1)
    else:
        term =   R(0) if w.lower() != var else -phi(var.upper())
    return term + phi(w) * fox_derivative(word[1:], phi, var)

def fox_derivative_r(words, phi, var, gens):
    R = phi.range()
    lookup_table = dict()
    # Fill out lookup table base cases
    for gen in gens:
        lookup_table[gen] = R(0)
        lookup_table[gen.upper()] = R(0)

    lookup_table[var] = R(1)
    lookup_table[var.upper()] = -phi(var.upper())

    output = []
    for word in words:
        def fox_derivative_recursion(i, j):
            if word[i:j] in lookup_table:
                return lookup_table[word[i:j]]
    
            mid = (i + j) // 2
            new_fox = fox_derivative_recursion(i, mid) + phi(word[i:mid]) * fox_derivative_recursion(mid, j)
            if j - i < 10:
                 lookup_table[word[i:j]] = new_fox
            return new_fox

        output.append(fox_derivative_recursion(0, len(word)))
    
    return output

# It's best to deal with matrixes of polynomials rather than Laurent
# polynomials, so we need to be able to clear denominators.  This add
# to the complexity of the code.  

def join_lists(list_of_lists):
    ans = []
    for L in list_of_lists:
        ans += L
    return ans

def uniform_poly_exponents(poly):
    return [list(e) if hasattr(e, "__getitem__") else (e,) for e in poly.exponents()]
    
def minimum_exponents(elts):
    A =  matrix(ZZ, join_lists([ uniform_poly_exponents(p) for p in elts]))
    return vector( [min(row) for row in A.transpose()] )

def convert_laurent_to_poly(elt, expshift, P):
   if elt == 0:
       return P(0)
   return sum( [  c*prod([g**e for g, e in zip(P.gens(), vector(exps) + expshift)]) for c, exps in zip(elt.coefficients(), uniform_poly_exponents(elt))])

def alexander_presentation(knot):
    manifold = knot.exterior()
    G = manifold.fundamental_group()
    
    # Fox derivative calculation of Alexandar Module
    phi = MapToGroupRingOfFreeAbelianization(G)
    R = phi.range()
    P = PolynomialRing(R.base_ring(), list(R.gens_dict().keys()))
    M = [fox_derivative_r(G.relators(), phi, var, G.generators()) for  var in G.generators()]
    #M = [[fox_derivative(rel, phi, var)  for rel in G.relators()] for  var in G.generators()]
    expshift = -minimum_exponents(join_lists(M))
    
    return matrix(P, [[ convert_laurent_to_poly(p, expshift, P) for p in row] for row in M])


def alexander_nullity(M):
    # Compute Alexandar Nullity
    n = M.nrows()
    F = M[0][0].parent().fraction_field()
    M_frac = M.change_ring(F)
    return n - M_frac.rank() - 1

def alexander_polynomial(M, d):
    n = M.nrows()
    alex_poly = gcd(M.minors(n - d - 1))
    
    # Normalize Alexandar Polynomial
    return convert_laurent_to_poly(alex_poly, -minimum_exponents( [alex_poly] ), M[0][0].parent())
