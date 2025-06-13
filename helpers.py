import csv
import snappy
from snappy import *
from alexander_poly import alexander_presentation, alexander_nullity, alexander_polynomial

def sig_zero(L):
    return (L.signature() == 0)

def linking_matrix_zero(L):
    output = True
    M = L.linking_matrix()
    for i in M:
        for j in i:
            if (j != 0):
                output = False
    return output

def fox_milnor(L):
    n = len(L.link_components)
    # Check Alexander nullity
    M = alexander_presentation(L)
    if alexander_nullity(M) != n - 1:
        return (False, False)

    # Check polynomial condition TODO
    a_poly = alexander_polynomial(M, n - 1)
    return (True, is_SplitPoly(a_poly))

#helper function for eisermann
def signed_determinant(K):
    V = K.seifert_matrix()
    return (-I * (V + V.transpose())).determinant()

#returns a pair of booleans. The first one is true if the divisibility condition of eisermann is satisfied,
#and the second one is true if the mod 32 condition is satisfied.
def eisermann(L):
    n = len(L.link_components)
    jones_L = L.jones_polynomial()
    q = jones_L.parent().gen()
    jones_unlink = (q + q**(-1))**(n-1)
    if not jones_unlink.divides(jones_L):
        print("no divides")
        return (False, False)
    quotient = jones_L / jones_unlink
    determinant_product = 1
    for component in L.link_components:
        K = L.sublink(component)
        if (len(K.link_components) != 0):
            determinant_product *= real(signed_determinant(K))
    quotientI = quotient(I)
    if (quotientI.imag_part() == 0 and quotientI.real() % 32 == determinant_product % 32):
        return (True, True)
    return (True, False)

def zeroify_simple(L):
    crossings = L.crossings
    link_mat = L.linking_matrix()
    for crossing in crossings:
        component1 = crossing.strand_components[0]
        component2 = crossing.strand_components[1]
        if (link_mat[component1][component2] > 0 and crossing.sign == 1):
            crossing.rotate_by_90()
            link_mat[component1][component2] -= 1
            link_mat[component2][component1] -= 1
        if (link_mat[component1][component2] < 0 and crossing.sign == -1):
            crossing.rotate_by_90()
            link_mat[component1][component2] += 1
            link_mat[component2][component1] += 1
    L._rebuild()
    L.simplify("global")
    if (len(L.link_components) >= 1):
        X = L.split_link_diagram()[0]
        L = X

# helper for zeroify_twist
# modifies the link diagram of L by giving a half-twist to the specified crossing.
# the orientation of this half-twist depends on whether "direction" is 1 or -1.
# after this function is run, L must be rebuilt.
def half_twist_crossing(L, crossing, direction):
    if (crossing.sign == 1):
        new_crossing = Crossing()
        for i in [1,2]:
            new_crossing[i] = crossing.adjacent[i]
        crossing[1] = new_crossing[0]
        crossing[2] = new_crossing[3]
        if (direction == -1):
            new_crossing.rotate_by_90()
    if (crossing.sign == -1):
        new_crossing = Crossing()
        for i in [3,2]:
            new_crossing[i] = crossing.adjacent[i]
        crossing[3] = new_crossing[0]
        crossing[2] = new_crossing[1]
        if (direction == 1):
            new_crossing.rotate_by_90()
    L.crossings.append(new_crossing)
def twist_crossing(L, crossing, direction):
    half_twist_crossing(L, crossing, direction)
    half_twist_crossing(L, crossing, direction)

def zeroify_twist(L):
    crossings = L.crossings
    link_mat = L.linking_matrix()
    num_crossings = len(crossings)
    for i in range(num_crossings):
        crossing = crossings[i]
        component1 = crossing.strand_components[0]
        component2 = crossing.strand_components[1]
        if (link_mat[component1][component2] > 0 and crossing.sign == -1):
            for i in range(link_mat[component1][component2]):
                twist_crossing(L, crossing, -1)
            link_mat[component1][component2] = 0
            link_mat[component2][component1] = 0
        if (link_mat[component1][component2] < 0 and crossing.sign == 1):
            for i in range(-link_mat[component1][component2]):
                twist_crossing(L, crossing, 1)
            link_mat[component1][component2] = 0
            link_mat[component2][component1] = 0
    for i in range(num_crossings):
        crossing = crossings[i]
        component1 = crossing.strand_components[0]
        component2 = crossing.strand_components[1]
        if (link_mat[component1][component2] > 0):
            for i in range(link_mat[component1][component2]):
                twist_crossing(L, crossing, -1)
            link_mat[component1][component2] = 0
            link_mat[component2][component1] = 0
        if (link_mat[component1][component2] < 0):
            for i in range(-link_mat[component1][component2]):
                twist_crossing(L, crossing, 1)
            link_mat[component1][component2] = 0
            link_mat[component2][component1] = 0
    L._rebuild()
    L.simplify("global")
    if (len(L.link_components) >= 1):
        X = L.split_link_diagram()[0]
        L = X

def test(L, writer):
    a = sig_zero(L)
    if not a:
        return

    b = linking_matrix_zero(L)
    if not b:
        return

    d = eisermann(L)
    if d[1]:
        return

    c = fox_milnor(L)

    if c[0] or (len(L.link_components) == 2 and (not d[0])):
        print("hit! ", len(L.link_components), len(L.crossings))
        writer.writerow([
            L.PD_code(), 
            len(L.link_components),
            len(L.crossings), 
            a, 
            b, 
            c, 
            d
        ])

def is_SplitPoly(poly1):
    #Check if poly2 has even degree
    variables = poly1.parent().gens()
    variableArray = []
    for var in variables:
        variableArray.append(var)
    for n in range(len(variableArray)):
        if (poly1.degree(variableArray[n]) % 2 != 0):
            return False

    #Check if poly2 has even number of factors
    factor_list = []
    for factor, multiplicity in poly1.factor():
        for i in range(multiplicity): 
            factor_list.append(factor)
    if (len(factor_list) % 2 != 0):
        return False

    #Cancel factors that satisfy Fox-Milnor
    for i in range(len(factor_list)): 
        if factor_list[i] == 'y':
            continue

        factor1 = factor_list[i]
        factor2 = factor1
        inverseVariables = []
        for v in range(len(variableArray)):
            inverseVariables.append(variableArray[v]**(-1))
        factor2 = factor2(inverseVariables)
        degrees = []
        for n in range(len(variableArray)):
            degrees.append(factor1.degree(variableArray[n]))
        for d in range(len(degrees)):
            factor2 = factor2 * (variableArray[d]**degrees[d])
        matched = False
        for j in range(i+1,len(factor_list)):
            if (factor_list[j] != 'y' and (factor2 == factor_list[j] or factor2 == -(factor_list[j]))):
                factor_list[i] = 'y'
                factor_list[j] = 'y'
                matched = True
                break
        if not matched:
            return False
    # If all factors cancel
    return True

def nullity(L):
    A = L.seifert_matrix()
    M = A + A.transpose()
    return M.nullity()

def jones_nullity(L):
    jones_L = L.jones_polynomial()
    q = jones_L.parent().gen()
    jones_unknot = (q + q**(-1))
    output = 0
    while jones_unknot.divides(jones_L):
        output += 1
        jones_L = jones_L / jones_unknot
    return output
