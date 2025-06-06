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
    if alexander_nullity(L) != n - 1:
        return False

    # Check polynomial condition TODO
    a_poly = alexander_polynomial(M, n - 1)
    return True

#helper function for eisermann
def signed_determinant(K):
    V = K.seifert_matrix()
    return (-I * (V + V.transpose())).determinant()

def eisermann(L):
    n = len(L.link_components)
    jones_L = L.jones_polynomial()
    q = jones_L.parent().gen()
    jones_unlink = (q + q^(-1))^(n-1)
    if not jones_unlink.divides(jones_L):
        return False
    quotient = jones_L / jones_unlink
    determinant_product = 1
    for component in L.link_components:
        K = L.sublink(component)
        if (len(K.link_components) != 0):
            determinant_product *= real(signed_determinant(K))
    quotientI = quotient(I)
    if (quotientI.imag_part() == 0 and quotientI.real() % 32 == determinant_product % 32):
        return True
    return False

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
    L = L.split_link_diagram()[0]

def test(L, writer):
    a = sig_zero(L)
    if not a:
        return

    b = linking_matrix_zero(L)
    if not b:
        return

    c = fox_milnor(L)
    if not c:
        return

    d = eisermann(L)

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
