#!/usr/bin/env python3
"""Exact checks for the displayed algebra, using only the Python standard library.

Lower-case letters denote generators; upper-case letters their inverses.
In these encodings a and b denote the manuscript's A and B.
No quotient search or group decision is performed. The geometric proofs
and the HNN normal-form arguments are not assertions of this program.
"""
from fractions import Fraction
import json


def red(word):
    out = []
    for letter in word:
        if out and out[-1] == letter.swapcase():
            out.pop()
        else:
            out.append(letter)
    return ''.join(out)


def inv(word):
    return word[::-1].swapcase()


def mul(*words):
    return red(''.join(words))


def conj(a, b):
    return mul(a, b, inv(a))


def comm(a, b):
    return mul(a, b, inv(a), inv(b))


def subst(word, images):
    return mul(*(images[t] if t.islower() else inv(images[t.lower()])
                 for t in word))


def exponents(word, alphabet):
    return [word.count(t) - word.count(t.upper()) for t in alphabet]


def rank(matrix):
    m = [list(map(Fraction, row)) for row in matrix]
    k = 0
    for col in range(len(m[0])):
        pivot = next((i for i in range(k, len(m)) if m[i][col]), None)
        if pivot is None:
            continue
        m[k], m[pivot] = m[pivot], m[k]
        lead = m[k][col]
        m[k] = [x / lead for x in m[k]]
        for i in range(len(m)):
            if i != k:
                lead = m[i][col]
                m[i] = [x - lead*y for x, y in zip(m[i], m[k])]
        k += 1
        if k == len(m):
            break
    return k


def determinant(matrix):
    m = [list(map(Fraction, row)) for row in matrix]
    result = Fraction(1)
    for col in range(len(m)):
        pivot = next((i for i in range(col, len(m)) if m[i][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            m[col], m[pivot] = m[pivot], m[col]
            result = -result
        lead = m[col][col]
        result *= lead
        for i in range(col + 1, len(m)):
            q = m[i][col] / lead
            m[i] = [x - q*y for x, y in zip(m[i], m[col])]
    return result


def matmul(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)]
            for row in a]


def main():
    # The fifteen geometric relations, independently transcribed.
    core, e, k0, ell = inv('rx'), 'sRS', 'Rmr', 'ax'
    d = mul('b', inv(k0))
    original = {
        'surface': mul(comm('x', 'y'), comm('r', 's')),
        'MC': comm('m', core), 'NE': comm('n', e),
        'Ax': 'axAR', 'Ay': 'ayAS', 'Ar': 'arAX', 'As': 'aNsAY',
        'Bx': 'bxBy', 'By': 'bMyBXY', 'Br': 'brBR',
        'Bs': mul('b', inv(k0), 'sBS'),
        'alpha_meridian': mul('amA', inv(k0)),
        'beta_meridian': mul('b', inv(k0), 'n', k0, 'BN'),
        'alpha_filling': mul('m', ell),
        'beta_filling': mul('n', inv(d)),
    }
    assert len(original) == 15

    r, s, z, j = 'axA', 'ayA', 'aa', comm('x', 'y')
    m = inv(ell)
    n = comm(s, inv(z))
    b = conj(inv(s), n)
    images = dict(a='a', x='x', y='y', r=r, s=s, m=m, n=n, b=b)
    assert b == 'AYaayA'
    mapped = {name: subst(word, images) for name, word in original.items()}
    R, T = mul(j, conj('a', j)), comm(z, 'x')
    X = mul(conj(b, 'x'), 'y')
    Y = mul(b, 'axy', inv(b), 'XY')
    S = mul(b, 'xa', s, inv(b), inv(s))
    k = conj(inv(r), m)
    E = conj(s, inv(r))
    D = mul(b, inv(k))
    S0 = mul(b, inv(k), s, inv(b), inv(s))

    # Reversible beta row replacement before B is eliminated.
    pre = images | {'b': 'b'}
    beta_pre = subst(original['beta_filling'], pre)
    bs_pre = subst(original['Bs'], pre)
    replacement = mul(inv(s), beta_pre, bs_pre, s)
    assert replacement == mul(b, 'B')
    assert beta_pre == mul(s, replacement, inv(s), inv(bs_pre))

    g = mul(inv(ell), 'X')
    bstar = conj(s, inv(z))
    u = mul(bstar, inv(r), 'a')
    v = mul(s, inv(z), inv(r), 'a')
    ne = mul(conj(u, inv(T)), conj(u, inv(R)),
             conj(mul(u, j, 'a'), R), conj(v, T))
    w = mul(inv(s), inv(E))
    assert k == mul(inv('xa'), inv(T))
    assert S0 == mul(conj(b, T), S)
    assert S0 == mul(D, inv(n))
    expected = {
        'surface': R, 'MC': conj(g, T), 'NE': ne,
        'Ax': '', 'Ay': '', 'Ar': T, 'As': '', 'Bx': X, 'By': Y,
        'Br': conj(w, inv(ne)), 'Bs': S0, 'alpha_meridian': '',
        'beta_meridian': mul(S0, conj(n, inv(S0))),
        'alpha_filling': '', 'beta_filling': inv(S0),
    }
    assert mapped == expected, [k for k in mapped if mapped[k] != expected[k]]
    # Reverse products in the mapped original rows.
    assert R == mapped['surface'] and T == mapped['Ar']
    assert X == mapped['Bx'] and Y == mapped['By']
    assert S == mul(conj(b, inv(mapped['Ar'])), mapped['Bs'])

    # Free surface-transport identity, with all original generators free.
    f0 = original['surface']
    left = mul(comm('x', 'My'), comm('r', mul(inv(k0), 's')))
    right = mul(conj('xMX', f0), conj('xM', inv(original['MC'])))
    assert left == right

    # Both directions of the last two-row replacement, before substitutions.
    xrow, yrow, srow = 'bxBy', 'baxyBXY', 'bxasBS'
    crow, astar, little_m = 'bYXsBSx', 'ByxbYX', 'yBXYb'
    lam = mul('a', inv(astar))
    assert little_m == inv(mul(astar, 'x'))
    assert lam == mul('ax', little_m)
    assert red(yrow) == conj('b', lam)
    assert red(srow) == mul(conj('bxB', yrow), xrow, conj('x', crow))
    assert conj('x', crow) == mul(inv(xrow), inv(conj('bxB', yrow)), srow)

    identity = {t: t for t in 'axyrsbmn'}
    for suffix in ('b', 'B'):
        theta = identity | {'s': 's'+suffix}
        reverse = identity | {'s': 's'+inv(suffix)}
        assert subst(crow, theta) == crow
        assert subst(xrow, theta) == xrow
        assert subst(f0, theta) == mul(f0, conj('s', comm('r', suffix)))
        for t in identity:
            assert subst(subst(t, theta), reverse) == t

    # Integral homology: a unimodular eight-row minor proves H1(V)=0.
    rows = ['Bs', 'By', 'Bx', 'Ax', 'Ay', 'As',
            'alpha_filling', 'beta_filling']
    minor = determinant([exponents(original[k], 'xyrsabmn') for k in rows])
    assert abs(minor) == 1
    source = [exponents(t, 'xyrs') for t in ('r', 'x', 'YXs')]
    target = [exponents(t, 'xyrs') for t in ('r', 'Y', 'Xs')]
    assert rank(source) == rank(target) == 3
    assert exponents(crow, 'xyrsb') == [0, -1, 0, 0, 0]
    assert exponents(mul(astar, 'x', inv(astar), 'R'), 'rsb') == [-1, 0, 0]
    assert exponents(lam, 'a') == [1]

    # Normal coordinates use (t, kappa*(w-1/2)) or (t,kappa*(v-1/2)).
    normal_alpha = [[-1, 0], [Fraction(1, 2), -1]]
    normal_beta = [[-1, 0], [0, 1]]
    local_reflection = [[1, 0], [0, -1]]
    reflected_alpha = [[1, 0], [-1, -1]]
    reflected_beta = [[1, 0], [0, -1]]
    assert determinant(normal_alpha) == 1 and determinant(normal_beta) == -1
    for chart, reflection in ((normal_alpha, reflected_alpha),
                              (normal_beta, reflected_beta)):
        assert matmul(chart, local_reflection) == matmul(reflection, chart)
        assert determinant(reflection) == -1
        assert matmul(reflection, reflection) == [[1, 0], [0, 1]]

    print(json.dumps(dict(
        verification='PASS', original_relations=15,
        forward_original_row_checks=len(mapped), reverse_retained_row_checks=5,
        reversible_beta_elimination=True, exact_final_two_row_replacement=True,
        theta_and_inverse_checked=True, abelianization_minor=int(minor),
        rank_three_source=source, rank_three_target=target,
        cotangent_meridian_signs=[1, -1], normal_reflections_checked=2,
        geometric_inclusions_proved_by_program=False,
        full_group_decided=False, J_in_final_group='unknown'), indent=2))


if __name__ == '__main__':
    main()
