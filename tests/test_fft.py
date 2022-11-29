from pypairing import ZR, G1, blsfft, blsmultiexp, robustblsfft
from adkg.polynomial import polynomials_over, fnt_decode_step1, fnt_decode_step2
import random
from random import shuffle

def get_omega(field, n, seed=None):
    """
    Given a field, this method returns an n^th root of unity.
    If the seed is not None then this method will return the
    same n'th root of unity for every run with the same seed

    This only makes sense if n is a power of 2!
    """
    assert n & n - 1 == 0, "n must be a power of 2"
    x = field.rand(seed)
    y = x**int(field(-1)/n)
    if y == 1 or y**(n//2) == 1:
        return get_omega(field, n)
    assert y**n == 1, "omega must be 2n'th root of unity"
    assert y**(n // 2) != 1, "omega must be primitive 2n'th root of unity"
    return y


def test_fft2():
    t = 1
    logq = 5
    n = 3*t + 1
    omega2 = get_omega(ZR, 2*n)
    omega = omega2**2
    g = G1.rand(b'')

    p = polynomials_over(ZR)
    q = 2**logq
    m = q//(t+1)
    m = 2
    omega = omega2**2

    zs = list(range(n))
    shuffle(zs)
    zs = zs[:t+1]
    polys = [p.random(t) for _ in range(m)]
    all_ys = [polys[i].evaluate_fft(omega, n) for i in range(m)]
    ys = []
    for i in range(m):
        for j in zs:
            ys.append(g**all_ys[i][j])

    fft_rep = robustblsfft(zs, ys, omega2, n)
    ccs = [[g**val for val in polys[i]] for i in range(m)]
    assert fft_rep == ccs

def test_fft():
    for i in range(5): 
        n = 2**(i+2)
        t = random.randint(1, n-1)
        g = G1.rand(b'g')

        omega = get_omega(ZR, n)

        p = polynomials_over(ZR)
        poly = p.random(t)

        eval_scalar = poly.evaluate_fft(omega, n)
        evals = [g**eval for eval in eval_scalar]

        coeffs = [g**x for x in poly.coeffs]
        fft_rep = blsfft(coeffs, omega, n)

        assert fft_rep == evals
        
        fft_inv = blsfft(evals, omega**(-1), n)
        n_zr = ZR(n)
        res = [x**(n_zr**(-1)) for x in fft_inv]

        coeffs += [g**0]*(n-len(coeffs))
        assert res == coeffs

def test_batch_fft():
    for i in range(5): 
        n = 2**(i+2)
        t = random.randint(1, n-1)
        g = G1.rand(b'g')

        omega2 = get_omega(ZR, 2 * n)
        omega = omega2 ** 2

        p = polynomials_over(ZR)
        poly = p.random(t)

        zs = list(range(n))
        shuffle(zs)
        zs = zs[:t+1]
        ys = list(poly.evaluate_fft(omega, n))
        evals = [g**ys[i] for i in range(n)]
        ys_ = [ys[i] for i in zs]
        ys = [g**ys[i] for i in zs]
        ys = ys + ys

        coeffs = [g**x for x in poly.coeffs]
        fft_rep = robustblsfft(zs, ys, omega2, n)

        # as_, ais_ = fnt_decode_step1(p, zs, omega2, n)

        # fft_inv = blsfft(evals, omega**(-1), n)
        # n_zr = ZR(n)
        # res = [x**(n_zr**(-1)) for x in fft_inv]

        # res = fnt_decode_step2(p, zs, ys_, as_, ais_, omega2, n)

        # if fft_rep[0] != coeffs:
        #     print(as_)
        #     print(ais_)
        assert fft_rep[0] == coeffs
        assert fft_rep[1] == coeffs