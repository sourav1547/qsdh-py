from pypairing import ZR, G1, blsfft, blsmultiexp, robustblsfft
from adkg.polynomial import polynomials_over
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
        ys = [g**ys[i] for i in zs]

        coeffs = [g**x for x in poly.coeffs]
        fft_rep = robustblsfft(zs, [ys], omega2, n)

        assert fft_rep[0] == coeffs