from pypairing import ZR, G1, blsfft, blsmultiexp, robustblsfft
from adkg.polynomial import polynomials_over
from adkg.utils.poly_misc import prep_for_fft, prep_for_fft_batch
from random import shuffle
from pytest import mark
import asyncio

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


@mark.parametrize("n, t, logq", [(16, 10, 5)])
@mark.benchmark(
    min_rounds=1,
    max_time=0.0005,
)
# @mark.parametrize("n, t, logq", [(64, 42, 10)])
def test_benchmark_qsdh_base(test_router, benchmark, n, t, logq):
    loop = asyncio.get_event_loop()
    omega2 = get_omega(ZR, 2*n)
    omega = omega2**2
    g = G1.rand(b'')

    p = polynomials_over(ZR)
    q = 2**logq
    m = q//(t+1)
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

    def _prog():
        loop.run_until_complete(run_base(zs, ys, omega2, n))

    benchmark(_prog)

async def run_base(zs, ys, omega2, n):
    robustblsfft(zs, ys, omega2, n)

# @mark.parametrize("n, t, logq", [(16, 10, 10), (16, 10, 15), (32, 21, 10), (32, 21, 15), (64, 42, 10), (64, 42, 15), (128, 85, 10), (32, 85, 15)])
@mark.parametrize("n, t, logq", [(16, 10, 10), (32, 21, 10)])
# @mark.parametrize("n, t, logq", [(64, 42, 10)])
def test_benchmark_qsdh_naive(test_router, benchmark, n, t, logq):
    loop = asyncio.get_event_loop()
    omega = get_omega(ZR, n)
    g = G1.rand(b'')

    p = polynomials_over(ZR)
    q = 2**logq
    m = q//(t+1)

    zs = list(range(n))
    shuffle(zs)
    zs = zs[:t+1]
    polys = [p.random(t) for _ in range(m)]
    all_ys = [polys[i].evaluate_fft(omega, n) for i in range(m)]
    all_evals = {}
    for node in zs:
        all_evals[node] = []
    
    for i in range(m):
        for node in zs:
            all_evals[node].append(g**all_ys[i][node])

    def _prog():
        loop.run_until_complete(run_basic_interpolate(all_evals, m, omega, t, n))
    benchmark(_prog)

async def run_basic_interpolate(all_evals, ell, omega, t, n):
    omegainv = omega**(-1)
    ninv = ZR(n)**(-1)

    for i in range(ell):
        coords = [None]*n
        for node, evals in all_evals.items():
            coords[node] = evals[i]
        coords = prep_for_fft(coords, omega, n, blsmultiexp, ZR)
        
        fft_inv = blsfft(coords, omegainv, n)
        [x**ninv for x in fft_inv[:t+1]]


# @mark.parametrize("n, t, logq", [(64, 42, 10)])
@mark.parametrize("n, t, logq", [(16, 10, 10), (32, 21, 10)])
def test_benchmark_qsdh_naive_batch(test_router, benchmark, n, t, logq):
    loop = asyncio.get_event_loop()
    omega = get_omega(ZR, n)
    g = G1.rand(b'')

    p = polynomials_over(ZR)
    q = 2**logq
    m = q//(t+1)

    zs = list(range(n))
    shuffle(zs)
    zs = zs[:t+1]
    polys = [p.random(t) for _ in range(m)]
    all_ys = [polys[i].evaluate_fft(omega, n) for i in range(m)]
    all_evals = {}
    for node in zs:
        all_evals[node] = []
    
    for i in range(m):
        for node in zs:
            all_evals[node].append(g**all_ys[i][node])

    def _prog():
        loop.run_until_complete(run_batch_interpolate(all_evals, m, omega, t, n))
    benchmark(_prog)

async def run_batch_interpolate(all_evals, ell, omega, t, n):
    omegainv = omega**(-1)
    ninv = ZR(n)**(-1)

    xs = all_evals.keys()
    ys = list(all_evals.values())
    coords = prep_for_fft_batch(xs, ys, omega, ell, n, blsmultiexp, ZR)

    for i in range(ell):
        fft_inv = blsfft(coords[i], omegainv, n)
        [x**ninv for x in fft_inv[:t+1]]