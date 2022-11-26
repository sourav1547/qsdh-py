from adkg.poly_commit_hybrid import PolyCommitHybrid
from pytest import mark
from adkg.polynomial import polynomials_over
from adkg.utils.poly_misc import get_omega
from adkg.adkg import ADKG
import asyncio
import numpy as np
import uvloop
asyncio.set_event_loop_policy(uvloop.EventLoopPolicy())
from pypairing import ZR, G1, G2, blsmultiexp as multiexp, dotprod, blsfft

def get_avss_params(n, G1, G2):
    g, h = G1.rand(b'g'), G1.rand(b'h')
    g2 = G2.rand(b'g')
    public_keys, private_keys = [None] * n, [None] * n
    for i in range(n):
        private_keys[i] = ZR.hash(str(i).encode())
        public_keys[i] = g**private_keys[i]
    return g, h, g2, public_keys, private_keys

@mark.parametrize("t, logq", [(5, 5), (5,10), (5,15)])
def test_benchmark_qsdh_base(test_router, benchmark, t, logq):
    loop = asyncio.get_event_loop()
    n = 3*t + 1
    g, _, g2, _, _ = get_avss_params(n, G1, G2)
    params = (g, g2, n, logq)

    def _prog():
        loop.run_until_complete(run_qsdh_base(params))

    benchmark(_prog)

async def run_qsdh_base(params):
    g, g2, n, logq, = params
    q = 2**logq
    output_g1 = {0:[g]*q}
    # output_g2 = {0:[g2]*q}
    for i in range(1, n+1):
        alpha = ZR.rand()
        cur_alpha = alpha
        output_g1[i] = [None]*q
        for ii in range(q):
            output_g1[i][ii] = output_g1[i-1][ii]**(cur_alpha)
            cur_alpha = cur_alpha*alpha
            # output_g2[i][ii] = output_g2[i-1][ii]**(alpha**ii)


async def run_qsdh_base_ver(params):
    g, g2, n, logq, = params
    q = 2**logq
    output_g1 = {0:[g]*q}
    # output_g2 = {0:[g2]*q}
    for i in range(1, n+1):
        alpha = ZR.rand()
        cur_alpha = alpha
        output_g1[i] = [None]*q
        for ii in range(q):
            output_g1[i][ii] = output_g1[i-1][ii]**(cur_alpha)
            cur_alpha = cur_alpha*alpha
            # output_g2[i][ii] = output_g2[i-1][ii]**(alpha**ii)