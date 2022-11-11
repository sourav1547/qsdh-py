from adkg.poly_commit_hybrid import PolyCommitHybrid
from pytest import mark
from adkg.polynomial import polynomials_over
from adkg.adkg import ADKG
import asyncio
import numpy as np
import math
import uvloop
asyncio.set_event_loop_policy(uvloop.EventLoopPolicy())
from pypairing import ZR, G1, blsmultiexp as multiexp, dotprod
# from pypairing import Curve25519ZR as ZR, Curve25519G as G1, curve25519multiexp as multiexp, curve25519dotprod as dotprod
    
import time

def get_avss_params(n, logq, G1):
    # g, h = G1.rand(b'g'), G1.rand(b'h')
    gs = [G1.rand(str(i).encode()) for i in range(logq)]
    h = G1.rand(b'h')
    public_keys, private_keys = [None] * n, [None] * n
    for i in range(n):
        private_keys[i] = ZR.hash(str(i).encode())
        public_keys[i] = gs[0]**private_keys[i]
    return gs, h, public_keys, private_keys


def gen_vector(t, n):
    deg = 2*t
    coeff_1 = np.array([[ZR(i+1)**j for j in range(t+1)] for i in range(n)])
    coeff_2 = np.array([[ZR(i+1)**j for j in range(t+1, deg+1)] for i in range(n)])
    hm_1 = np.array([[ZR(i+1)**j for j in range(n)] for i in range(t+1)])
    hm_2 = np.array([[ZR(i+1)**j for j in range(n)] for i in range(deg-t)])
    rm_1 = np.matmul(coeff_1, hm_1)
    rm_2 = np.matmul(coeff_2, hm_2)

    return (rm_1.tolist(), rm_2.tolist())

@mark.asyncio
async def test_adkg(test_router):
    t = 2
    logq = 5
    q = math.pow(2,logq)
    n = 3 * t + 1
    gs, h, pks, sks = get_avss_params(n, logq, G1)
    sends, recvs, _ = test_router(n, maxdelay=0.01)
    pc = PolyCommitHybrid(gs, h, ZR, multiexp)
    mat1, mat2 = gen_vector(t, n)

    dkg_tasks = [None] * n # async task for adkg
    dkg_list = [None] * n #

    start_time = time.time()
    curve_params = (ZR, G1, multiexp, dotprod)

    for i in range(n):
        dkg = ADKG(pks, sks[i], gs, h, n, t, logq, i, sends[i], recvs[i], pc, curve_params, (mat1, mat2))
        dkg_list[i] = dkg
        dkg_tasks[i] = asyncio.create_task(dkg.run_adkg(start_time))
    
    outputs = await asyncio.gather(
        *[dkg_list[i].output_queue.get() for i in range(n)]
    )
    for dkg in dkg_list:
        dkg.kill()
    for task in dkg_tasks:
        task.cancel()
    
    
    shares = []
    low_doubles = [[] for _ in range(logq)]
    high_doubles = [[] for _ in range(logq)]

    i = 1
    for _, _, sk, _, doubles in outputs:
        shares.append([i, sk])
        l_shares, h_shares, _, _ = doubles
        for ii in range(logq):
            low_doubles[ii].append([i, l_shares[ii]])
            high_doubles[ii].append([i, h_shares[ii]])
        i = i + 1

    poly = polynomials_over(ZR)
    for ii in range(logq):
        l_const = poly.interpolate_at(low_doubles[ii],0)
        h_const = poly.interpolate_at(high_doubles[ii],0)
        assert l_const == h_const
    
    # msk = poly.interpolate_at(shares,0)
    # mpk = gs[0]**msk

    # for i in range(n):
    #     assert(mpk == outputs[i][3])

    # mks_set = outputs[0][1]
    # for i in range(1, n):
    #     assert mks_set == outputs[i][1]

    # mks_sum = ZR(0)
    # for node in mks_set:
    #     mks_sum = mks_sum + outputs[node][0]
    # assert msk == mks_sum

    def check_degree(claimed_degree, points):
        dual_code = gen_dual_code(n, claimed_degree, poly)
        check = dot(points, dual_code)
        return check == ZR(0)

    def gen_dual_code(n, degree, poly):
        def get_vi(i, n):
            out = ZR(1)
            for j in range(1, n+1):
                if j != i:
                    out = out / (i-j)
            return out
        q = poly.random(n -degree -2)
        q_evals = [q(i+1) for i in range(n)]
        return [q_evals[i] * get_vi(i+1, n) for i in range(n)]
    

    def dot(a, b):
        res = ZR(0)
        for i in range(len(a)):
            res = res + a[i][1]*b[i]
        return res
    
    for ii in range(logq):
        assert not check_degree(2*t-1, high_doubles[ii])
        assert check_degree(2*t, high_doubles[ii])

        assert not check_degree(t-1, low_doubles[ii])
        assert check_degree(t, low_doubles[ii])

    # assert not check_degree(deg-1, shares)
    # assert check_degree(deg, shares)
    assert True