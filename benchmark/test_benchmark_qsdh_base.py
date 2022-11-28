from pytest import mark
import asyncio
from adkg.extra_proofs import PoK

import uvloop
asyncio.set_event_loop_policy(uvloop.EventLoopPolicy())
from pypairing import ZR, G1, G2, blsmultiexp as multiexp, pair
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


async def run_qsdh_base_serial(params):
    g, g2, n, logq, = params
    q = 2**logq
    output_g1 = {0:[g]*q}
    output_g2 = {0:[g2]*2}
    proofs = [None]*(n-1)
    
    def verify_prev(idx):
        cur_g = g
        for node in range(1,idx):
            pok = PoK(cur_g, ZR, multiexp)
            if not pok.verify(output_g1[node], proofs[node]):
                exit()
            cur_g = output_g1[node]
            sig1 = ZR.rand()
            sig2 = ZR.rand()
            sig1_vec = [sig1**i for i in range(q)]
            l1 = multiexp(output_g1[node], sig1_vec)
            r1 = g2*(output_g2[0]**sig2)
            l2 = g*multiexp(output_g1[node][1:], sig1_vec[1:])
            r2 = g2*(output_g2[1]**sig2)
            if not pair(l1,r1) == pair(l2,r2):
                exit()

       
    for i in range(1, n+1):
        # Verify previoius updates
        verify_prev(i)
        # Update the parameters
        alpha = ZR.rand()
        cur_alpha = alpha
        output_g1[i] = [None]*q
        output_g2[i] = [output_g2[i-1]**alpha**j for j in range(2)]
        for ii in range(q):
            output_g1[i][ii] = output_g1[i-1][ii]**(cur_alpha)
            cur_alpha = cur_alpha*alpha

        # Generate proof of knowledge of alpha
        g = output_g1[i-1][1]
        h = g**alpha
        pok = PoK(g, ZR, multiexp)
        pf = pok.prove(alpha, h)
        proofs[i] = pf
        
async def run_qsdh_base_pipe(params):
    g, g2, n, logq, = params
    q = 2**logq
    output_g1 = {0:[g]*q}
    output_g2 = {0:[g2]*2}
    proofs = [None]*(n-1)
    
    def verify(idx):
        cur_g = g
        for node in range(1,idx):
            sig_node = ZR.rand()
            sigma_vec = []