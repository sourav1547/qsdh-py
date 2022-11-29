import asyncio
from adkg.extra_proofs import PoK
import uvloop
from pypairing import ZR, G1, G2, blsmultiexp as multiexp, pair


asyncio.set_event_loop_policy(uvloop.EventLoopPolicy())

def test_qsdh_base_no_verf():
    t = 1
    n = 3*t + 1
    logq = 5
    g = G1.rand(b'g')
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

def verify_prev(idx, output_g1, output_g2, q, g1, g2, proofs):
    cur_g = g1
    for node in range(1,idx):
        pok = PoK(cur_g, ZR, multiexp)
        if not pok.verify(output_g1[node][0], proofs[node]):
            exit()
        cur_g = output_g1[node]
        sig1 = ZR.rand()
        sig2 = ZR.rand()
        sig1_vec = [sig1**i for i in range(q)]
        l1 = multiexp(output_g1[node], sig1_vec)
        r1 = g2*(output_g2[node][0]**sig2)
        l2 = g1*multiexp(output_g1[node][1:], sig1_vec[1:])
        r2 = g2*(output_g2[node][1]**sig2)
        if not pair(l1,r1) == pair(l2,r2):
            print("Pairing test failed")


def test_qsdh_base_serial():
    t = 1
    n = 3*t + 1
    logq = 5 
    g1, g2 = G1.rand(b'g'), G2.rand(b'g')
    q = 2**logq

    output_g1 = {0:[g1]*q}
    output_g2 = {0:[g2]*2}
    proofs = [None]*(n-1)
    for i in range(1, n+1):
        # Verify previoius updates
        verify_prev(i, output_g1, output_g2, q, g1, g2, proofs)
        # Update the parameters
        alpha = ZR.rand()
        cur_alpha = alpha
        output_g1[i] = [None]*q
        output_g2[i] = [output_g2[i-1][j]**(alpha**j) for j in range(2)]
        for ii in range(q):
            output_g1[i][ii] = output_g1[i-1][ii]**(cur_alpha)
            cur_alpha = cur_alpha*alpha

        # Generate proof of knowledge of alpha
        g = output_g1[i-1][0]
        h = g**alpha
        pok = PoK(g, ZR, multiexp)
        pf = pok.prove(alpha, h)
        proofs[i] = pf
        
def test_qsdh_base_pipe():
    t = 1
    n = 3*t + 1
    logq = 5 
    g, g2 = G1.rand(b'g'), G2.rand(b'g')
    q = 2**logq


    output_g1 = {0:[g]*q}
    output_g2 = {0:[g2]*2}
    proofs = [None]*(n-1)