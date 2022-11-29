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

def verify_prev(start, cnode, output_g1, output_g2, q, g1, g2, proofs):
    cur_g = output_g1[start-1][0]
    for pnode in range(start, cnode):
        pok = PoK(cur_g, ZR, multiexp)
        assert pok.verify(output_g1[pnode][0], proofs[pnode])
        cur_g = output_g1[pnode][0]
        sig1 = ZR.rand()
        sig2 = ZR.rand()
        sig1_vec = [sig1**i for i in range(q)]
        l1 = multiexp(output_g1[pnode], sig1_vec)
        r1 = g2*(output_g2[pnode][0]**sig2)
        l2 = g1*multiexp(output_g1[pnode][:q-1], sig1_vec[1:])
        r2 = output_g2[pnode][1]*(output_g2[pnode][1]**sig2)
        assert pair(l1,r1) == pair(l2,r2)


def test_qsdh_base_serial():
    t = 1
    n = 3*t + 1
    logq = 5 
    g1, g2 = G1.rand(b'g'), G2.rand(b'g')
    q = 2**logq

    output_g1 = {0:[g1]*q}
    output_g2 = {0:[g2]*2}
    proofs = [None]*(n+1)
    for node in range(1, n+1):
        # Verify previoius updates
        verify_prev(1, node, output_g1, output_g2, q, g1, g2, proofs)
        # Update the parameters
        alpha = ZR.rand()
        cur_alpha = alpha
        output_g1[node] = [None]*q
        output_g2[node] = [output_g2[node-1][j]**(alpha**j) for j in range(2)]
        for ii in range(q):
            output_g1[node][ii] = output_g1[node-1][ii]**(cur_alpha)
            cur_alpha = cur_alpha*alpha

        # Generate proof of knowledge of alpha
        g = output_g1[node-1][0]
        h = g**alpha
        pok = PoK(g, ZR, multiexp)
        proofs[node] = pok.prove(alpha, h)
        
def test_qsdh_base_pipe_verf():
    t = 1
    n = 3*t + 1
    logq = 5 
    g1, g2 = G1.rand(b'g'), G2.rand(b'g')
    q = 2**logq

    output_g1 = {0:[g1]*q}
    output_g2 = {0:[g2]*2}
    proofs = [None]*(n+1)
    for node in range(1, n+1):
        # Verify previoius updates
        start = node-1
        if node == 1:
            start = 1
        verify_prev(start, node, output_g1, output_g2, q, g1, g2, proofs)
        # Update the parameters
        alpha = ZR.rand()
        cur_alpha = alpha
        output_g1[node] = [None]*q
        output_g2[node] = [output_g2[node-1][j]**(alpha**j) for j in range(2)]
        for ii in range(q):
            output_g1[node][ii] = output_g1[node-1][ii]**(cur_alpha)
            cur_alpha = cur_alpha*alpha

        # Generate proof of knowledge of alpha
        g = output_g1[node-1][0]
        h = g**alpha
        pok = PoK(g, ZR, multiexp)
        proofs[node] = pok.prove(alpha, h)