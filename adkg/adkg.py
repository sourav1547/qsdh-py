from adkg.polynomial import polynomials_over
from adkg.utils.poly_misc import interpolate_g1_at_x
from adkg.utils.misc import wrap_send, subscribe_recv
import asyncio
import hashlib, time
from math import ceil, pow
import logging
from adkg.utils.bitmap import Bitmap
from adkg.acss_ht import ACSS_HT
from adkg.squaring import SQUARE

from adkg.broadcast.tylerba import tylerba
from adkg.broadcast.optqrbc import optqrbc

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.NOTSET)

class ADKGMsgType:
    ACSS = "A"
    RBC = "R"
    ABA = "B"
    PREKEY = "P"
    KEY = "K"
    SQ = "S"
    
class CP:
    def __init__(self, g, h, ZR):
        self.g  = g
        self.h = h
        self.ZR = ZR

    def dleq_derive_chal(self, x, y, a1, a2):
        commit = str(x)+str(y)+str(a1)+str(a2)
        try:
            commit = commit.encode()
        except AttributeError:
            pass 
        hs =  hashlib.sha256(commit).digest() 
        return self.ZR.hash(hs)

    def dleq_verify(self, x, y, chal, res):
        a1 = self.multiexp([x, self.g],[chal, res])
        a2 = self.multiexp([y, self.h],[chal, res])

        eLocal = self.dleq_derive_chal(x, a1, y, a2)
        return eLocal == chal

    def dleq_prove(self, alpha, x, y):
        w = self.ZR.random()
        a1 = self.g**w
        a2 = self.h**w
        e = self.dleq_derive_chal(x, a1, y, a2)
        return  e, w - e*alpha # return (challenge, response)

class PoK:
    def __init__(self, g, ZR, multiexp):
        self.g  = g
        self.ZR = ZR
        self.multiexp = multiexp

    def pok_derive_chal(self, x, a):
        commit = str(x)+str(a)
        try:
            commit = commit.encode()
        except AttributeError:
            pass 
        hs =  hashlib.sha256(commit).digest() 
        return self.ZR.hash(hs)

    def pok_verify(self, x, chal, res):
        a = self.multiexp([x, self.g],[chal, res])
        eLocal = self.pok_derive_chal(x, a)
        return eLocal == chal

    def pok_prove(self, alpha, x):
        w = self.ZR.rand()
        a = self.g**w
        e = self.pok_derive_chal(x, a)
        return  e, w - e*alpha # return (challenge, response)
    
class ADKG:
    def __init__(self, public_keys, private_key, gs, h, n, t, logq, my_id, send, recv, pc, curve_params, matrices):
        self.public_keys, self.private_key, self.gs, self.h = (public_keys, private_key, gs, h)
        self.n, self.t, self.logq, self.my_id = (n, t, logq, my_id)
        self.q = pow(2, self.logq)
        # Total number of secrets: 1 for ACS, 1 for tau, logq*(1+2) for random double sharing
        self.sc = 3*self.logq + 2 
        self.send, self.recv, self.pc = (send, recv, pc)
        self.ZR, self.G1, self.multiexp, self.dotprod = curve_params
        self.poly = polynomials_over(self.ZR)
        self.poly.clear_cache() #FIXME: Not sure why we need this.
        # Create a mechanism to split the `recv` channels based on `tag`
        self.subscribe_recv_task, self.subscribe_recv = subscribe_recv(recv)
        self.matrix = matrices

        # Create a mechanism to split the `send` channels based on `tag`
        def _send(tag):
            return wrap_send(tag, send)
        self.get_send = _send
        self.output_queue = asyncio.Queue()


        self.benchmark_logger = logging.LoggerAdapter(
            logging.getLogger("benchmark_logger"), {"node_id": self.my_id}
        )
            
    def kill(self):
        try:
            self.subscribe_recv_task.cancel()
            for task in self.acss_tasks:
                task.cancel()
            self.acss.kill()
            self.acss_task.cancel()
        except Exception:
            logging.info("ADKG task finished")
        

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        return self

    async def acss_step(self, outputs, values, acss_signal):
        acsstag = ADKGMsgType.ACSS
        acsssend, acssrecv = self.get_send(acsstag), self.subscribe_recv(acsstag)
        self.acss = ACSS_HT(self.public_keys, self.private_key, self.gs, self.h, self.n, self.t, self.sc, self.my_id, acsssend, acssrecv, self.pc, self.ZR, self.G1)
        self.acss_tasks = [None] * self.n
        for i in range(self.n):
            if i == self.my_id:
                self.acss_tasks[i] = asyncio.create_task(self.acss.avss(0, values=values))
            else:
                self.acss_tasks[i] = asyncio.create_task(self.acss.avss(0, dealer_id=i))

        while True:
            (dealer, _, shares, commitments) = await self.acss.output_queue.get()
            outputs[dealer] = {'shares':shares, 'commits':commitments}
            if len(outputs) >= self.n - self.t:
                # print("Player " + str(self.my_id) + " Got shares from: " + str([output for output in outputs]))
                acss_signal.set()

            if len(outputs) == self.n:
                return    

    async def commonsubset(self, rbc_out, acss_outputs, acss_signal, rbc_signal, rbc_values, coin_keys, aba_in, aba_out):
        assert len(rbc_out) == self.n
        assert len(aba_in) == self.n
        assert len(aba_out) == self.n

        aba_inputted = [False]*self.n
        aba_values = [0]*self.n

        async def _recv_rbc(j):
            # rbc_values[j] = await rbc_out[j]
            rbcl = await rbc_out[j].get()
            rbcb = Bitmap(self.n, rbcl)
            rbc_values[j] = []
            for i in range(self.n):
                if rbcb.get_bit(i):
                    rbc_values[j].append(i)
                    
            if not aba_inputted[j]:
                aba_inputted[j] = True
                aba_in[j](1)
            
            subset = True
            while True:
                acss_signal.clear()
                for k in rbc_values[j]:
                    if k not in acss_outputs.keys():
                        subset = False
                if subset:
                    coin_keys[j]((acss_outputs, rbc_values[j]))
                    return
                await acss_signal.wait()

        r_threads = [asyncio.create_task(_recv_rbc(j)) for j in range(self.n)]

        async def _recv_aba(j):
            aba_values[j] = await aba_out[j]()  # May block

            if sum(aba_values) >= 1:
                # Provide 0 to all other aba
                for k in range(self.n):
                    if not aba_inputted[k]:
                        aba_inputted[k] = True
                        aba_in[k](0)
        
        await asyncio.gather(*[asyncio.create_task(_recv_aba(j)) for j in range(self.n)])
        # assert sum(aba_values) >= self.n - self.t  # Must have at least N-f committed
        assert sum(aba_values) >= 1  # Must have at least N-f committed

        # Wait for the corresponding broadcasts
        for j in range(self.n):
            if aba_values[j]:
                await r_threads[j]
                assert rbc_values[j] is not None
            else:
                r_threads[j].cancel()
                rbc_values[j] = None

        rbc_signal.set()

    async def agreement(self, key_proposal, acss_outputs, acss_signal):
        aba_inputs = [asyncio.Queue() for _ in range(self.n)]
        aba_outputs = [asyncio.Queue() for _ in range(self.n)]
        rbc_outputs = [asyncio.Queue() for _ in range(self.n)]
        
        coin_keys = [asyncio.Queue() for _ in range(self.n)]

        async def predicate(_key_proposal):
            kp = Bitmap(self.n, _key_proposal)
            kpl = []
            for ii in range(self.n):
                if kp.get_bit(ii):
                    kpl.append(ii)
            if len(kpl) <= self.t:
                return False
        
            while True:
                subset = True
                for kk in kpl:
                    if kk not in acss_outputs.keys():
                        subset = False
                if subset:
                    acss_signal.clear()    
                    return True
                acss_signal.clear()
                await acss_signal.wait()

        async def _setup(j):
            
            # starting RBC
            rbctag =ADKGMsgType.RBC + str(j) # (R, msg)
            rbcsend, rbcrecv = self.get_send(rbctag), self.subscribe_recv(rbctag)

            rbc_input = None
            if j == self.my_id: 
                riv = Bitmap(self.n)
                for k in key_proposal: 
                    riv.set_bit(k)
                rbc_input = bytes(riv.array)

            # rbc_outputs[j] = 
            asyncio.create_task(
                optqrbc(
                    rbctag,
                    self.my_id,
                    self.n,
                    self.t,
                    j,
                    predicate,
                    rbc_input,
                    rbc_outputs[j].put_nowait,
                    rbcsend,
                    rbcrecv,
                )
            )

            abatag = ADKGMsgType.ABA + str(j) # (B, msg)
            # abatag = j # (B, msg)
            abasend, abarecv =  self.get_send(abatag), self.subscribe_recv(abatag)

            def bcast(o):
                for i in range(self.n):
                    abasend(i, o)
                
            aba_task = asyncio.create_task(
                tylerba(
                    abatag,
                    self.my_id,
                    self.n,
                    self.t,
                    coin_keys[j].get,
                    aba_inputs[j].get,
                    aba_outputs[j].put_nowait,
                    bcast,
                    abarecv,
                )
            )
            return aba_task

        work_tasks = await asyncio.gather(*[_setup(j) for j in range(self.n)])
        rbc_signal = asyncio.Event()
        rbc_values = [None for i in range(self.n)]

        return (
            self.commonsubset(
                rbc_outputs,
                acss_outputs,
                acss_signal,
                rbc_signal,
                rbc_values,
                [_.put_nowait for _ in coin_keys],
                [_.put_nowait for _ in aba_inputs],
                [_.get for _ in aba_outputs],
            ),
            self.derive_key(
                acss_outputs,
                acss_signal,
                rbc_values,
                rbc_signal,
            ),
            work_tasks,
        )

    async def derive_key(self, acss_outputs, acss_signal, rbc_values, rbc_signal):
        await rbc_signal.wait()
        rbc_signal.clear()

        self.mks = set() # master key set
        for ks in  rbc_values:
            if ks is not None:
                self.mks = self.mks.union(set(list(ks)))
                if len(self.mks) >= self.n-self.t:
                    break
        
        # Waiting for all ACSS to terminate
        for k in self.mks:
            if k not in acss_outputs:
                await acss_signal.wait()
                acss_signal.clear()

        # Generating secret shares of tau
        t_secrets = [self.ZR(0)]*self.n
        t_randomness = [self.ZR(0)]*self.n
        t_commits = [self.G1.identity()]*self.n

        for node in range(self.n):
            t_secrets[node] = acss_outputs[node]['shares']['msg'][2]
            t_randomness[node] = acss_outputs[node]['shares']['rand'][1]
            t_commits[node] = acss_outputs[node]['commits'][2][0]
        
        tz_shares = [self.ZR(0)]*self.n
        tr_shares = [self.ZR(0)]*self.n
        for i in range(self.n):
            tz_shares[i] = self.dotprod(self.matrix[0][i], t_secrets)
            tr_shares[i] = self.dotprod(self.matrix[0][i], t_randomness)


        z_shares = [[self.ZR(0)]*self.logq for _ in range(self.n)]
        r_shares = [[self.ZR(0)]*self.logq for _ in range(self.n)]
        dz_shares = [[self.ZR(0)]*self.logq for _ in range(self.n)]
        dr_shares = [[self.ZR(0)]*self.logq for _ in range(self.n)]
        
        for count in range(self.logq):
            # Generating degree t shares 
            idx = count*3 + 1

            # For (n,t+1) shares of random values
            secrets = [self.ZR(0)]*self.n
            randomness = [self.ZR(0)]*self.n 
            commits = [self.G1.identity()]*self.n
            low_const, low_const_r = self.ZR(0), self.ZR(0)
            high_const, high_const_r = self.ZR(0), self.ZR(0)

            # Reordering (n,t+1) shares of random values.
            for node in range(self.n):
                if node in self.mks:
                    secrets[node] = acss_outputs[node]['shares']['msg'][idx+1]
                    randomness[node] = acss_outputs[node]['shares']['rand'][idx]
                    commits[node] = acss_outputs[node]['commits'][idx+1][0]
                    # FIXME: To update the corresponding randomness appropriately
                    low_const = low_const + secrets[node]
                    low_const_r = low_const_r + randomness[node]

            # For (n,2t+1) shares of random values
            d_secrets = [[self.ZR(0)]*self.n for _ in range(2)]
            d_randomness = [[self.ZR(0)]*self.n for _ in range(2)]
            d_commits = [[self.G1.identity()]*self.n for _ in range(2)]

            # TODO: FIXME: Currently, the secret in double sharing do no match. 
            # Without this fix, the implementation is incorrect.

            # Generating degree t shares
            for lidx in range(2):
                for node in range(self.n):
                    if node in self.mks:
                        d_secrets[lidx][node] = acss_outputs[node]['shares']['msg'][idx+lidx+1]
                        d_randomness[lidx][node] = acss_outputs[node]['shares']['rand'][idx+lidx]
                        d_commits[lidx][node] = acss_outputs[node]['commits'][idx+lidx+1][0]
                        # The constant term only depends on the first set of inputs.
                        if lidx == 0:
                            high_const = high_const + d_secrets[lidx][node]
                            high_const_r = high_const_r + d_randomness[lidx][node]
            
            for i in range(self.n):
                z_shares[i][count] = self.dotprod(self.matrix[0][i], secrets)
                r_shares[i][count] = self.dotprod(self.matrix[0][i], randomness)
                for sec in range(2):
                    dz_shares[i][count] = dz_shares[i][count] + self.dotprod(self.matrix[sec][i], d_secrets[sec])
                    dr_shares[i][count] = r_shares[i][count] + self.dotprod(self.matrix[sec][i], d_randomness[sec])
                
                # Updating the term to ensure that both low-degree and high-degree have matching constant term
                # FIXME: To update the commitments accordingly
                dz_shares[i][count] = dz_shares[i][count] + low_const - high_const
                dr_shares[i][count] = dr_shares[i][count] + low_const_r - high_const_r

            # TODO: Probably we will need to do something similar for the commitments as well.

        # Sending PREKEY messages
        keytag = ADKGMsgType.PREKEY
        send, recv = self.get_send(keytag), self.subscribe_recv(keytag)

        for i in range(self.n):
            send(i, (tz_shares[i], tr_shares[i], z_shares[i], r_shares[i], dz_shares[i], dr_shares[i]))
        
        tz_share_shares = []
        tr_share_shares = []

        # For double shares
        shares_shares = [[] for _ in range(self.logq)]
        randoms_shares = [[] for _ in range(self.logq)]
        d_shares_shares = [[] for _ in range(self.logq)]
        d_randoms_shares= [[] for _ in range(self.logq)]

        shares, randoms = [None]*self.logq, [None]*self.logq
        d_shares, d_randoms = [None]*self.logq, [None]*self.logq
        while True:
            (sender, msg) = await recv()
            t_share_p, t_random_p, shares_p, randoms_p, d_shares_p, d_randoms_p = msg

            tz_share_shares.append([sender+1, t_share_p])
            tr_share_shares.append([sender+1, t_random_p])

            for ii in range(self.logq):
                shares_shares[ii].append([sender+1, shares_p[ii]])
                randoms_shares[ii].append([sender+1, randoms_p[ii]])
                d_shares_shares[ii].append([sender+1, d_shares_p[ii]])
                d_randoms_shares[ii].append([sender+1, d_randoms_p[ii]])


            # Interpolating the share
            if len(tz_share_shares) >= self.t+1:    
                t_share =  self.poly.interpolate_at(tz_share_shares, 0)
                t_random =  self.poly.interpolate_at(tr_share_shares, 0)

                for ii in range(self.logq):
                    shares[ii] = self.poly.interpolate_at(shares_shares[ii], 0)
                    randoms[ii] = self.poly.interpolate_at(randoms_shares[ii], 0)
                    d_shares[ii] = self.poly.interpolate_at(d_shares_shares[ii], 0)
                    d_randoms[ii] = self.poly.interpolate_at(d_randoms_shares[ii], 0)
                
                break

                # TODO(@sourav): Implement the verification check and also the fallback path
                # commit = self.G1.identity()
                # for sec in range(self.sc-1):
                #     commit = commit*self.multiexp(commits[sec], self.matrix[sec][self.my_id])
                # if self.multiexp([self.g, self.h],[secret, random]) == commit:
                #     break
                

        mt = self.gs[0]**t_share
        mtr = self.h**t_random
        gpok = PoK(self.gs[0], self.ZR, self.multiexp)
        hpok = PoK(self.h, self.ZR, self.multiexp)
        gchal, gres = gpok.pok_prove(t_share, mt)
        hchal, hres = hpok.pok_prove(t_random, mtr)

        # TODO: We will need to generate proof of knowledge for these values as well.
        low_commits = [self.gs[0]**shares[ii] for ii in range(self.logq)]
        low_r_commits = [self.gs[0]**randoms[ii] for ii in range(self.logq)]
        high_commits = [self.gs[0]**d_shares[ii] for ii in range(self.logq)]
        high_r_commits = [self.gs[0]**d_randoms[ii] for ii in range(self.logq)]

        keytag = ADKGMsgType.KEY
        send, recv = self.get_send(keytag), self.subscribe_recv(keytag)

        for i in range(self.n):
            send(i, (mt, mtr, gchal, gres, hchal, hres, low_commits, low_r_commits, high_commits, high_r_commits))
        
        pk_shares = [[self.my_id+1, mt]]
        rk_shares = [[self.my_id+1, mtr]]

        low_f_commits = [[] for _ in range(self.logq)]
        high_f_commits = [[] for _ in range(self.logq)]

        for ii in range(self.logq):
            low_f_commits[ii].append([self.my_id+1, low_commits[ii]])
            high_f_commits[ii].append([self.my_id+1, high_commits[ii]])


        while True:
            (sender, msg) = await recv()
            if sender != self.my_id:
                x, y, gchal, gres, hchal, hres, low_commits_p, _, high_commits_p, _ = msg
                valid_pok = gpok.pok_verify(x, gchal, gres) and hpok.pok_verify(y, hchal, hres)
                if valid_pok:
                    pk_shares.append([sender+1, x])
                    rk_shares.append([sender+1, y])

                    for ii in range(self.logq):
                        low_f_commits[ii].append([sender+1, low_commits_p[ii]])
                        high_f_commits[ii].append([sender+1, high_commits_p[ii]])

            if len(pk_shares) > 2*self.t:
                break
        pk =  interpolate_g1_at_x(pk_shares, 0, self.G1, self.ZR)
        rk =  interpolate_g1_at_x(rk_shares, 0, self.G1, self.ZR)
        com0 = self.multiexp(t_commits, [self.ZR(1)]*self.n)
        # TODO:(@sourav) FIXME! To do FFT in the exponent here
        # TODO:(@sourav) FIXME! Add the fallback path
        assert pk*rk == com0
        return (self.mks, t_share, pk, pk_shares, shares, d_shares, low_f_commits, high_f_commits)

    async def squaring(self, t_share, t_pk, t_commits, params):
        sqtag = ADKGMsgType.SQ
        sqsend, sqrecv = self.get_send(sqtag), self.subscribe_recv(sqtag)
        self.sq = SQUARE(self.gs[0], self.h, self.n, self.t, self.logq, self.my_id, sqsend, sqrecv, (self.ZR, self.G1, self.multiexp, self.dotprod))
        self.sq_task = asyncio.create_task(self.sq.square(t_share, t_pk, t_commits, params))
        # powers-of-two shares, powers-of-two commits, and already computed powers
        pt_shares, pt_commits, powers  = await self.sq.output_queue.get()
        return (pt_shares, pt_commits, powers)
        
    async def run_adkg(self, start_time):
        logging.info(f"Run ADKG called")
        acss_outputs = {}
        acss_signal = asyncio.Event()

        acss_start_time = time.time()
        values =[self.ZR.rand() for _ in range(self.sc)]
        self.acss_task = asyncio.create_task(self.acss_step(acss_outputs, values, acss_signal))
        await acss_signal.wait()
        acss_signal.clear()
        acss_time = time.time() - acss_start_time
        self.benchmark_logger.info(f"ACSS time: {(acss_time)}")
        key_proposal = list(acss_outputs.keys())
        create_acs_task = asyncio.create_task(self.agreement(key_proposal, acss_outputs, acss_signal))
        acs, key_task, work_tasks = await create_acs_task
        await acs
        output = await key_task
        await asyncio.gather(*work_tasks)
        mks, t_share, t_pk, pks, shares, d_shares, low_f_commits, high_f_commits = output
        params = (shares, d_shares, low_f_commits, high_f_commits)
        sq_task = asyncio.create_task(self.squaring(t_share, t_pk, pks, params))
        pt_shares, pt_commits, powers = await sq_task
        # create_qsdh_task = asyncio.create_task(self.qsdh(output))
        # params = await create_qsdh_task
        # qsdh_time = time.time()-start_time
        # logging.info("ADKG time 2: %f", qsdh_time)
        # self.benchmark_logger.info("ADKG time: %f", qsdh_time)
        self.output_queue.put_nowait((values[1], mks, t_share, t_pk, params, pt_shares, pt_commits, powers))