from adkg.polynomial import polynomials_over
from adkg.utils.misc import wrap_send, subscribe_recv
import asyncio
import time
import logging
from adkg.utils.bitmap import Bitmap
from adkg.acss_ht import ACSS_HT
from adkg.squaring import SQUARE
from adkg.all_powers import ALL_POWERS
from adkg.randousha import RANDOUSHA
from adkg.g2_powers import G2_POWERS

from adkg.broadcast.tylerba import tylerba
from adkg.broadcast.optqrbc import optqrbc

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.NOTSET)

class ADKGMsgType:
    ACSS = "A"
    RBC = "R"
    ABA = "B"
    SQ = "S"
    AP = "X"
    RDS = "D"
    G = "G"
    
class ADKG:
    def __init__(self, public_keys, private_key, g, h, g2, n, t, logq, my_id, omega2, send, recv, pc, curve_params, matrices):
        self.public_keys, self.private_key, self.g, self.h = (public_keys, private_key, g, h)
        self.g2 = g2
        self.n, self.t, self.logq, self.my_id, self.omega2 = (n, t, logq, my_id, omega2)
        self.q = 2**self.logq
        # Total number of secrets: 1 for ACS, 1 for tau, logq*(1+2) for random double sharing
        self.sc = 3*self.logq + 2 
        self.send, self.recv, self.pc = (send, recv, pc)
        self.curve_params = curve_params
        self.ZR, self.G1, self.multiexp, self.dotprod, self.blsfft = self.curve_params
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
            self.rds.kill()
            self.acss_task.cancel()
            self.rds_task.cancel()
        except Exception:
            logging.info("ADKG task finished")
        

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        return self

    async def acss_step(self, outputs, values, acss_signal):
        acsstag = ADKGMsgType.ACSS
        acsssend, acssrecv = self.get_send(acsstag), self.subscribe_recv(acsstag)
        self.acss = ACSS_HT(self.public_keys, self.private_key, self.g, self.h, self.n, self.t, self.sc, self.my_id, acsssend, acssrecv, self.pc, self.ZR, self.G1)
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

        # Invoke the protocol for randomness extraction
        rdstag = ADKGMsgType.RDS
        rdssend, rdsrecv = self.get_send(rdstag), self.subscribe_recv(rdstag)
        self.rds = RANDOUSHA(self.g, self.h, self.n, self.t, self.logq, self.my_id, rdssend, rdsrecv, (self.ZR, self.G1, self.multiexp, self.dotprod), self.matrix)
        self.rds_task = asyncio.create_task(self.rds.randousha(self.mks, acss_outputs))
        outputs  = await self.rds.output_queue.get()
        return outputs

    async def squaring(self, t_share, t_pk, t_commits, params):
        sqtag = ADKGMsgType.SQ
        sqsend, sqrecv = self.get_send(sqtag), self.subscribe_recv(sqtag)
        self.sq = SQUARE(self.g, self.n, self.t, self.logq, self.my_id, sqsend, sqrecv, (self.ZR, self.G1, self.multiexp, self.dotprod))
        self.sq_task = asyncio.create_task(self.sq.square(t_share, t_pk, t_commits, params))
        # powers-of-two shares, powers-of-two commits, and already computed powers
        pt_shares, pt_commits, powers  = await self.sq_task
        return (pt_shares, pt_commits, powers)

    
    async def g2_powers(self, t_shares, t_commits):
        gtag = ADKGMsgType.G
        gsend, grecv = self.get_send(gtag), self.subscribe_recv(gtag)
        self.gp = G2_POWERS(self.g, self.g2, self.n, self.t, self.logq, self.my_id, gsend, grecv, (self.ZR, self.G1, self.multiexp, self.dotprod))
        self.gp_task = asyncio.create_task(self.gp.powers(t_shares, t_commits))
        g2powers  = await self.gp_task
        return g2powers
    
    async def all_powers(self, t_shares, t_commits, t_powers, g2powers):
        aptag = ADKGMsgType.AP
        apsend, aprecv = self.get_send(aptag), self.subscribe_recv(aptag)
        self.ap = ALL_POWERS(self.g, self.g2, self.n, self.t, self.logq, self.my_id, self.omega2, apsend, aprecv, self.curve_params)
        self.ap_task = asyncio.create_task(self.ap.powers(t_shares, t_commits, t_powers, g2powers))
        powers  = await self.ap_task
        return powers
        
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
        setup_time = time.time() - acss_start_time
        self.benchmark_logger.info(f"Setup time: {(setup_time)}")

        mks, t_share, t_pk, pks, shares, d_shares, low_f_commits, high_f_commits = output
        params = (shares, d_shares, low_f_commits, high_f_commits)
        sq_task = asyncio.create_task(self.squaring(t_share, t_pk, pks, params))
        pt_shares, pt_commits, t_powers = await sq_task
        squaring_time = time.time() - acss_start_time
        self.benchmark_logger.info(f"Squaring time: {(squaring_time)}")

        g2_task = asyncio.create_task(self.g2_powers(pt_shares, pt_commits))
        g2powers = await g2_task
        g2_time = time.time() - acss_start_time
        self.benchmark_logger.info(f"G2 time: {(g2_time)}")

        ap_task = asyncio.create_task(self.all_powers(pt_shares, pt_commits, t_powers, g2powers))
        powers = await ap_task
        all_powers_time = time.time() - acss_start_time
        self.benchmark_logger.info(f"All powers time: {(all_powers_time)}")

        self.output_queue.put_nowait((values[1], mks, t_share, t_pk, params, pt_shares, pt_commits, powers))
        return 