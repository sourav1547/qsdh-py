from adkg.polynomial import polynomials_over
from adkg.utils.misc import wrap_send, subscribe_recv
from adkg.utils.serilization import Serial
from adkg.utils.poly_misc import prep_for_fft, interpolate_g1_batch_at, prep_for_fft_batch
from adkg.extra_proofs import CP
from pypairing import pair, robustblsfft
import time
import hashlib


import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

# Uncomment this when you want logs from this file.
# logger.setLevel(logging.DEBUG)
class APMsgType:
    SHARE = "S"
    EVAL = "E"

class ALL_POWERS:
    #@profile
    def __init__(
            self, g, g2, n, t, logq, my_id, omega2, send, recv, curve_params
    ):  # (# noqa: E501)
        self.n, self.t, self.logq, self.my_id = n, t, logq, my_id
        self.g, self.g2 = g, g2
        self.omega2 =  omega2
        self.omega = omega2**2
        self.omegainv =  self.omega**(-1)
        self.ZR, self.G1, self.multiexp, self.dotprod, self.blsfft = curve_params
        self.sr = Serial(self.G1)
        self.ninv = self.ZR(n)**(-1)
        self.rands = [self.ZR.rand() for _ in range(2**self.logq)]

        self.benchmark_logger = logging.LoggerAdapter(
            logging.getLogger("benchmark_logger"), {"node_id": self.my_id}
        )

        # Create a mechanism to split the `recv` channels based on `tag`
        self.subscribe_recv_task, self.subscribe_recv = subscribe_recv(recv)

        # Create a mechanism to split the `send` channels based on `tag`
        def _send(tag):
            return wrap_send(tag, send)

        self.get_send = _send
        self.poly = polynomials_over(self.ZR)
        self.poly.clear_cache()

    def __enter__(self):
        return self

    def kill(self):
        self.subscribe_recv_task.cancel()

    def batch_mul_shares(self, idx, cur_share):
        start = time.time()
        # m = len(cur_powers)
        m = 2**idx
        padding = 0 + (m%(self.deg) != 0)
        ell = m//(self.deg) + padding
        polys = [None]*ell
        for i in range(ell):
            if i == ell-1 and padding:
                polys[i] = self.cur_powers[i*(self.deg):m]
            else:
                polys[i] = self.cur_powers[i*(self.deg): (i+1)*(self.deg)]
            
        self.evals = [self.blsfft(polys[i], self.omega, self.n) for i in range(ell)] 
        out_evals = [[self.evals[i][node]**cur_share for i in range(ell)] for node in range(self.n)]
        
        # computing the seed for random linear combination
        # Only compute once while sending the messages
        datas = str(self.evals)
        try:
            commit = datas.encode()
        except AttributeError:
            pass 
        hs =  hashlib.sha256(commit).digest() 
        seed = self.ZR.hash(hs)

        self.rands2[idx] = [self.ZR.hash(str(seed+i).encode()) for i in range(ell)]
        proofs = [None]*self.n
        x = self.g**cur_share
        self.h_vector[idx] = [self.G1.identity()]*self.n
        for node in range(self.n):
            h = self.multiexp([self.evals[i][node] for i in range(ell)], self.rands2[idx])
            if node == self.my_id:
                self.h_vector[idx] = h
            y = self.multiexp(out_evals[node], self.rands2[idx])
            cp = CP(self.g, h, self.ZR, self.multiexp)
            proofs[node] = cp.prove(cur_share, x, y)
        time_taken = time.time() - start
        self.benchmark_logger.info(f"SAHRE gen, ell:{ell}, time:{time_taken}")
        return (out_evals, proofs)

    def gen_eval_shares(self, idx, all_shares):
        m = 2**idx
        ell = m//(self.deg) + (m%(self.deg) != 0)
        xs = []
        ys_shares = [[] for _ in range(ell)]
        for node, shares in all_shares.items():
            xs.append(node+1)
            for loc in range(ell):
                ys_shares[loc].append(shares[loc])
        
        return interpolate_g1_batch_at(xs, ys_shares, 0, self.multiexp, self.ZR)    

    def update_powers(self, idx, all_evals):
        start = time.time()
        m = 2**idx
        padding  = 0 + (m%(self.deg) != 0)
        ell = m//(self.deg) + padding
        
        zs = list(all_evals.keys())
        assert len(zs) == self.deg
        ys = [None]*(self.deg*ell)

        for i in range(ell):
            ii = 0
            for _, evals in all_evals.items():
                ys[i*self.deg+ii] = evals[i]
                ii = ii + 1
        fft_rep = robustblsfft(zs, ys, self.omega2, self.n)
        for i in range(ell):    
            if i == ell-1 and padding:
                self.cur_powers[m+i*(self.deg):m+i*(self.deg)+m] = fft_rep[i][:m%self.deg]
            else:
                self.cur_powers[m+i*(self.deg):m+(i+1)*(self.deg)] = fft_rep[i]

        time_taken = time.time() - start
        self.benchmark_logger.info(f"Next message gen, ell:{ell}, time:{time_taken}")

    def update_powers_naive(self, idx, all_evals):
        start = time.time()
        m = 2**idx
        padding  = 0 + (m%(self.deg) != 0)
        ell = m//(self.deg) + padding
        outputs = [None]*m

        xs = all_evals.keys()
        ys = list(all_evals.values())
        coords = prep_for_fft_batch(xs, ys, self.omega, ell, self.n, self.multiexp, self.ZR)
        for i in range(ell):
            if i == ell-1 and padding:
                fft_inv = self.blsfft(coords[i], self.omegainv, self.n)
                outputs[i*(self.deg):i*(self.deg)+m] = [x**self.ninv for x in fft_inv[:m%(self.deg)]]
            else:
                fft_inv = self.blsfft(coords[i], self.omegainv, self.n)
                outputs[i*(self.deg):(i+1)*(self.deg)] = [x**self.ninv for x in fft_inv[:self.deg]]

        time_taken = time.time() - start
        self.benchmark_logger.info(f"Next message gen, ell:{ell}, time:{time_taken}")
        self.cur_powers[2**idx:2**(idx+1)] = outputs
        return


    def verify_share(self, sender, s_idx, powers, pf):
        x = self.t_commits[s_idx][sender+1]
        y = self.multiexp(powers, self.rands2[s_idx])
        h = self.h_vector[s_idx]
        cp = CP(self.g, h, self.ZR, self.multiexp)
        return cp.verify(x, y, pf)
    
    def verify_eval(self, sender, idx, powers):
        ell = len(powers)
        g1rlc = self.multiexp(powers, self.rands[:ell])
        l_g1rlc = self.multiexp([self.evals[i][sender] for i in range(ell)], self.rands[:ell])
        return pair(g1rlc, self.g2) == pair(l_g1rlc, self.g2powers[idx])

    def verify_msgs(self, type, s_idx):
        v_msgs = {}
        msgs = self.msg_buff[s_idx][type]
        if type == APMsgType.SHARE:
            for sender, s_msg in msgs.items():
                powers, pf = self.decode(s_msg)
                if self.verify_share(sender, s_idx, powers, pf):
                    v_msgs[sender] = powers
                if len(v_msgs) == self.t+1:
                    break
        elif type == APMsgType.EVAL:
            for sender, s_msg in msgs.items():
                s_powers = self.decode(s_msg)
                if self.verify_eval(sender, s_idx, s_powers):
                    v_msgs[sender] = s_powers
                if len(v_msgs) == self.deg:
                    break
        return v_msgs

    def encode(self, msgs, pf):
        chal, res = pf
        datab = self.sr.serialize_gs(msgs)
        datab.extend(self.sr.serialize_f(chal))
        datab.extend(self.sr.serialize_f(res))
        return datab
    
    def decode(self, data):
        g_size = 48
        f_size = 32
        data = bytes(data)
        ell = (len(data)-2*f_size)//g_size
        shares = self.sr.deserialize_gs(data[:ell*g_size])
        chal = self.sr.deserialize_f(data[ell*g_size: ell*g_size+f_size])
        resp = self.sr.deserialize_f(data[ell*g_size+f_size:])
        return shares, (chal, resp)

    #@profile
    async def powers(self, t_shares, t_commits, t_powers, g2powers):
        """
        Generating all powers protocol
        """
        aptag = f"AP"        
        self.start_time  = time.time()          
        send, recv = self.get_send(aptag), self.subscribe_recv(aptag)
        logger.debug("[%d] Starting all powers phase", self.my_id)
        self.t_shares, self.t_commits, self.t_powers, self.g2powers = t_shares, t_commits, t_powers, g2powers
        self.deg = self.n-self.t

        self.msg_buff = [{APMsgType.SHARE:{}, APMsgType.EVAL:{}} for _ in range(self.logq)]
        self.cur_powers = [self.G1.identity()]*(2**self.logq)
        self.cur_powers[0] = self.g
        self.rands2 = {}
        self.h_vector = {}

        for idx in range(self.logq):
            elapsed_time = time.time() - self.start_time
            self.benchmark_logger.info(f"AP start index:{idx}, elapsed time {(elapsed_time)}")
            cur_share = self.t_shares[idx]
            share_msgs, share_pfs = self.batch_mul_shares(idx, cur_share)
            for node in range(self.n):
                datab = self.encode(share_msgs[node], share_pfs[node])
                send(node, (APMsgType.SHARE, idx, datab))

            # Processing old buffered messages
            temp_shares = self.verify_msgs(APMsgType.SHARE, idx)
            if len(temp_shares) > self.t:
                eval_powers = self.gen_eval_shares(idx, temp_shares)
                eval_msg = self.sr.serialize_gs(eval_powers)
                for i in range(self.n):
                    send(i, (APMsgType.EVAL, idx, eval_msg))
            temp_evals = self.verify_msgs(APMsgType.EVAL, idx)
            if len(temp_evals) == self.deg:
                self.update_powers(idx, temp_evals)
                continue
            
            while True:
                sender, msg = await recv()
                s_type, s_idx, s_msg = msg
                
                # Buffering future messages
                if s_idx > idx:
                    self.msg_buff[s_idx][s_type][sender] = msg
                elif s_idx == idx:
                    # Only makes sense to process if have not processed already
                    if s_type == APMsgType.SHARE and len(temp_shares) <= self.t:    
                        powers, pfs = self.decode(s_msg)
                        if self.verify_share(sender, s_idx, powers, pfs):
                            temp_shares[sender] = powers
                        if len(temp_shares) == self.t+1:
                            eval_powers = self.gen_eval_shares(idx, temp_shares)
                            eval_msg = self.sr.serialize_gs(eval_powers)
                            for i in range(self.n):
                                send(i, (APMsgType.EVAL, idx, eval_msg))
                    
                    elif s_type == APMsgType.EVAL:
                        s_powers = self.sr.deserialize_gs(bytes(s_msg))
                        if self.verify_eval(sender, s_idx, s_powers):
                            temp_evals[sender] = s_powers
                        if len(temp_evals) == self.deg:
                            # Updating cur_powers with newly computed values
                            self.update_powers(idx, temp_evals)
                            break

        return self.cur_powers