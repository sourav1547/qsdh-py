import asyncio
from adkg.polynomial import polynomials_over
from adkg.utils.misc import wrap_send, subscribe_recv
from adkg.utils.serilization import Serial
from adkg.utils.poly_misc import interpolate_g1_at_x, interpolate_g1_at
from adkg.extra_proofs import CP


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
            self, g, h, n, t, logq, my_id, omega, send, recv, curve_params
    ):  # (# noqa: E501)
        self.n, self.t, self.logq, self.my_id = n, t, logq, my_id
        self.q = 2**self.logq
        self.g, self.h = g, h
        self.omega = omega
        self.ZR, self.G1, self.multiexp, self.dotprod, self.blsfft = curve_params
        self.sr = Serial(self.G1)

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
        self.output_queue = asyncio.Queue()

    def __enter__(self):
        return self

    def kill(self):
        self.subscribe_recv_task.cancel()

    def batch_mul_shares(self, cur_share, cur_powers):
        m = len(cur_powers)
        padding = 0 + (m%(self.t+1) != 0)
        ell = m//(self.t+1) + padding
        polys = [None]*ell
        for i in range(ell):
            if i == ell-1 and padding:
                polys[i] = cur_powers[i*(self.t+1):m]
            else:
                polys[i] = cur_powers[i*(self.t+1): (i+1)*(self.t+1)]
        
        cp = CP(self.g, self.g**cur_share, self.ZR, self.multiexp)
        evals = [self.blsfft(polys[i], self.omega, self.n) for i in range(ell)] 
        out_evals = [[evals[i][node]**cur_share for i in range(ell)] for node in range(self.n)]
        proofs = [[cp.prove(cur_share, evals[i][node], out_evals[node][i]) for i in range(ell)] for node in range(self.n)]
        return (out_evals, proofs)
    
    def gen_eval_shares(self, idx, all_shares):
        m = 2**idx
        ell = m//(self.t+1) + (m%(self.t+1) != 0)
        outputs = [self.G1.identity()]*ell
        for loc in range(ell):
            coords = []
            for node, shares in all_shares.items():
                coords.append([node+1, shares[loc]])
            outputs[loc] = interpolate_g1_at_x(coords, 0, self.G1, self.ZR)            
        return outputs
        

    def gen_next_powers(self, idx, all_evals):
        m = 2**idx
        padding  = 0 + (m%(self.t+1) != 0)
        ell = m//(self.t+1) + padding
        outputs = [self.G1.identity()]*m
        n_zr_inv = self.ZR(self.n)**(-1)

        for i in range(ell):
            coords = [None]*self.n
            for node, evals in all_evals.items():
                coords[node] = evals[i]

            if i == ell-1 and padding:
                fft_inv = self.blsfft(coords, self.omega**(-1), self.n)
                outputs[i*(self.t+1):i*(self.t+1)+m] = [x**n_zr_inv for x in fft_inv[:m%(self.t+1)]]
            else:
                fft_inv = self.blsfft(coords, self.omega**(-1), self.n)
                outputs[i*(self.t+1):(i+1)*(self.t+1)] = [x**n_zr_inv for x in fft_inv[:self.t+1]]
        return outputs

    def verify_share(self, sender, s_idx, s_msg):
        powers, pfs = s_msg
        return True
    
    def verify_eval(self, sender, s_idx, s_msg):
        return True

    #@profile
    async def powers(self, t_shares, t_commits, t_powers):
        """
        Generating all powers protocol
        """
        aptag = f"AP"                    
        send, recv = self.get_send(aptag), self.subscribe_recv(aptag)
        logger.debug("[%d] Starting all powers phase", self.my_id)
        self.t_shares, self.t_commits, self.t_powers = t_shares, t_commits, t_powers

        msg_buff = [{APMsgType.SHARE:{}, APMsgType.EVAL:{}} for _ in range(self.logq)]
        cur_powers = [self.g]

        
        for idx in range(self.logq):
            cur_share = self.t_shares[idx]
            share_msgs, share_pfs = self.batch_mul_shares(cur_share, cur_powers)
            for node in range(self.n):
                send(node, (APMsgType.SHARE, idx, (share_msgs[node], share_pfs[node])))
        
            temp_shares = {}
            temp_evals = {}
            while True:
                (sender, msg) = await recv()
                s_type, s_idx, s_msg = msg
                
                # Buffering future messages
                if s_idx > idx:
                    msg_buff[s_type][sender] = msg
                elif s_idx == idx:
                    # Only makes sense to process if have not processed already
                    if s_type == APMsgType.SHARE and len(temp_shares) <= self.t:    
                        if self.verify_share(sender, s_idx, s_msg):
                            powers, _ = s_msg # second coodinates are proofs
                            temp_shares[sender] = powers
                        if len(temp_shares) == self.t+1:
                            eval_msg = self.gen_eval_shares(idx, temp_shares)
                            for i in range(self.n):
                                send(i, (APMsgType.EVAL, idx, eval_msg))
                    
                    elif s_type == APMsgType.EVAL:
                        if self.verify_eval(sender, s_idx, s_msg):
                            temp_evals[sender] = s_msg
                        if len(temp_evals) == self.n:
                            next_powers = self.gen_next_powers(idx, temp_evals)
                            # Updating cur_powers as [cur_powers, next_powers]
                            cur_powers = cur_powers + next_powers
                            break

        self.output_queue.put_nowait(cur_powers)
        return