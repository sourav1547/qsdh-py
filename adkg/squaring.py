import asyncio
from adkg.polynomial import polynomials_over
from adkg.utils.misc import wrap_send, subscribe_recv
from adkg.utils.serilization import Serial
from adkg.extra_proofs import CP


import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

# Uncomment this when you want logs from this file.
# logger.setLevel(logging.DEBUG)

class SQUARE:
    #@profile
    def __init__(
            self, g, n, t, logq, my_id, send, recv, curve_params
    ):  # (# noqa: E501)
        self.n, self.t, self.logq, self.my_id = n, t, logq, my_id
        self.g = g
        self.ZR, self.G1, self.multiexp, _ = curve_params
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

    def __enter__(self):
        return self

    def kill(self):
        self.subscribe_recv_task.cancel()
        
    def verify_sq(self, idx, sender, s_reveal, s_y, s_pf):
        s_x = self.th_powers[idx-1][sender+1]
        cp = CP(self.g, s_x, self.ZR, self.multiexp)
        if cp.verify(s_x, s_y, s_pf):
            return self.multiexp([self.high_f_commits[idx-1][sender+1], self.g], [self.ZR(1), s_reveal]) == s_y
        return False

    def process(self, idx, sq_shares):
        # interpolating revealed message and next share
        reveal = self.poly.interpolate_at(sq_shares, 0)
        cur_share = reveal + self.shares[idx-1]
        self.out_shares[idx] = cur_share
        
        # Computing tau^{2^idx}.g and threshold public keys
        self.powers[idx] = (self.g**reveal)*self.high_f_commits[idx-1][0]
        self.th_powers[idx] = [(self.g**reveal)*self.low_f_commits[idx-1][i] for i in range(self.n+1)]
        self.inc_idx = True
        return cur_share

    
    #@profile
    async def square(self, t_share, t_pk, tpks, params):
        """
        Squaring protocol
        """
        sqtag = f"SQ"                    
        send, recv = self.get_send(sqtag), self.subscribe_recv(sqtag)
        logger.info("[%d] Starting squaring phase", self.my_id)
        self.out_shares = {0: t_share}
        self.powers = {0: t_pk}
        self.th_powers = {0:tpks}
        self.shares, self.d_shares, self.low_f_commits, self.high_f_commits = params
        
        msg_buff = [[] for _ in range(self.logq)]
        cur_share = t_share
        for idx in range(1, self.logq+1):
            # Computing the reveal message
            self.inc_idx = False
            cur_share_sq = cur_share*cur_share
            reveal_msg =  cur_share_sq - self.d_shares[idx-1]
            # Proof of equality of discrete logarithm
            cp = CP(self.g, self.g**cur_share, self.ZR, self.multiexp)
            y = self.g**cur_share_sq
            pf = cp.prove(cur_share, self.g**cur_share, y)
            for i in range(self.n):
                send(i, (idx, reveal_msg, y, pf))

            sq_shares = []
            past_msgs = msg_buff[idx-1]
            for msg in past_msgs:
                sender, s_reveal, s_y, s_pf = msg
                if sender == self.my_id:
                    sq_shares.append([sender+1, s_reveal])
                elif self.verify_sq(idx, sender, s_reveal, s_y, s_pf):
                    sq_shares.append([sender+1, s_reveal])
                    if len(sq_shares) > 2*self.t:
                        cur_share = self.process(idx, sq_shares)
                        break
            
            if self.inc_idx:
                continue

            while True:
                (sender, msg) = await recv()
                s_idx, s_reveal, s_y, s_pf = msg

                if s_idx > idx:
                    msg_buff[s_idx-1].append((sender, s_reveal, s_y, s_pf))
                elif (s_idx == idx):
                    if sender == self.my_id:
                        sq_shares.append([sender+1, s_reveal])
                    elif self.verify_sq(idx, sender, s_reveal, s_y, s_pf):
                        sq_shares.append([sender+1, s_reveal])
                        if len(sq_shares) > 2*self.t:
                            cur_share = self.process(idx, sq_shares)
                            break

        assert len(self.out_shares) == self.logq+1
        assert len(self.th_powers) == self.logq+1
        assert len(self.powers ) == self.logq+1
        
        return (self.out_shares, self.th_powers, self.powers)
        