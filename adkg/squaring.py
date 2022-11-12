import asyncio
from collections import defaultdict
import math
from adkg.polynomial import polynomials_over
from adkg.utils.misc import wrap_send, subscribe_recv
from adkg.utils.serilization import Serial
from adkg.utils.poly_misc import interpolate_g1_at_x


import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

# Uncomment this when you want logs from this file.
# logger.setLevel(logging.DEBUG)


class SquareMessageType:
    DATA = 1
    
class SQUARE:
    #@profile
    def __init__(
            self, g, h, n, t, logq, my_id, send, recv, curve_params
    ):  # (# noqa: E501)
        self.n, self.t, self.logq, self.my_id = n, t, logq, my_id
        self.q = math.pow(2,self.logq)
        self.g, self.h = g, h
        self.ZR, self.G1, self.multiexp, self.dotprod = curve_params
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

        self.sq_status = defaultdict(lambda: True)
        self.poly = polynomials_over(self.ZR)
        self.poly.clear_cache()
        self.output_queue = asyncio.Queue()
        self.tagvars = {}
        self.tasks = []
        self.data = {}

    def __enter__(self):
        return self

    def kill(self):
        self.subscribe_recv_task.cancel()
        for task in self.tasks:
            task.cancel()
        for key in self.tagvars:
            for task in self.tagvars[key]['tasks']:
                task.cancel()

    # async def _process_prv_msg(self, resutls):
    #     return
    #     for idx, values in self.unvalidated_list.items():
    #         if len(values) > self.t:
    #             match_count = 0
    #             match_ids = {}
    #             for sender, _ in values:
    #                 if self.square_tpks[idx][sender]:
    #                     match_count = match_count + 1
    #                     match_ids[idx] = True
    #             if match_count <= self.t:
    #                 continue
            
    #             for sender, sender_reveal in values:
    #                 if sender not in match_ids:
    #                     continue
    #                 sender_ptpk = self.square_tpks[idx][sender]
    #                 gr_t = self.doubles[2][sender]
    #                 gr_2t = self.doubles[3][sender]
    #                 valid = self.pair(g**sender_reveal, self.G2) == self.pair(sender_ptpk, sender_ptpk)*self.pair(gr_2t, self.G2)
    #                 if valid:
    #                     # Store in validated lists
    #                 if self.square_tpks[idx][sender]:
    #     return results

    #@profile    
    async def _process_msg(self, sender, sender_idx, sender_reveal):
        results  = []
        if not self.square_tpks[sender_idx][sender]:
            self.unvalidated_list[sender_idx].append((sender, sender_reveal))
            return None
        sender_ptpk = self.square_tpks[sender_idx][sender]
        gr_t = self.doubles[2][sender]
        gr_2t = self.doubles[3][sender]
        valid = self.pair(g**sender_reveal, self.G2) == self.pair(sender_ptpk, sender_ptpk)*self.pair(gr_2t, self.G2)

        if valid:
            self.square_tpks[sender_idx]['data'][sender] = self.G1**sender_reveal*gr_t
            self.square_tpks[sender_idx]['count'] = self.square_tpks[sender_idx]['count'] + 1
            if self.square_tpks[sender_idx]['count'] > self.t:
                # TODO: Need to do a FFT in the exponent to fill the threshold public keys
                results.append(sender_idx)
                self._process_prv_msg(results)
        return results
        
    #@profile
    async def square(self, t_share, t_pk, tpks, params):
        """
        Squaring protocol
        """
        sqtag = f"SQ"                    
        send, recv = self.get_send(sqtag), self.subscribe_recv(sqtag)
        logger.debug("[%d] Starting squaring phase", self.my_id)
        out_shares = {0: t_share}
        powers = {0: t_pk}
        th_powers = {0:tpks}

        # TODO: The current implementation is very sequential and proceeds in rounds 
        self.shares, self.d_shares, self.low_f_commits, self.high_f_commits = params
        cur_share = t_share
        for idx in range(self.logq):
            reveal_msg = cur_share*cur_share - self.d_shares[idx]
            for i in range(self.n):
                send(i, (idx, reveal_msg))
            
            sq_shares = []
            while True:
                (sender, msg) = await recv()
                s_idx, s_reveal = msg

                if s_idx == idx:
                    sq_shares.append([sender+1, s_reveal])
                    if len(sq_shares) > 2*self.t:
                        reveal = self.poly.interpolate_at(sq_shares, 0)
                        ran_commit = interpolate_g1_at_x(self.high_f_commits[idx], 0, self.G1, self.ZR)
                        cur_share = reveal + self.shares[idx]
                        out_shares[idx+1] = cur_share
                        powers[idx+1] = (self.g**reveal)*ran_commit
                        th_powers_idx = []
                        for item in self.low_f_commits[idx]:
                            node, value = item[0], item[1]
                            th_powers_idx.append([node, (self.g**reveal)*value])
                        th_powers[idx+1] = th_powers_idx
                        break
        self.output_queue.put_nowait((out_shares, th_powers, powers))
        return

                # next = self._process_msg(sender, s_idx, s_reveal, idx)
                # if next:
                #     share = self.shares[idx]
                #     reveal_msg = share*share + self.doubles[idx-1][1]
                #     for i in range(self.n):
                #         send(i, (SquareMessageType.DATA, idx, reveal_msg))
                #     idx=idx+1

                # if idx == self.logq:
                #     self.output_queue.put_nowait((self.shares, self.square_tpks, self.square_pks))
                #     break



        # self.shares = [None*self.logq]
        # self.square_tpks = {1:tpk}
        # self.square_pks = {1:pk}
        # self.shares[0] = share


        
        # Protocol Idea:
        # 1. Store all SQ message: one per sender per index
        # 2. Maintain the list of valid messages
        # 2.1 Delete the invalid shares
        # 3. Upon receiving enough valid messages for current index,
        # 3.1 Compute the next index
        # 3.2 Send the next message
        # 4. Check if everything is complete: clean and return.
        