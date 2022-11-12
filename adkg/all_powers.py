import asyncio
from collections import defaultdict
from adkg.polynomial import polynomials_over
from adkg.utils.misc import wrap_send, subscribe_recv
from adkg.utils.serilization import Serial
from adkg.utils.poly_misc import interpolate_g1_at_x


import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

# Uncomment this when you want logs from this file.
# logger.setLevel(logging.DEBUG)


class ALL_POWERS:
    #@profile
    def __init__(
            self, g, h, n, t, logq, my_id, send, recv, curve_params
    ):  # (# noqa: E501)
        self.n, self.t, self.logq, self.my_id = n, t, logq, my_id
        self.q = 2**self.logq
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

      
    #@profile
    async def powers(self, t_shares, t_commits, t_powers):
        """
        Generating all powers protocol
        """
        aptag = f"AP"                    
        send, recv = self.get_send(aptag), self.subscribe_recv(aptag)
        logger.debug("[%d] Starting all powers phase", self.my_id)
        self.t_shares, self.t_commmits, self.t_powers = t_shares, t_commits, t_powers

        # TODO: The current implementation sends q group elements in each round
        bitlen  = '{0:0'+str(self.logq)+'b}'
        powers = [self.g]*self.q
        for idx in range(self.logq):
            out_msg = [None]*self.q
            cur_share = t_shares[self.logq-idx-1]
            bit_idx = [1 if int(bitlen.format(a)[idx]) else 0 for a in range(self.q)]
            out_msg = [powers[a]**cur_share if bit_idx[a] else None for a in range(self.q)]
            
            for i in range(self.n):
                send(i, (idx, out_msg))
            
            in_msgs, in_count = [[] for _ in range(self.q)], 0
            while True:
                (sender, msg) = await recv()
                s_idx, s_powers = msg
                assert len(s_powers) == self.q
                
                if s_idx == idx:
                    in_count = in_count+1
                    for a in range(self.q):
                        if bit_idx[a]:
                            in_msgs[a].append([sender+1, s_powers[a]])
                            
                if in_count > self.t:
                    for a in range(self.q):
                        if bit_idx[a]:
                            powers[a] = interpolate_g1_at_x(in_msgs[a], 0, self.G1, self.ZR)
                    break
        self.output_queue.put_nowait(powers)
        return