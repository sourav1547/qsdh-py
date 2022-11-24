import asyncio
from adkg.utils.misc import wrap_send, subscribe_recv
from adkg.utils.serilization import Serial
from adkg.utils.poly_misc import interpolate_g1_at_x
from pypairing import G2, pair, blsmultiexp2


import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

# Uncomment this when you want logs from this file.
# logger.setLevel(logging.DEBUG)

class G2_POWERS:
    #@profile
    def __init__(
            self, g, g2, n, t, logq, my_id, send, recv, curve_params
    ):  # (# noqa: E501)
        self.n, self.t, self.logq, self.my_id = n, t, logq, my_id
        self.g, self.g2 = g, g2
        self.ZR, self.G1, self.multiexp, _ = curve_params
        self.sr = Serial(self.G1)
        self.rands = [self.ZR.rand() for _ in range(self.logq+1)]

        self.benchmark_logger = logging.LoggerAdapter(
            logging.getLogger("benchmark_logger"), {"node_id": self.my_id}
        )

        # Create a mechanism to split the `recv` channels based on `tag`
        self.subscribe_recv_task, self.subscribe_recv = subscribe_recv(recv)

        # Create a mechanism to split the `send` channels based on `tag`
        def _send(tag):
            return wrap_send(tag, send)

        self.get_send = _send
        self.output_queue = asyncio.Queue()

    def __enter__(self):
        return self

    def kill(self):
        self.subscribe_recv_task.cancel()

    def verify_share(self, sender, s_powers):                
        g1_pows = [self.t_commits[i][sender+1] for i in range(self.logq+1)]
        g1_rlc = self.multiexp(g1_pows, self.rands) 
        g2_rlc = blsmultiexp2(s_powers, self.rands)
        self.benchmark_logger.info(f"G2 phase pairing verification: {pair(g1_rlc, self.g2) == pair(self.g, g2_rlc)}")  
        return True
    
    #@profile
    async def powers(self, t_shares, t_commits):
        """
        Generating g2 powers of tau squaring phase
        """
        gtag = f"G"                    
        send, recv = self.get_send(gtag), self.subscribe_recv(gtag)
        logger.debug("[%d] Starting all powers phase", self.my_id)
        self.t_shares, self.t_commits, = t_shares, t_commits

        self.benchmark_logger.info(f"G2 phase started")   
        g2powers  = [G2.identity()]*(self.logq+1)

        g2_th_powers = [self.g2**t_shares[i] for i in range(self.logq+1)]
        for node in range(self.n):
            send(node, g2_th_powers)

        self.benchmark_logger.info(f"G2 phase message sent")   
        g2_temps = [[] for _ in range(self.logq+1)]

        while True:
            (sender, s_powers) = await recv()
            self.benchmark_logger.info(f"G2 phase message received from:{sender}, myid: {self.my_id}")   
            if sender == self.my_id:
                for i in range(self.logq+1):
                    g2_temps[i].append([sender+1, s_powers[i]])
                continue
 
            if self.verify_share(sender, s_powers):
                for i in range(self.logq+1):
                    g2_temps[i].append([sender+1, s_powers[i]])
                self.benchmark_logger.info(f"G2 phase len check: {len(g2_temps[0])}")   

                if len(g2_temps[0]) > self.t:
                    for i in range(self.logq+1):
                        g2powers[i] = interpolate_g1_at_x(g2_temps[i], 0, G2, self.ZR)
                    self.benchmark_logger.info(f"G2 powers computed")   
                    break
        
        self.output_queue.put_nowait(g2powers)
        return