from adkg.config import HbmpcConfig
from adkg.ipc import ProcessProgramRunner
from adkg.adkg import ADKG
from adkg.poly_commit_hybrid import PolyCommitHybrid
from adkg.utils.poly_misc import get_omega
from pypairing import ZR, G1, G2, blsmultiexp as multiexp, dotprod, blsfft
import asyncio
import time
import logging
import uvloop
import numpy as np

logger = logging.getLogger("benchmark_logger")
logger.setLevel(logging.ERROR)
# Uncomment this when you want logs from this file.
logger.setLevel(logging.NOTSET)

def get_avss_params(n, G1, G2):
    g, h = G1.rand(b'g'), G1.rand(b'h')
    g2 = G2.rand(b'g')
    public_keys, private_keys = [None] * n, [None] * n
    for i in range(n):
        private_keys[i] = ZR.hash(str(i).encode())
        public_keys[i] = g**private_keys[i]
    return g, h, g2, public_keys, private_keys

def gen_vector(t, n):
    deg = 2*t
    coeff_1 = np.array([[ZR(i+1)**j for j in range(t+1)] for i in range(n)])
    coeff_2 = np.array([[ZR(i+1)**j for j in range(t+1, deg+1)] for i in range(n)])
    hm_1 = np.array([[ZR(i+1)**j for j in range(n)] for i in range(t+1)])
    hm_2 = np.array([[ZR(i+1)**j for j in range(n)] for i in range(deg-t)])
    rm_1 = np.matmul(coeff_1, hm_1)
    rm_2 = np.matmul(coeff_2, hm_2)

    return (rm_1.tolist(), rm_2.tolist())

async def _run(peers, n, t, k, my_id, start_time):
    g, h, g2, pks, sks = get_avss_params(n, G1, G2)
    pc = PolyCommitHybrid(g, h, ZR, multiexp)
    logq = k
    omega2 = get_omega(ZR, 2*n, 1729)
    mat1, mat2 = gen_vector(t, n)
    curve_params = (ZR, G1, multiexp, dotprod, blsfft)

    async with ProcessProgramRunner(peers, n, t, my_id) as runner:
        send, recv = runner.get_send_recv("")
        logging.debug(f"Starting ADKG: {(my_id)}")
        logging.debug(f"Start time: {(start_time)}, diff {(start_time-int(time.time()))}")
        benchmark_logger = logging.LoggerAdapter(
           logging.getLogger("benchmark_logger"), {"node_id": my_id}
        )
        with ADKG(pks, sks[my_id], g, h, g2, n, t, logq, my_id, omega2, send, recv, pc, curve_params, (mat1, mat2)) as adkg:
            while True:
                if time.time() > start_time:
                    break
                time.sleep(0.1)
            begin_time = time.time()
            logging.info(f"ADKG start time: {(begin_time)}")
            adkg_task = asyncio.create_task(adkg.run_adkg(begin_time))
            await adkg_task
            adkg.kill()
            adkg_task.cancel()
        bytes_sent = runner.node_communicator.bytes_sent
        for k,v in runner.node_communicator.bytes_count.items():
            logging.info(f"[{my_id}] Bytes Sent: {k}:{v} which is {round((100*v)/bytes_sent,3)}%")
        logging.info(f"[{my_id}] Total bytes sent out aa: {bytes_sent}")

if __name__ == "__main__":
    from adkg.config import HbmpcConfig
    logging.info("Running ADKG ...")
    HbmpcConfig.load_config()

    loop = uvloop.new_event_loop()
    asyncio.set_event_loop(loop)
    try:
        loop.run_until_complete(
            _run(
                HbmpcConfig.peers,
                HbmpcConfig.N,
                HbmpcConfig.t,
                HbmpcConfig.k,
                HbmpcConfig.my_id,
                HbmpcConfig.time,
            )
        )
    finally:
        loop.close()
