import asyncio
from adkg.polynomial import polynomials_over
from adkg.utils.misc import wrap_send, subscribe_recv
from adkg.utils.poly_misc import interpolate_g1_at_x
from adkg.extra_proofs import PoK


import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

class MsgType:
    PREKEY = "P"
    KEY = "K"

# Move the code for random double sharing to here
class RANDOUSHA:
    def __init__(self, gs, h, n, t, logq, my_id, send, recv, curve_params, matrices):
        self.gs, self.h = gs, h
        self.n, self.t, self.logq, self.my_id = (n, t, logq, my_id)
        self.q = 2**self.logq
        # Total number of secrets: 1 for ACS, 1 for tau, logq*(1+2) for random double sharing
        self.sc = 3*self.logq + 2 
        self.send, self.recv = send, recv
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

    async def randousha(self, mks, acss_outputs):
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
                if node in mks:
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
                    if node in mks:
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
        keytag = MsgType.PREKEY
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

        keytag = MsgType.KEY
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
        self.output_queue.put_nowait((mks, t_share, pk, pk_shares, shares, d_shares, low_f_commits, high_f_commits))
        return