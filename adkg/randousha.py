import asyncio
from adkg.polynomial import polynomials_over
from adkg.utils.misc import wrap_send, subscribe_recv
from adkg.utils.poly_misc import interpolate_g1_at_x, interpolate_g1_at_all, evaluate_g1_at_all
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

    # Randomness extraction for shares of tau
    def gen_tau_shares(self, acss_outputs):
        secrets = [self.ZR(0)]*self.n
        randomness = [self.ZR(0)]*self.n
        commits = [self.G1.identity()]*self.n

        for node in range(self.n):
            secrets[node] = acss_outputs[node]['shares']['msg'][2]
            randomness[node] = acss_outputs[node]['shares']['rand'][1]
            commits[node] = acss_outputs[node]['commits'][2][0]
        
        z_shares = [self.ZR(0)]*self.n
        r_shares = [self.ZR(0)]*self.n
        for i in range(self.n):
            z_shares[i] = self.dotprod(self.matrix[0][i], secrets)
            r_shares[i] = self.dotprod(self.matrix[0][i], randomness)
        
        return (z_shares, r_shares, commits)
    
    # Randomness extraction for double shares
    def gen_double_shares(self, mks, acss_outputs):

        z_shares = [[self.ZR(0)]*self.logq for _ in range(self.n)]
        r_shares = [[self.ZR(0)]*self.logq for _ in range(self.n)]
        dz_shares = [[self.ZR(0)]*self.logq for _ in range(self.n)]
        dr_shares = [[self.ZR(0)]*self.logq for _ in range(self.n)]
        self.com_z_0 = [None]*self.logq
        
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
                    
                    low_const = low_const + secrets[node]
                    low_const_r = low_const_r + randomness[node]
            
            self.com_z_0[count] = self.multiexp(commits,[self.ZR(1)]*self.n)

            # For (n,2t+1) shares of random values
            d_secrets = [[self.ZR(0)]*self.n for _ in range(2)]
            d_randomness = [[self.ZR(0)]*self.n for _ in range(2)]
            d_commits = [[self.G1.identity()]*self.n for _ in range(2)]

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

        return (z_shares, r_shares, dz_shares, dr_shares)

    async def pre_key(self, z_shares, r_shares, dz_shares, dr_shares):
        # Sending PREKEY messages
        keytag = MsgType.PREKEY
        send, recv = self.get_send(keytag), self.subscribe_recv(keytag)

        for i in range(self.n):
            send(i, (self.tz_shares[i], self.tr_shares[i], z_shares[i], r_shares[i], dz_shares[i], dr_shares[i]))
        
        # To store shares of tau
        tz_share_shares = []
        tr_share_shares = []

        # To store shares of double shares
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

                # TODO: Can we do some kind of batch implementation to speed up things?
                for ii in range(self.logq):
                    shares[ii] = self.poly.interpolate_at(shares_shares[ii], 0)
                    randoms[ii] = self.poly.interpolate_at(randoms_shares[ii], 0)
                    d_shares[ii] = self.poly.interpolate_at(d_shares_shares[ii], 0)
                    d_randoms[ii] = self.poly.interpolate_at(d_randoms_shares[ii], 0)
                
                break
                # TODO(@sourav): 
                # 1. Implement the verification to check correctness of the reconstructed  value
                # 2. Implement the OEC
        return t_share, t_random, shares, randoms, d_shares, d_randoms

    async def derive_key(self, shares, randoms, d_shares, d_randoms):
        mt = self.gs[0]**self.t_share
        mtr = self.h**self.t_random
        gpok = PoK(self.gs[0], self.ZR, self.multiexp)
        hpok = PoK(self.h, self.ZR, self.multiexp)
        gpf = gpok.prove(self.t_share, mt)
        hpf = hpok.prove(self.t_random, mtr)

        # TODO: We will need to generate proof of knowledge for these values as well.
        # Can we possibly prove it in batch?
        low_commits = [self.gs[0]**shares[ii] for ii in range(self.logq)]
        low_r_commits = [self.h**randoms[ii] for ii in range(self.logq)]
        high_commits = [self.gs[0]**d_shares[ii] for ii in range(self.logq)]
        high_r_commits = [self.h**d_randoms[ii] for ii in range(self.logq)]

        low_pfs = [None]*self.logq
        low_r_pfs = [None]*self.logq
        high_pfs = [None]*self.logq
        high_r_pfs = [None]*self.logq
        
        for ii in range(self.logq):
            low_pfs[ii] = gpok.prove(shares[ii], self.gs[0]**shares[ii])
            low_r_pfs[ii] = hpok.prove(randoms[ii],self.h**randoms[ii])
            high_pfs[ii] = gpok.prove(d_shares[ii], self.gs[0]**d_shares[ii])
            high_r_pfs[ii] = hpok.prove(d_randoms[ii], self.h**d_randoms[ii])

        all_commits = (low_commits, low_r_commits, high_commits, high_r_commits)
        all_proofs = (low_pfs, low_r_pfs, high_pfs, high_r_pfs)

        keytag = MsgType.KEY
        send, recv = self.get_send(keytag), self.subscribe_recv(keytag)

        for i in range(self.n):
            send(i, (mt, mtr, gpf, hpf, all_commits, all_proofs))
        
        pk_shares = [[self.my_id+1, mt]]
        rk_shares = [[self.my_id+1, mtr]]

        low_f_commits = [[] for _ in range(self.logq)]
        low_f_randoms = [[] for _ in range(self.logq)]
        high_f_commits = [[] for _ in range(self.logq)]
        high_f_randoms = [[] for _ in range(self.logq)]

        for ii in range(self.logq):
            low_f_commits[ii].append([self.my_id+1, low_commits[ii]])
            low_f_randoms[ii].append([self.my_id+1, low_r_commits[ii]])
            high_f_commits[ii].append([self.my_id+1, high_commits[ii]])
            high_f_randoms[ii].append([self.my_id+1, high_r_commits[ii]])

        while True:
            (sender, msg) = await recv()
            if sender != self.my_id:
                x, y, gpf_s, hpf_s, all_s_commits, all_s_proofs = msg
                low_com_s, low_r_com_s, high_com_s, high_r_com_s = all_s_commits 
                low_pfs_s, low_r_pfs_s, high_pfs_s, high_r_pfs_s = all_s_proofs

                valid_pok = gpok.verify(x, gpf_s) and hpok.verify(y, hpf_s)
                if valid_pok:
                    pk_shares.append([sender+1, x])
                    rk_shares.append([sender+1, y])

                    for ii in range(self.logq):
                        valid_z_pok = gpok.verify(low_com_s[ii], low_pfs_s[ii])
                        valid_z_pok = valid_z_pok and hpok.verify(low_r_com_s[ii], low_r_pfs_s[ii])
                        if valid_z_pok:
                            low_f_commits[ii].append([sender+1, low_com_s[ii]])
                            low_f_randoms[ii].append([sender+1, low_r_com_s[ii]])


                        valid_d_pok = gpok.verify(high_com_s[ii], high_pfs_s[ii])
                        valid_d_pok = valid_d_pok and hpok.verify(high_r_com_s[ii], high_r_pfs_s[ii])
                        if valid_d_pok:
                            high_f_commits[ii].append([sender+1, high_com_s[ii]])
                            high_f_randoms[ii].append([sender+1, high_r_com_s[ii]])


            if len(pk_shares) > 2*self.t:
                break
        pk =  interpolate_g1_at_x(pk_shares, 0, self.G1, self.ZR)
        rk =  interpolate_g1_at_x(rk_shares, 0, self.G1, self.ZR)
        
        com0 = self.multiexp(self.t_commits, [self.ZR(1)]*self.n)
        assert pk*rk == com0

        for i in range(self.logq):
            high_pk = interpolate_g1_at_x(high_f_commits[i], 0, self.G1, self.ZR)
            high_rk = interpolate_g1_at_x(high_f_randoms[i], 0, self.G1, self.ZR)
            assert high_pk*high_rk == self.com_z_0[i]

            low_pk = interpolate_g1_at_x(low_f_commits[i], 0, self.G1, self.ZR)
            low_rk = interpolate_g1_at_x(low_f_randoms[i], 0, self.G1, self.ZR)
            assert low_pk*low_rk == self.com_z_0[i]
        # FIXME! Add the fallback path for the case when this assertion fails
        
        # TODO: To optimize this using NTT
        pk_out_shares = interpolate_g1_at_all(pk_shares, self.n, self.G1, self.ZR)
        low_out_commits =[None]*self.logq
        high_out_commits=[None]*self.logq
        for ii in range(self.logq):
            low_out_commits[ii] = interpolate_g1_at_all(low_f_commits[ii], self.n, self.G1, self.ZR)
            high_out_commits[ii] = interpolate_g1_at_all(high_f_commits[ii], self.n, self.G1, self.ZR)
            
        return (pk, pk_out_shares, low_out_commits, high_out_commits)

    async def randousha(self, mks, acss_outputs):
        # Generating messages for shares of tau
        self.tz_shares, self.tr_shares, self.t_commits = self.gen_tau_shares(acss_outputs)
        
        # Generating messages for double shaers
        z_shares, r_shares, dz_shares, dr_shares  = self.gen_double_shares(mks, acss_outputs)
        
        # Running the Randomness Extraction Phase
        # TODO: Probably we will need to do something similar for the commitments as well.
        pre_key_task = asyncio.create_task(self.pre_key(z_shares, r_shares, dz_shares, dr_shares))
        self.t_share, self.t_random, shares, randoms, d_shares, d_randoms = await pre_key_task

        # Deriving public commitments
        key_task = asyncio.create_task(self.derive_key(shares, randoms, d_shares, d_randoms))
        pk, pk_shares, low_f_commits, high_f_commits = await key_task

        self.output_queue.put_nowait((mks, self.t_share, pk, pk_shares, shares, d_shares, low_f_commits, high_f_commits))
        return