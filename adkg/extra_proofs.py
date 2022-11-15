import hashlib

# Chaum Pedersend protocol for equality of discrete logarithms
class CP:
    def __init__(self, g, h, ZR, multiexp):
        self.g  = g
        self.h = h
        self.ZR = ZR
        self.multiexp = multiexp

    def dleq_derive_chal(self, x, y, a1, a2):
        commit = str(x)+str(y)+str(a1)+str(a2)
        try:
            commit = commit.encode()
        except AttributeError:
            pass 
        hs =  hashlib.sha256(commit).digest() 
        return self.ZR.hash(hs)

    def dleq_verify(self, x, y, chal, res):
        a1 = self.multiexp([x, self.g],[chal, res])
        a2 = self.multiexp([y, self.h],[chal, res])

        eLocal = self.dleq_derive_chal(x, a1, y, a2)
        return eLocal == chal

    def dleq_prove(self, alpha, x, y):
        w = self.ZR.random()
        a1 = self.g**w
        a2 = self.h**w
        e = self.dleq_derive_chal(x, a1, y, a2)
        return  e, w - e*alpha # return (challenge, response)

# Schnorr's sigma protocol for proof of knowledge
class PoK:
    def __init__(self, g, ZR, multiexp):
        self.g  = g
        self.ZR = ZR
        self.multiexp = multiexp

    def pok_derive_chal(self, x, a):
        commit = str(x)+str(a)
        try:
            commit = commit.encode()
        except AttributeError:
            pass 
        hs =  hashlib.sha256(commit).digest() 
        return self.ZR.hash(hs)

    def pok_verify(self, x, chal, res):
        a = self.multiexp([x, self.g],[chal, res])
        eLocal = self.pok_derive_chal(x, a)
        return eLocal == chal

    def pok_prove(self, alpha, x):
        w = self.ZR.rand()
        a = self.g**w
        e = self.pok_derive_chal(x, a)
        return  e, w - e*alpha # return (challenge, response)
    