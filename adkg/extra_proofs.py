import hashlib

# Chaum Pedersend protocol for equality of discrete logarithms
class CP:
    def __init__(self, g, h, ZR, multiexp):
        self.g  = g
        self.h = h
        self.ZR = ZR
        self.multiexp = multiexp

    def derive_chal(self, x, y, a1, a2):
        commit = str(x)+str(y)+str(a1)+str(a2)
        try:
            commit = commit.encode()
        except AttributeError:
            pass 
        hs =  hashlib.sha256(commit).digest() 
        return self.ZR.hash(hs)

    def verify(self, x, y, proof):
        chal, res = proof
        a1 = self.multiexp([x, self.g],[chal, res])
        a2 = self.multiexp([y, self.h],[chal, res])

        eLocal = self.derive_chal(x, a1, y, a2)
        return eLocal == chal

    def prove(self, alpha, x, y):
        w = self.ZR.random()
        a1 = self.g**w
        a2 = self.h**w
        e = self.derive_chal(x, a1, y, a2)
        return  e, w - e*alpha # return (challenge, response)

# Schnorr's sigma protocol for proof of knowledge
class PoK:
    def __init__(self, g, ZR, multiexp):
        self.g  = g
        self.ZR = ZR
        self.multiexp = multiexp

    def derive_chal(self, x, a):
        commit = str(x)+str(a)
        try:
            commit = commit.encode()
        except AttributeError:
            pass 
        hs =  hashlib.sha256(commit).digest() 
        return self.ZR.hash(hs)

    def verify(self, x, proof):
        chal, res = proof
        a = self.multiexp([x, self.g],[chal, res])
        eLocal = self.derive_chal(x, a)
        return eLocal == chal

    def prove(self, alpha, x):
        w = self.ZR.rand()
        a = self.g**w
        e = self.derive_chal(x, a)
        return  e, w - e*alpha # return (challenge, response)
    