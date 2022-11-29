# This function outputs roots of unity
from functools import reduce
import operator

def get_omega(field, n, seed=None):
    """
    Given a field, this method returns an n^th root of unity.
    If the seed is not None then this method will return the
    same n'th root of unity for every run with the same seed

    This only makes sense if n is a power of 2!
    """
    assert n & n - 1 == 0, "n must be a power of 2"
    x = field.rand(str(seed).encode())
    y = x**int(field(-1)/n)
    if y == 1 or y**(n//2) == 1:
        return get_omega(field, n, seed+1)
    assert y**n == 1, "omega must be 2n'th root of unity"
    assert y**(n // 2) != 1, "omega must be primitive 2n'th root of unity"
    return y

# Helper Functions
def lagrange_at_x(s, j, x, ZR):
    s = sorted(s)
    assert j in s
    l1 = [x - jj for jj in s if jj != j]
    l2 = [j - jj for jj in s if jj != j]
    (num, den) = (ZR(1), ZR(1))
    for item in l1:
        num *= item
    for item in l2:
        den *= item
    return num / den


def interpolate_g1_at_x(coords, x, G, ZR, order=-1):
    if order == -1:
        order = len(coords)
    xs = []
    sortedcoords = sorted(coords, key=lambda x: x[0])
    for coord in sortedcoords:
        xs.append(coord[0])
    s = set(xs[0:order])
    out = G.identity()
    for i in range(order):
        out *= (sortedcoords[i][1] ** (lagrange_at_x(s, xs[i], x, ZR)))
    return out

def interpolate_g1_at(shares, x_recomb, multiexp, ZR):
    # shares are in the form (x, y=f(x))
    if type(x_recomb) is int:
        x_recomb = ZR(x_recomb)
    xs, ys = zip(*shares)
    vector = []
    for i, x_i in enumerate(xs):
        factors = [
            #(x_k - x_recomb) / (x_k - x_i) for k, x_k in enumerate(xs) if k != i
            (x_recomb - x_k) / (x_i - x_k) for k, x_k in enumerate(xs) if k != i
        ]
        vector.append(reduce(operator.mul, factors))
    #return sum(map(operator.mul, ys, vector))
    #sum = field(0)
    #for i in map(operator.mul, vector, ys):
        #sum += i
    return multiexp(list(ys), vector)

def interpolate_g1_batch_at(xs, shares, x_recomb, multiexp, ZR):
    if type(x_recomb) is int:
        x_recomb = ZR(x_recomb)
    vector = []
    for i, x_i in enumerate(xs):
        factors = [
            #(x_k - x_recomb) / (x_k - x_i) for k, x_k in enumerate(xs) if k != i
            (x_recomb - x_k) / (x_i - x_k) for k, x_k in enumerate(xs) if k != i
        ]
        vector.append(reduce(operator.mul, factors))
    #return sum(map(operator.mul, ys, vector))
    #sum = field(0)
    #for i in map(operator.mul, vector, ys):
        #sum += i
    return [multiexp(share, vector) for share in shares]


# To optimize this using NTT
def interpolate_g1_at_all(coords, n, G1, ZR, order=-1):
    outputs = [None]*(n+1)
    for (idx, value) in coords:
        outputs[idx] = value
    
    for idx in range(n+1): 
        if outputs[idx] == None:
            outputs[idx] = interpolate_g1_at_x(coords, idx, G1, ZR, order)
    return outputs


# To optimize this using NTT
def prep_for_fft(coords, omega, n, multiexp, ZR):
    lag_coords = []
    for idx in range(n):
        if coords[idx] is not None:
            lag_coords.append([omega**idx, coords[idx]])        
    
    for idx in range(n): 
        if coords[idx] == None:
            coords[idx] = interpolate_g1_at(lag_coords, omega**idx, multiexp, ZR)
    return coords

# To optimize this using NTT
def prep_for_fft_batch(xs, ys, omega, ell, n, multiexp, ZR):
    deg = len(xs)
    outputs = [[None]*n for _ in range(ell)]
    missing_cords = []
    for i in range(0, n):
        if i not in xs:
            missing_cords.append(i)
    
    zs = [omega**(i+1) for i in xs]
    
    lags = {}
    for node in missing_cords:
        x_recomb = omega**(node+1)
        vector = []
        for i, x_i in enumerate(zs):
            factors = [
                (x_recomb - x_k) / (x_i - x_k) for k, x_k in enumerate(zs) if k != i
            ]
            vector.append(reduce(operator.mul, factors))
        lags[node] = vector
    
    for i in range(ell):
        for k, x_k in enumerate(xs):
            outputs[i][x_k] = ys[k][i]

        ys_i = [ys[k][i] for k in range(deg)]
        for node in missing_cords:
            outputs[i][node] = multiexp(ys_i, lags[node])
    
    return outputs

# To optimize this using NTT
def evaluate_g1_at_all(coeffs, n, ZR, multiexp):
    return [evaluate_g1_at_x(coeffs, x, ZR, multiexp) for x in range(1,n+1)]

def evaluate_g1_at_x(coeffs, x, ZR, multiexp):
    powers = [ZR(x**j) for j in range(len(coeffs))]
    return multiexp(coeffs, powers)
