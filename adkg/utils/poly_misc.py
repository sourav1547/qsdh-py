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


def interpolate_g1_at_x(coords, x, G1, ZR, order=-1):
    if order == -1:
        order = len(coords)
    xs = []
    sortedcoords = sorted(coords, key=lambda x: x[0])
    for coord in sortedcoords:
        xs.append(coord[0])
    s = set(xs[0:order])
    out = G1.identity()
    for i in range(order):
        out *= (sortedcoords[i][1] ** (lagrange_at_x(s, xs[i], x, ZR)))
    return out


def get_g1_coeffs(coords, t, n, G1, ZR, order=-1):
    coeffs = [G1.identity()]*(t+1)
    # FIXME: This is an incorrect implementation.
    # Currently this is just a place holder
    for i in range(t+1):
        coeffs[i] = coords[i][1]
    return coeffs

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
def evaluate_g1_at_all(coeffs, n, ZR, multiexp):
    return [evaluate_g1_at_x(coeffs, x, ZR, multiexp) for x in range(1,n+1)]

def evaluate_g1_at_x(coeffs, x, ZR, multiexp):
    powers = [ZR(x**j) for j in range(len(coeffs))]
    return multiexp(coeffs, powers)
