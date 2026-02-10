import math
EPS = 1e-9
class gint:
    def __init__(self, r, i=0):
        self.r = int(round(r))
        self.i = int(round(i))
    def is_zero(self):
        return self.r == 0 and self.i == 0
    def __add__(self, o):
        return gint(self.r + o.r, self.i + o.i)
    def __sub__(self, o):
        return gint(self.r - o.r, self.i - o.i)
    def __mul__(self, o):
        return gint(self.r*o.r - self.i*o.i, self.r*o.i + self.i*o.r)
    def __floordiv__(self, o):
        d = o.r**2 + o.i**2
        return gint(round((self.r*o.r + self.i*o.i)/d), round((self.i*o.r - self.r*o.i)/d))
    def __mod__(self, o):
        q = self // o
        return self - q * o
    def __repr__(self):
        return "(" + str(self.r) + "+" + str(self.i) + "j)"
class poly:
    def __init__(self, coeffs):
        self.c = []
        for x in coeffs:
            self.c.append(complex(x))
        while len(self.c) > 1 and abs(self.c[-1]) < EPS:
            self.c.pop()
    def is_zero(self):
        return len(self.c) == 1 and abs(self.c[0]) < EPS
    def degree(self):
        return len(self.c) - 1
    def __add__(self, o):
        L = max(len(self.c), len(o.c))
        res = []
        for i in range(L):
            v1 = 0j
            if i < len(self.c): v1 = self.c[i]
            v2 = 0j
            if i < len(o.c): v2 = o.c[i]
            res.append(v1 + v2)
        return poly(res)
    def __sub__(self, o):
        L = max(len(self.c), len(o.c))
        res = []
        for i in range(L):
            v1 = 0j
            if i < len(self.c): v1 = self.c[i]
            v2 = 0j
            if i < len(o.c): v2 = o.c[i]
            res.append(v1 - v2)
        return poly(res)
    def __mul__(self, o):
        L = len(self.c) + len(o.c) - 1
        res = [0j] * L
        for i in range(len(self.c)):
            for j in range(len(o.c)):
                res[i+j] = res[i+j] + self.c[i] * o.c[j]
        return poly(res)
    def divmod(self, o):
        n = []
        for x in self.c: n.append(x)
        d = o.c
        q = [0j] * max(0, len(n) - len(d) + 1)
        for i in range(len(q)-1, -1, -1):
            q[i] = n[i+len(d)-1] / d[-1]
            for j in range(len(d)):
                n[i+j] = n[i+j] - q[i] * d[j]
        return poly(q), poly(n)
    def __repr__(self):
        if self.is_zero(): return "0"
        s = ""
        for i in range(len(self.c)):
            val = self.c[i]
            if abs(val) > EPS:
                if s != "": s += " + "
                c_s = str(val)
                if abs(val.imag) < EPS: c_s = str(val.real)
                if i == 0: s += c_s
                elif i == 1: s += c_s + "x"
                else: s += c_s + "x^" + str(i)
        return s
class rtnlf:
    def __init__(self, num, den=None):
        if isinstance(num, poly): self.num = num
        else: self.num = poly([num])
        if den is None: self.den = poly([1.0])
        else: self.den = den
    def is_zero(self): return self.num.is_zero()
    def __add__(self, o): return rtnlf(self.num*o.den + o.num*self.den, self.den*o.den)
    def __sub__(self, o): return rtnlf(self.num*o.den - o.num*self.den, self.den*o.den)
    def __mul__(self, o): return rtnlf(self.num*o.num, self.den*o.den)
    def inv(self): return rtnlf(self.den, self.num)
    def __repr__(self):
        if self.is_zero(): return "0"
        if self.den.degree() == 0 and abs(self.den.c[0] - 1.0) < EPS: return str(self.num)
        return "(" + str(self.num) + ") / (" + str(self.den) + ")"
class intd:
    def is_zero(self, a): return a == 0
    def add(self, a, b): return a + b
    def sub(self, a, b): return a - b
    def mul(self, a, b): return a * b
    def div(self, a, b): return a // b
    def divides(self, a, b):
        if a == 0: return b == 0
        return b % a == 0
    def bezout(self, a, b):
        def egcd(a, b):
            if a == 0: return b, 0, 1
            g, y, x = egcd(b % a, a)
            return g, x - (b // a) * y, y
        d, x1, x2 = egcd(a, b)
        return d, x1, x2, a//d, b//d
class gssd:
    def is_zero(self, a): return a.is_zero()
    def add(self, a, b): return a + b
    def sub(self, a, b): return a - b
    def mul(self, a, b): return a * b
    def div(self, a, b): return a // b
    def divides(self, a, b):
        if a.is_zero(): return b.is_zero()
        return (b % a).is_zero()
    def bezout(self, a, b):
        def egcd(a, b):
            if a.is_zero(): return b, gint(0), gint(1)
            q = b // a
            g, y, x = egcd(b % a, a)
            return g, x - q * y, y
        d, x1, x2 = egcd(a, b)
        return d, x1, x2, a // d, b // d
class pld:
    def is_zero(self, a): return a.is_zero()
    def add(self, a, b): return a + b
    def sub(self, a, b): return a - b
    def mul(self, a, b): return a * b
    def div(self, a, b): return a.divmod(b)[0]
    def divides(self, a, b):
        q, r = b.divmod(a)
        return r.is_zero()
    def bezout(self, a, b):
        def egcd(p1, p2):
            if p1.is_zero(): return p2, poly([0]), poly([1])
            q, r = p2.divmod(p1)
            g, y, x = egcd(r, p1)
            return g, x - (q * y), y
        d, x1, x2 = egcd(a, b)
        inv_lead = poly([1.0/d.c[-1]])
        d = d * inv_lead
        x1 = x1 * inv_lead
        x2 = x2 * inv_lead
        return d, x1, x2, a.divmod(d)[0], b.divmod(d)[0]
class rtnd:
    def is_zero(self, a): return a.is_zero()
    def add(self, a, b): return a + b
    def sub(self, a, b): return a - b
    def mul(self, a, b): return a * b
    def val(self, f): return f.den.degree() - f.num.degree()
    def divides(self, a, b):
        if a.is_zero(): return b.is_zero()
        return self.val(b) >= self.val(a)
    def div(self, a, b): return a * b.inv()
    def bezout(self, a, b):
        if self.val(a) <= self.val(b): return a, rtnlf(1), rtnlf(0), rtnlf(1), b * a.inv()
        return b, rtnlf(0), rtnlf(1), a * b.inv(), rtnlf(1)
def snf(A_raw, dom):
    M = len(A_raw)
    N = len(A_raw[0])
    A = []
    for i in range(M):
        row = []
        for j in range(N): row.append(A_raw[i][j])
        A.append(row)
    k = 0
    while k < M and k < N:
        piv = None
        for r in range(k, M):
            for c in range(k, N):
                if not dom.is_zero(A[r][c]):
                    piv = (r, c)
                    break
            if piv: break
        if not piv: break
        r, c = piv
        A[k], A[r] = A[r], A[k]
        for i in range(M): A[i][k], A[i][c] = A[i][c], A[i][k]
        loop = True
        while loop:
            loop = False
            for j in range(k + 1, N):
                if not dom.is_zero(A[k][j]):
                    if not dom.divides(A[k][k], A[k][j]):
                        d, x, y, u, v = dom.bezout(A[k][k], A[k][j])
                        for r_ in range(M):
                            vk, vj = A[r_][k], A[r_][j]
                            A[r_][k] = dom.add(dom.mul(x, vk), dom.mul(y, vj))
                            A[r_][j] = dom.sub(dom.mul(u, vj), dom.mul(v, vk))
                        loop = True
            for i in range(k + 1, M):
                if not dom.is_zero(A[i][k]):
                    if not dom.divides(A[k][k], A[i][k]):
                        d, x, y, u, v = dom.bezout(A[k][k], A[i][k])
                        for j in range(k, N):
                            vk, vi = A[k][j], A[i][j]
                            A[k][j] = dom.add(dom.mul(x, vk), dom.mul(y, vi))
                            A[i][j] = dom.sub(dom.mul(u, vi), dom.mul(v, vk))
                        loop = True
            for j in range(k + 1, N):
                if not dom.is_zero(A[k][j]):
                    q = dom.div(A[k][j], A[k][k])
                    for r_ in range(M): A[r_][j] = dom.sub(A[r_][j], dom.mul(q, A[r_][k]))
            for i in range(k + 1, M):
                if not dom.is_zero(A[i][k]):
                    q = dom.div(A[i][k], A[k][k])
                    for j in range(k, N): A[i][j] = dom.sub(A[i][j], dom.mul(q, A[k][j]))
        k = k + 1
    return A