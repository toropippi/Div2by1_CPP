//https://www.hanshq.net/files/eprime.cより拝借
#include <iostream>

#ifndef NUM_DECIMALS
#define NUM_DECIMALS 10000ULL
#endif
#define NUM_DIGITS (1 + NUM_DECIMALS)
#define NUM_BITS (NUM_DIGITS * 3322 / 1000 + 999)


uint64_t umullo(uint64_t a, uint64_t b)
{
        return a * b;
}

#if defined(__GNUC__)
void umul(uint64_t a, uint64_t b, uint64_t *p1, uint64_t *p0)
{
        __uint128_t p = (__uint128_t)a * b;
        *p1 = (uint64_t)(p >> 64);
        *p0 = (uint64_t)p;
}

uint64_t umulhi(uint64_t a, uint64_t b)
{
        return (uint64_t)(((__uint128_t)a * b) >> 64);
}
#elif defined(_MSC_VER) && defined(_M_X64)
#include <intrin.h>
void umul(uint64_t a, uint64_t b, uint64_t *p1, uint64_t *p0)
{
        *p0 = _umul128(a, b, p1);
}

uint64_t umulhi(uint64_t a, uint64_t b)
{
        return __umulh(a, b);
}
#else
void umul(uint64_t x, uint64_t y, uint64_t *p1, uint64_t *p0)
{
        uint32_t x0 = (uint32_t)x, x1 = (uint32_t)(x >> 32);
        uint32_t y0 = (uint32_t)y, y1 = (uint32_t)(y >> 32);
        uint64_t p;
        uint32_t res0, res1, res2, res3;

        p = (uint64_t)x0 * y0;
        res0 = (uint32_t)p;
        res1 = (uint32_t)(p >> 32);

        p = (uint64_t)x0 * y1;
        res1 = (uint32_t)(p += res1);
        res2 = (uint32_t)(p >> 32);

        p = (uint64_t)x1 * y0;
        res1 = (uint32_t)(p += res1);
        p >>= 32;
        res2 = (uint32_t)(p += res2);
        res3 = (uint32_t)(p >> 32);

        p = (uint64_t)x1 * y1;
        res2 = (uint32_t)(p += res2);
        res3 += (uint32_t)(p >> 32);

        *p0 = ((uint64_t)res1 << 32) | res0;
        *p1 = ((uint64_t)res3 << 32) | res2;
}

uint64_t umulhi(uint64_t a, uint64_t b)
{
        uint64_t p0, p1;
        umul(a, b, &p1, &p0);
        return p1;
}
#endif

/* Algorithm 2 from Mﾃｶller and Granlund
   "Improved division by invariant integers". */
uint64_t reciprocal_word(uint64_t d)
{
        uint64_t d0, d9, d40, d63, v0, v1, v2, ehat, v3, v4, hi, lo;

        static const uint64_t table[] = {
        /* Generated with:
           for (int i = (1 << 8); i < (1 << 9); i++)
                   printf("0x%03x,\n", ((1 << 19) - 3 * (1 << 8)) / i); */
        0x7fd, 0x7f5, 0x7ed, 0x7e5, 0x7dd, 0x7d5, 0x7ce, 0x7c6, 0x7bf, 0x7b7,
        0x7b0, 0x7a8, 0x7a1, 0x79a, 0x792, 0x78b, 0x784, 0x77d, 0x776, 0x76f,
        0x768, 0x761, 0x75b, 0x754, 0x74d, 0x747, 0x740, 0x739, 0x733, 0x72c,
        0x726, 0x720, 0x719, 0x713, 0x70d, 0x707, 0x700, 0x6fa, 0x6f4, 0x6ee,
        0x6e8, 0x6e2, 0x6dc, 0x6d6, 0x6d1, 0x6cb, 0x6c5, 0x6bf, 0x6ba, 0x6b4,
        0x6ae, 0x6a9, 0x6a3, 0x69e, 0x698, 0x693, 0x68d, 0x688, 0x683, 0x67d,
        0x678, 0x673, 0x66e, 0x669, 0x664, 0x65e, 0x659, 0x654, 0x64f, 0x64a,
        0x645, 0x640, 0x63c, 0x637, 0x632, 0x62d, 0x628, 0x624, 0x61f, 0x61a,
        0x616, 0x611, 0x60c, 0x608, 0x603, 0x5ff, 0x5fa, 0x5f6, 0x5f1, 0x5ed,
        0x5e9, 0x5e4, 0x5e0, 0x5dc, 0x5d7, 0x5d3, 0x5cf, 0x5cb, 0x5c6, 0x5c2,
        0x5be, 0x5ba, 0x5b6, 0x5b2, 0x5ae, 0x5aa, 0x5a6, 0x5a2, 0x59e, 0x59a,
        0x596, 0x592, 0x58e, 0x58a, 0x586, 0x583, 0x57f, 0x57b, 0x577, 0x574,
        0x570, 0x56c, 0x568, 0x565, 0x561, 0x55e, 0x55a, 0x556, 0x553, 0x54f,
        0x54c, 0x548, 0x545, 0x541, 0x53e, 0x53a, 0x537, 0x534, 0x530, 0x52d,
        0x52a, 0x526, 0x523, 0x520, 0x51c, 0x519, 0x516, 0x513, 0x50f, 0x50c,
        0x509, 0x506, 0x503, 0x500, 0x4fc, 0x4f9, 0x4f6, 0x4f3, 0x4f0, 0x4ed,
        0x4ea, 0x4e7, 0x4e4, 0x4e1, 0x4de, 0x4db, 0x4d8, 0x4d5, 0x4d2, 0x4cf,
        0x4cc, 0x4ca, 0x4c7, 0x4c4, 0x4c1, 0x4be, 0x4bb, 0x4b9, 0x4b6, 0x4b3,
        0x4b0, 0x4ad, 0x4ab, 0x4a8, 0x4a5, 0x4a3, 0x4a0, 0x49d, 0x49b, 0x498,
        0x495, 0x493, 0x490, 0x48d, 0x48b, 0x488, 0x486, 0x483, 0x481, 0x47e,
        0x47c, 0x479, 0x477, 0x474, 0x472, 0x46f, 0x46d, 0x46a, 0x468, 0x465,
        0x463, 0x461, 0x45e, 0x45c, 0x459, 0x457, 0x455, 0x452, 0x450, 0x44e,
        0x44b, 0x449, 0x447, 0x444, 0x442, 0x440, 0x43e, 0x43b, 0x439, 0x437,
        0x435, 0x432, 0x430, 0x42e, 0x42c, 0x42a, 0x428, 0x425, 0x423, 0x421,
        0x41f, 0x41d, 0x41b, 0x419, 0x417, 0x414, 0x412, 0x410, 0x40e, 0x40c,
        0x40a, 0x408, 0x406, 0x404, 0x402, 0x400
        };

        d0 = d & 1;
        d9 = d >> 55;
        d40 = (d >> 24) + 1;
        d63 = (d >> 1) + d0;
        v0 = table[d9 - (1 << 8)];
        v1 = (v0 << 11) - (umullo(umullo(v0, v0), d40) >> 40) - 1;
        v2 = (v1 << 13) + (umullo(v1, (1ULL << 60) - umullo(v1, d40)) >> 47);
        ehat = (v2 >> 1) * d0 - umullo(v2, d63);
        v3 = (v2 << 31) + (umulhi(v2, ehat) >> 1);
        umul(v3, d, &hi, &lo);
        v4 = v3 - (hi + d + (lo + d < lo));
        return v4;
}

/* Algorithm 4 from Mﾃｶller and Granlund
   "Improved division by invariant integers".
   Divide u1:u0 by d, returning the quotient and storing the remainder in r.
   v is the approximate reciprocal of d, as computed by reciprocal_word. */
uint64_t div2by1(uint64_t u1, uint64_t u0, uint64_t d, uint64_t *r, uint64_t v)
{
        uint64_t q0, q1;
        umul(v, u1, &q1, &q0);
        q0 = q0 + u0;
        q1 = q1 + u1 + (q0 < u0);
        q1++;
        *r = u0 - umullo(q1, d);
        q1 = (*r > q0) ? q1 - 1 : q1;
        *r = (*r > q0) ? *r + d : *r;
        if (*r >= d) {
                q1++;
                *r -= d;
        }
        return q1;
}
/* Count leading zeros. */

int clz(uint64_t x)
{
        int n = 0;
        while ((x << n) <= UINT64_MAX / 2) n++;
        return n;
}


/* Right-shift that also handles the 64 case. */
uint64_t shr(uint64_t x, int n)
{
        return n < 64 ? (x >> n) : 0;
}

/* Divide n-place integer u by d, yielding n-place quotient q. */
void divnby1(int n, const uint64_t *u, uint64_t d, uint64_t *q)
{
        uint64_t v, k, ui;
        int l, i;
        /* Normalize d, storing the shift amount in l. */
        l = clz(d);
        d <<= l;
        /* Compute the reciprocal. */
        v = reciprocal_word(d);
        /* Perform the division. */
        k = shr(u[n - 1], 64 - l);
        for (i = n - 1; i >= 1; i--) {
                ui = (u[i] << l) | shr(u[i - 1], 64 - l);
                q[i] = div2by1(k, ui, d, &k, v);
        }
        q[0] = div2by1(k, u[0] << l, d, &k, v);
}

/* Multiply n-place integer u by x in place, returning the overflow word. */
uint64_t mulnby1(int n, uint64_t *u, uint64_t x)
{
        uint64_t k, p1, p0;
        int i;
        k = 0;
        for (i = 0; i < n; i++) {
                umul(u[i], x, &p1, &p0);
                u[i] = p0 + k;
                k = p1 + (u[i] < p0);
        }
        return k;
}

/* Compute x * y mod n, where n << s is normalized and
   v is the approximate reciprocal of n << s. */
uint64_t mulmodn(uint64_t x, uint64_t y, uint64_t n, int s, uint64_t v)
{
        uint64_t hi, lo, r;
        umul(x, y, &hi, &lo);
        div2by1((hi << s) | shr(lo, 64 - s), lo << s, n << s, &r, v);
        return r >> s;
}

/* Compute x^p mod n by means of left-to-right binary exponentiation. */
uint64_t powmodn(uint64_t x, uint64_t p, uint64_t n)
{
        uint64_t res, v;
        int i, l, s;
        s = clz(n);
        v = reciprocal_word(n << s);
        res = x;
        l = 63 - clz(p);
        for (i = l - 1; i >= 0; i--) {
                res = mulmodn(res, res, n, s, v);
                if (p & (1ULL << i)) {
                        res = mulmodn(res, x, n, s, v);
                }
        }
        return res;
}













int main()
{
        int n=3;
        //1423890609600*(1<<128)をdで割りたい
        uint64_t u[n]={0,0,1423890609600};
        uint64_t q[n];
        uint64_t d=1423890609601;

        divnby1(n,u,d,q);
        for(int i=0;i<n;i++){
            std::cout<<q[i]<<std::endl;
        }

        uint64_t aa=powmodn(1024,7128435678909,d);
        std::cout<<aa<<std::endl;
        
        return 0;
}