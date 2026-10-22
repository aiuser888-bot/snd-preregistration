/* Second implementation of protocol 3.1 (block angles), in C, written from the text of the protocol (items 1-8) and not translated from the Python code.
 *
 * Implements: protocol 3.1 (normative algorithm), 12.5.1 step 5 (E1: a bit-exact recomputation in another language)
 * Verified by: tests/test_snd.py::test_c_implementation_matches_python
 *
 * Build:  gcc -O1 -ffp-contract=off -o philox_theta philox_theta.c      (no fused multiply-add: the angle formula is defined with separate roundings)
 * Usage:  philox_theta K0 K1 BLOCK M L [raw]
 *           K0, K1: key words as hexadecimal; BLOCK, M, L: decimal.
 *         Without `raw`: line 1 = alpha as 16 hex digits (binary64 bit pattern), then L lines "perm_i theta_i-bits".
 *         With `raw`: the L angles as 8 little-endian bytes each (item 8), for piping into sha256sum.
 * Only unsigned 64-bit integers, one 128-bit product (unsigned __int128, gcc) and IEEE-754 binary64 arithmetic are used.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef unsigned __int128 u128;

static const uint64_t M0 = 0xD2E7470EE14C6C93ULL, M1 = 0xCA5A826395121157ULL;
static const uint64_t W0 = 0x9E3779B97F4A7C15ULL, W1 = 0xBB67AE8584CAA73BULL;

/* item 1: Philox4x64-10 */
static void philox4x64_10(const uint64_t ctr_in[4], const uint64_t key_in[2], uint64_t out[4]) {
    uint64_t c0 = ctr_in[0], c1 = ctr_in[1], c2 = ctr_in[2], c3 = ctr_in[3];
    uint64_t k0 = key_in[0], k1 = key_in[1];
    for (int r = 0; r < 10; r++) {
        u128 p0 = (u128)M0 * c0, p1 = (u128)M1 * c2;
        uint64_t hi0 = (uint64_t)(p0 >> 64), lo0 = (uint64_t)p0;
        uint64_t hi1 = (uint64_t)(p1 >> 64), lo1 = (uint64_t)p1;
        uint64_t n0 = hi1 ^ c1 ^ k0, n1 = lo1, n2 = hi0 ^ c3 ^ k1, n3 = lo0;
        c0 = n0; c1 = n1; c2 = n2; c3 = n3;
        k0 += W0; k1 += W1;                                   /* modulo 2^64 */
    }
    out[0] = c0; out[1] = c1; out[2] = c2; out[3] = c3;
}

/* item 3: the stream; counter C = (m, b, 0, 0) as a 256-bit little-endian integer, block j uses C + 1 + j (mod 2^256) */
typedef struct { uint64_t key[2]; uint64_t counter[4]; uint64_t buf[4]; int used; int have_reserve; uint32_t reserve; } stream_t;

static void add_small(uint64_t c[4], uint64_t v) {
    for (int i = 0; i < 4 && v; i++) { uint64_t old = c[i]; c[i] += v; v = (c[i] < old) ? 1 : 0; }
}

static void stream_init(stream_t *s, uint64_t k0, uint64_t k1, uint64_t b, uint64_t m) {
    s->key[0] = k0; s->key[1] = k1;
    s->counter[0] = m; s->counter[1] = b; s->counter[2] = 0; s->counter[3] = 0;
    s->used = 4; s->have_reserve = 0; s->reserve = 0;
}

static uint64_t next64(stream_t *s) {
    if (s->used == 4) {
        add_small(s->counter, 1);                              /* the counter is incremented before use */
        philox4x64_10(s->counter, s->key, s->buf);
        s->used = 0;
    }
    return s->buf[s->used++];
}

/* item 4: 32-bit draws, low half first, the high half is kept */
static uint32_t next32(stream_t *s) {
    if (s->have_reserve) { s->have_reserve = 0; return s->reserve; }
    uint64_t w = next64(s);
    s->reserve = (uint32_t)(w >> 32); s->have_reserve = 1;
    return (uint32_t)(w & 0xFFFFFFFFULL);
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: %s K0 K1 BLOCK M L [raw]\n", argv[0]); return 2; }
    uint64_t k0 = strtoull(argv[1], NULL, 16), k1 = strtoull(argv[2], NULL, 16);
    uint64_t b = strtoull(argv[3], NULL, 10), m = strtoull(argv[4], NULL, 10);
    uint64_t L = strtoull(argv[5], NULL, 10);
    int raw = argc > 6 && strcmp(argv[6], "raw") == 0;
    if (L == 0 || L >= (1ULL << 32)) { fprintf(stderr, "L must be in 1..2^32-1\n"); return 2; }
    stream_t s; stream_init(&s, k0, k1, b, m);
    const double two_pi = 6.283185307179586;                   /* the binary64 number nearest to 2 pi (item 5) */
    /* item 5: alpha = fl(2 pi * ((w >> 11) * 2^-53)) from the first word */
    uint64_t w = next64(&s);
    double u = (double)(w >> 11) * 0x1.0p-53;
    volatile double alpha = two_pi * u;
    /* item 6: Fisher-Yates with rejection sampling; the first word is used by alpha, 32-bit draws begin afterwards */
    uint32_t *perm = (uint32_t *)malloc(L * sizeof(uint32_t));
    if (!perm) return 3;
    for (uint64_t i = 0; i < L; i++) perm[i] = (uint32_t)i;
    for (uint64_t i = L - 1; i >= 1; i--) {
        int bits = 0; for (uint64_t t = i; t; t >>= 1) bits++;
        uint64_t mask = ((uint64_t)1 << bits) - 1, j;
        for (;;) { uint32_t r = next32(&s); j = r & mask; if (j <= i) break; }
        uint32_t tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
    }
    /* item 7: theta_i = alpha + (2 pi * pi_i) / L, in this order, each step rounded to nearest */
    uint64_t ab; memcpy(&ab, (const void *)&alpha, 8);
    if (!raw) printf("%016llx\n", (unsigned long long)ab);
    for (uint64_t i = 0; i < L; i++) {
        volatile double prod = two_pi * (double)perm[i];
        volatile double quo = prod / (double)L;
        volatile double th = alpha + quo;
        double thd = th;
        if (raw) {                                             /* item 8: 8 bytes, little endian, whatever the host order */
            uint64_t bits64; memcpy(&bits64, &thd, 8);
            unsigned char by[8];
            for (int k = 0; k < 8; k++) by[k] = (unsigned char)((bits64 >> (8 * k)) & 0xFF);
            fwrite(by, 1, 8, stdout);
        } else {
            uint64_t bits64; memcpy(&bits64, &thd, 8);
            printf("%u %016llx\n", perm[i], (unsigned long long)bits64);
        }
    }
    free(perm);
    return 0;
}
