"""Known-answer vectors.  theta_kat.json: protocol 3.1 (theta generation), computed with philox_ref.py only (Python integers, no numpy).
reference_kat.json: protocol 7(D) and 4 (the integer reference value: pi, cos/sin, enclosures of U and U_mon) computed with highprec.py, and protocol 2
(canonical serialisation and digest of f~) computed with hashlib and struct only.

Implements: protocol 3.1 (normative algorithm), 12.5.1 step 5 (E1: an independent implementation must reproduce these), 7(D) (reference value), 2 (H_f~).
Verified by: tests/test_snd.py::test_known_answer_vectors_file_is_current, test_reference_vectors_file_is_current

Usage:  python scripts/make_kat.py            write tests/vectors/theta_kat.json and tests/vectors/reference_kat.json
        python scripts/make_kat.py --check    exit 1 if a file differs from what this script generates
"""
import hashlib, json, math, os, struct, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from snd import philox_ref as P
from snd import highprec as HP

OUT = os.path.join(os.path.dirname(__file__), "..", "tests", "vectors", "theta_kat.json")
REF_OUT = os.path.join(os.path.dirname(__file__), "..", "tests", "vectors", "reference_kat.json")
H = lambda v: "0x%016x" % v
MAXW = (1 << 64) - 1


def generate():
    keys = [(0, 0), (MAXW, MAXW), (0x243F6A8885A308D3, 0x13198A2E03707344), (0x0123456789ABCDEF, 0xFEDCBA9876543210)]
    ctrs = [(0, 0, 0, 0), (MAXW, MAXW, MAXW, MAXW), (0x452821E638D01377, 0xBE5466CF34E90C6C, 0xC0AC29B7C97C50DD, 0x3F84D5B5B5470917),
            (1, 2, 3, 4)]
    doc = {"description": "Known-answer vectors for SND protocol 3.1. Words are unsigned 64-bit integers written as 0x%016x. "
                          "'cipher': Philox4x64-10 block function (Random123), counter -> 4 output words. "
                          "'stream': the first 8 words of the stream for (key, initial counter): the stream uses cipher(counter+1), cipher(counter+2) "
                          "(counter is a 256-bit little-endian integer of four words, +1 modulo 2^256). "
                          "'seed_to_key': SHA-256 of the seed material, key = digest bytes 0..7 and 8..15 as little-endian words. "
                          "'block_theta': alpha, permutation and angles (binary64 bit patterns) of block b at re-draw m. "
                          "'theta_digest': SHA-256 over the little-endian binary64 angles of the listed blocks in ascending order (protocol 12.5.1 step 4).",
           "generated_by": "scripts/make_kat.py using snd/philox_ref.py only"}
    # External value: ISO C++ std::philox4x64 (Philox4x64-10), 10000th output of the default-constructed engine, 3409172418970261260
    # (quoted in the libstdc++ patch discussion https://gcc.gnu.org/pipermail/libstdc++/2025-July/062703.html).  That engine has
    # key (20111115, 0), counter 0 and increments the counter AFTER each block, so output #10000 is word 3 of cipher((2499, 0, 0, 0)).
    ext = P.philox4x64_10((2499, 0, 0, 0), (20111115, 0))
    doc["external_cpp_standard"] = dict(key=[H(20111115), H(0)], counter=[H(2499), H(0), H(0), H(0)], word_index=3, value_dec=3409172418970261260,
                                        value_hex=H(3409172418970261260), computed_hex=H(ext[3]))
    doc["cipher"] = [dict(key=[H(k) for k in key], counter=[H(c) for c in ctr], out=[H(w) for w in P.philox4x64_10(ctr, key)])
                     for key in keys[:3] for ctr in ctrs[:3]]
    st = []
    for key, ctr in ((keys[0], ctrs[0]), (keys[3], ctrs[3]), (keys[2], (MAXW, MAXW, 0, 0)), (keys[1], ctrs[1])):
        s = P.PhiloxStream(key, ctr)
        st.append(dict(key=[H(k) for k in key], counter=[H(c) for c in ctr], words=[H(s.next_u64()) for _ in range(8)]))
    doc["stream"] = st
    sk = []
    for raw in (bytes(range(64)), bytes([0xA5] * 64), bytes(range(200, 256)) + bytes(range(8))):
        for dev, ver in ((1, "2.0.13"), (4, "2.0.13")):
            mat = raw + b"|dev=%d|ver=%s" % (dev, ver.encode())
            d = hashlib.sha256(mat).digest()
            sk.append(dict(seed_material_hex=mat.hex(), sha256=d.hex(), key=[H(k) for k in P.key_from_digest(d)]))
    doc["seed_to_key"] = sk
    bt = []
    for key, b, m, L in ((keys[3], 0, 0, 1), (keys[3], 0, 0, 2), (keys[3], 0, 0, 8), (keys[3], 1, 1, 16), (keys[2], 7, 0, 10),
                         (keys[2], 2 ** 20 + 3, 1999, 33)):
        th, alpha, perm = P.block_theta(key, b, m, L)
        bt.append(dict(key=[H(k) for k in key], block=b, m=m, L=L, alpha_bits=H(P.struct.unpack("<Q", P.struct.pack("<d", alpha))[0]),
                       perm=perm, theta_bits=[H(P.struct.unpack("<Q", P.struct.pack("<d", t))[0]) for t in th]))
    doc["block_theta"] = bt
    td = []
    for key, blocks, L, m in ((keys[3], [0, 1, 2], 10000, 0), (keys[3], [0, 1, 2], 10000, 1), (keys[2], [5, 9, 10, 4000000], 10000, 1999)):
        h = hashlib.sha256()
        for b in sorted(blocks):
            h.update(P.theta_bytes(P.block_theta(key, b, m, L)[0]))
        td.append(dict(key=[H(k) for k in key], blocks=blocks, L=L, m=m, sha256=h.hexdigest()))
    doc["theta_digest"] = td
    return doc


def _unit(tag, i):
    """Deterministic number in [0, 1) with a short binary expansion, from SHA-256 (exactly representable: a 53-bit integer times 2^-53)."""
    return (int.from_bytes(hashlib.sha256(b"%s-%d" % (tag, i)).digest()[:8], "big") >> 11) / float(1 << 53)


def _frac(q):
    return "%d/%d" % (q.numerator, q.denominator)


def generate_reference():
    doc = {"description": "Known-answer vectors for the integer reference value (SND protocol 7(D), 4) and the f~ digest (protocol 2). Floats are written with float.hex(). "
                          "'pi_fixed': floor-style fixed-point pi, value = integer / 2^bits (error <= 2 units). "
                          "'sincos': cos and sin of a binary64 angle as integers in units of 2^-P (error <= 1 unit). "
                          "'reference_U': enclosure (re, im, err) of sum_i w_i e^{i theta_i} (exact fractions; with 'classes': U_mon of protocol 12.2) for the listed w, theta (K x L, C order). "
                          "'ftilde_digest': SHA-256 of the little-endian binary64 serialisation of the rows of f~ in ascending order of the absolute block numbers (protocol 2).",
           "generated_by": "scripts/make_kat.py using snd/highprec.py (integers) and hashlib/struct"}
    doc["pi_fixed"] = [dict(bits=b, value=str(HP.pi_fixed(b))) for b in (64, 128, 256, 1024)]
    angles = [0.0, 1e-10, 0.5, 1.0, math.pi / 2, math.pi, 2 * math.pi, 3 * math.pi / 2, 4 * math.pi, -1.0, 12.566370614359172, 1e5]
    angles += [4 * math.pi * _unit(b"kat-angle", i) for i in range(10)]
    doc["sincos"] = [dict(theta=t.hex(), P=Pb, cos=str(HP.sincos_fixed(t, Pb)[0]), sin=str(HP.sincos_fixed(t, Pb)[1])) for t in angles for Pb in (128,)]
    doc["sincos"] += [dict(theta=t.hex(), P=300, cos=str(HP.sincos_fixed(t, 300)[0]), sin=str(HP.sincos_fixed(t, 300)[1])) for t in angles[:3] + angles[-2:]]
    K, L = 2, 6
    w = [4.0 * _unit(b"kat-w", i) - 2.0 for i in range(K * L)]
    th = [4 * math.pi * _unit(b"kat-theta", i) for i in range(K * L)]
    cls = [i % 3 for i in range(K * L)]
    import numpy as np
    ru = []
    for classes in (None, cls):
        re, im, err = HP.reference_U(np.array(w).reshape(K, L), np.array(th).reshape(K, L), 128, None if classes is None else (np.array(classes).reshape(K, L), 3))
        ru.append(dict(shape=[K, L], w=[x.hex() for x in w], theta=[x.hex() for x in th], classes=classes, P=128, re=_frac(re), im=_frac(im), err=_frac(err)))
    doc["reference_U"] = ru
    Kf, Lf, blocks = 3, 5, [7, 3, 5]
    ft = [4.0 * _unit(b"kat-ft", i) - 2.0 for i in range(Kf * Lf)]
    h = hashlib.sha256()
    for r in sorted(range(Kf), key=lambda r: blocks[r]):
        h.update(b"".join(struct.pack("<d", x) for x in ft[r * Lf:(r + 1) * Lf]))
    doc["ftilde_digest"] = [dict(shape=[Kf, Lf], blocks=blocks, ft=[x.hex() for x in ft], sha256=h.hexdigest())]
    return doc


def _texts():
    return {OUT: json.dumps(generate(), indent=1, sort_keys=True) + "\n", REF_OUT: json.dumps(generate_reference(), indent=1, sort_keys=True) + "\n"}


if __name__ == "__main__":
    texts = _texts()
    if "--check" in sys.argv:
        bad = 0
        for path, text in texts.items():
            ok = os.path.exists(path) and open(path).read() == text
            print(os.path.basename(path), "current" if ok else "STALE")
            bad += not ok
        sys.exit(1 if bad else 0)
    for path, text in texts.items():
        open(path, "w", newline="\n").write(text)
        print("wrote", os.path.normpath(path), len(text), "bytes")
