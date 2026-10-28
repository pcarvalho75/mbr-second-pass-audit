"""Readers for the DIA-NN binary predicted library (carried from validation_v5/candidate_library.py)."""
import struct
from pathlib import Path


def read_i32(f):
    b = f.read(4)
    if len(b) < 4:
        raise EOFError
    return struct.unpack('<i', b)[0]


def read_str(f):
    n = read_i32(f)
    if n <= 0:
        return ''
    return f.read(n).decode('utf-8', 'replace')


def read_i32_vec(f):
    n = read_i32(f)
    if n <= 0:
        return []
    raw = f.read(n * 4)
    return list(struct.unpack('<' + 'i' * n, raw))


def load_pg(path: Path) -> dict:
    with path.open('rb') as f:
        sentinel = read_i32(f)
        f.seek(12, 1)
        libname = read_str(f)
        fasta = read_str(f)
        protein_count = read_i32(f)
        for _ in range(protein_count):
            _flag = read_i32(f)
            modlist = read_i32(f)
            read_str(f)
            read_str(f)
            read_str(f)
            f.seek(8, 1)
            if modlist > 0:
                f.seek(modlist * 4, 1)
        pg_count = read_i32(f)
        idx = {}
        nonempty_s1 = 0
        s1_eq_s0 = 0
        s1_samples = []
        v0_nonzero = 0
        v1_nonzero = 0
        tree_nonzero = 0
        for _ in range(pg_count):
            tree = read_i32(f)
            s0 = read_str(f)
            s1 = read_str(f)
            s2 = read_str(f)
            v0 = read_i32_vec(f)
            v1 = read_i32_vec(f)
            v2 = read_i32_vec(f)
            if tree > 0:
                tree_nonzero += 1
                f.seek(tree * 4, 1)
            if s1:
                nonempty_s1 += 1
            if s1 == s0:
                s1_eq_s0 += 1
            if v0:
                v0_nonzero += 1
            if v1:
                v1_nonzero += 1
            if len(s1_samples) < 8:
                s1_samples.append({'s0': s0[:80], 's1': s1[:80], 's2': s2[:80], 'v0n': len(v0), 'v1n': len(v1), 'v2n': len(v2), 'tree': tree})
            for i in v2:
                if i not in idx:
                    idx[i] = (s0, s1, s2)
        return {'sentinel': sentinel, 'libname': libname, 'fasta': fasta, 'proteinCount': protein_count, 'pgCount': pg_count, 'mapped': len(idx), 'nonemptyS1': nonempty_s1, 's1EqS0': s1_eq_s0, 'v0Nonzero': v0_nonzero, 'v1Nonzero': v1_nonzero, 'treeNonzero': tree_nonzero, 's1Samples': s1_samples, 'idx': idx}


class Cur:

    def __init__(self, buf):
        self.b = buf
        self.p = 0

    def i32(self):
        v = struct.unpack_from('<i', self.b, self.p)[0]
        self.p += 4
        return v

    def f32(self):
        v = struct.unpack_from('<f', self.b, self.p)[0]
        self.p += 4
        return v

    def skip(self, n):
        self.p += n

    def sbytes(self):
        n = self.i32()
        if n <= 0:
            return b''
        s = self.b[self.p:self.p + n]
        self.p += n
        return s

    def sstr(self):
        return self.sbytes().decode('utf-8', 'replace')


def parse_header(c):
    version = c.i32()
    reserved = c.i32()
    flagA = c.i32()
    flagB = c.i32()
    libname = c.sstr()
    fasta = c.sstr()
    nprot = c.i32()
    for _ in range(nprot):
        c.skip(4)
        nmods = c.i32()
        c.sbytes()
        c.sbytes()
        c.sbytes()
        c.skip(8)
        if nmods > 0:
            c.skip(4 * nmods)
    npg = c.i32()
    for _ in range(npg):
        tree = c.i32()
        c.sbytes()
        c.sbytes()
        c.sbytes()
        for _ in range(3):
            nv = c.i32()
            if nv > 0:
                c.skip(4 * nv)
        if tree > 0:
            c.skip(4 * tree)
    for _ in range(3):
        nv = c.i32()
        for _ in range(nv):
            c.sbytes()
    c.skip(16)
    return (version, reserved, flagA, flagB, libname, fasta, nprot, npg)
