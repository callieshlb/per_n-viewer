#!/usr/bin/env python3
"""E95: lossless packing of the chart grids for download (finding A of SITE_REVIEW.md). Each part: the float32 grids viewed as
little-endian uint32 words, differenced along the grid's x direction per channel (word i minus word i - 4 inside each row of
S texels x 4 channels; the first texel of a row is kept), byte-shuffled (all byte 0, then all byte 1, ...), gzip (level 9,
mtime 0). The page undoes it (DecompressionStream + unshuffle + running sums) and gets the same bytes. Writes chart_grids_<i>.dxs
(<= 200 charts each) and adds parts_dxs / dxs to chart_meta.json; the raw parts stay (fallback for browsers without
DecompressionStream). Verifies the round trip byte for byte. Raw values only."""
import gzip, hashlib, io, json, os
import numpy as np
M = json.load(open('chart_meta.json')); S, NC = M['S'], M['NC']
raw = b''.join(open(f, 'rb').read() for f in M['parts']); per = S*S*16
assert len(raw) == NC*per
CPP = 200; out = []; tot = 0
for i, c0 in enumerate(range(0, NC, CPP)):
    blk = raw[c0*per:(c0 + CPP)*per]; n = len(blk)//per
    u = np.frombuffer(blk, '<u4').reshape(n, S, S, 4); d = u.copy(); d[:, :, 1:, :] -= u[:, :, :-1, :]
    sh = np.frombuffer(d.tobytes(), np.uint8).reshape(-1, 4).T.copy().tobytes()
    bio = io.BytesIO()
    with gzip.GzipFile(fileobj=bio, mode='wb', compresslevel=9, mtime=0) as g: g.write(sh)
    z = bio.getvalue(); fn = f'chart_grids_{i}.dxs'; open(fn, 'wb').write(z); out.append(fn); tot += len(z)
    # round trip
    b2 = np.frombuffer(gzip.decompress(z), np.uint8).reshape(4, -1).T.copy().view('<u4').reshape(n, S, S, 4).copy()
    np.cumsum(b2, axis=2, dtype=np.uint32, out=b2)
    assert b2.tobytes() == blk, 'round trip failed'
    print(f'   {fn}: charts {c0}..{c0 + n - 1}, raw {len(blk)/1e6:.1f} MB -> {len(z)/1e6:.1f} MB', flush=True)
M['parts_dxs'] = out; M['dxs'] = dict(charts_per_part=CPP, sha256_raw=hashlib.sha256(raw).hexdigest())
json.dump(M, open('chart_meta.json', 'w'))
smp = np.frombuffer(raw, '<f4')[::9973].astype(np.float64).sum()
print(f'packed {NC} charts: raw {len(raw)/1e6:.1f} MB in {len(M["parts"])} parts -> {tot/1e6:.1f} MB in {len(out)} parts (ratio {tot/len(raw):.3f}); round trip identical; '
      f'sha256 of the raw grids {M["dxs"]["sha256_raw"][:16]}...; sample sum (every 9973rd float) {smp!r}')
