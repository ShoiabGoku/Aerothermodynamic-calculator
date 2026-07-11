# Generate PWA icons: dark theme, blunt body + bow shock motif.
from PIL import Image, ImageDraw
import os

def make(size, maskable, fname):
    S = size
    img = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    if maskable:
        d.rectangle([0, 0, S, S], fill=(10, 14, 23, 255))          # full-bleed
        pad = int(S * 0.14)                                        # safe zone
    else:
        d.rounded_rectangle([0, 0, S - 1, S - 1], radius=int(S * 0.18), fill=(10, 14, 23, 255))
        pad = int(S * 0.10)
    w = S - 2 * pad
    lw = max(2, S // 48)

    cx = pad + int(w * 0.62)      # body nose x
    cy = S // 2
    R  = int(w * 0.16)            # nose radius

    # bow shock: family of hyperbola-like arcs upstream of the nose
    shock_x = cx - int(R * 1.9)
    for i, (col, alpha) in enumerate([((255, 120, 73), 255), ((61, 169, 252), 170), ((123, 223, 242), 110)]):
        off = i * int(S * 0.055)
        pts = []
        for t in range(-100, 101):
            yy = t / 100.0
            xx = (shock_x - off) + (yy * yy) * w * 0.34   # hyperbolic shock shape
            pts.append((xx, cy + yy * w * 0.46))
        d.line(pts, fill=col + (alpha,), width=lw + (2 if i == 0 else 0))

    # blunt body (capsule): nose circle + flared cone
    body_col = (230, 237, 246, 255)
    d.ellipse([cx - R, cy - R, cx + R, cy + R], fill=body_col)
    d.polygon([(cx, cy - R), (cx + int(w * 0.30), cy - int(R * 1.7)),
               (cx + int(w * 0.30), cy + int(R * 1.7)), (cx, cy + R)],
              fill=body_col)
    # heat glow on the nose
    d.ellipse([cx - R, cy - R, cx - R // 3, cy + R], fill=(255, 120, 73, 200))

    img.save(fname)
    print("wrote", fname, size)

here = os.path.dirname(os.path.abspath(__file__))
make(192, False, os.path.join(here, "icon-192.png"))
make(512, False, os.path.join(here, "icon-512.png"))
make(512, True,  os.path.join(here, "icon-512-maskable.png"))
