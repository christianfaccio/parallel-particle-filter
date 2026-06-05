#!/usr/bin/env python3
"""
Clean 2-panel previewer for the particle-filter run (no 3D clutter).

  TOP  panel: bird's-eye floor plan (X-Y)
  SIDE panel: elevation view        (X-Z)

The room and feature spheres are drawn as simple outlines (the "map" is just
these surfaces). On top we overlay one live LiDAR scan (cyan), the ground-truth
trajectory (green) and the estimated trajectory (red, should cover the green),
and the vehicle position at the snapshot step (white dot).

Dependency-free (stdlib only). Usage: python3 tools/render_scene.py [out.png]
"""
import math, struct, zlib, sys

# --- scene geometry, matching src/world.c -----------------------------------
HX, HY, ZLO, ZHI = 10.0, 10.0, 0.0, 3.0
SPHERES = [(4.0,3.0,1.0,0.8), (-5.0,-2.0,0.9,0.6), (2.0,-6.0,1.5,1.0),
           (-6.0,5.0,0.7,0.5), (6.5,-4.0,2.0,0.7)]            # cx,cy,cz,r

W, H = 900, 1180
BG = (20, 22, 28)

def read_pcd(path):
    pts, started = [], False
    try: f = open(path)
    except OSError: return pts
    for line in f:
        if not started:
            started = line.startswith("DATA"); continue
        s = line.split()
        if len(s) >= 3: pts.append((float(s[0]), float(s[1]), float(s[2])))
    f.close(); return pts

def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "scene.png"
    scan = read_pcd("scan.pcd")
    gt   = read_pcd("gt_traj.pcd")
    est  = read_pcd("est_traj.pcd")
    if not scan and not gt:
        print("No PCD files — run `PF_DUMP=1 ./bin/pf_serial` first."); return
    veh = est[len(est)//4] if est else (gt[len(gt)//4] if gt else (0,0,0))

    img = bytearray(BG * (W*H))
    def px(x, y, col, r=0):
        for dx in range(-r, r+1):
            for dy in range(-r, r+1):
                X, Y = x+dx, y+dy
                if 0 <= X < W and 0 <= Y < H:
                    o = (Y*W + X)*3; img[o],img[o+1],img[o+2] = col
    def line(x0,y0,x1,y1,col):
        dx,dy = abs(x1-x0), abs(y1-y0); sx = 1 if x0<x1 else -1; sy = 1 if y0<y1 else -1
        err = dx-dy
        while True:
            px(x0,y0,col)
            if x0==x1 and y0==y1: break
            e2 = 2*err
            if e2 > -dy: err -= dy; x0 += sx
            if e2 <  dx: err += dx; y0 += sy
    def circle(cx,cy,rad,col):
        for a in range(0,360,3):
            t = math.radians(a); px(cx+int(rad*math.cos(t)), cy+int(rad*math.sin(t)), col)

    # a panel maps world (a,b) -> pixels within [x0,x1]x[y0,y1] (b points up)
    def make_panel(amin,amax,bmin,bmax, x0,y0,x1,y1):
        sc = min((x1-x0)/(amax-amin), (y1-y0)/(bmax-bmin))
        cax, cbx = 0.5*(amin+amax), 0.5*(bmin+bmax)
        cx, cy = 0.5*(x0+x1), 0.5*(y0+y1)
        def T(a,b):
            return int(cx + (a-cax)*sc), int(cy - (b-cbx)*sc)
        return T

    # ---------- TOP panel (X-Y) ----------
    top = make_panel(-12,12, -12,12, 30, 40, W-30, 40 + (W-60))
    x0,y0 = top(-HX,-HY); x1,y1 = top(HX,HY)
    line(x0,y0,x1,y0,(120,130,150)); line(x1,y0,x1,y1,(120,130,150))
    line(x1,y1,x0,y1,(120,130,150)); line(x0,y1,x0,y0,(120,130,150))   # room walls
    for cx,cy,cz,r in SPHERES:
        a,b = top(cx,cy); circle(a,b,int(r* (x1-x0)/(2*HX)),(230,200,60))
    sc_top = (x1-x0)/(2*HX)
    for p in scan:
        a,b = top(p[0],p[1]); px(a,b,(40,200,230))
    for traj,col in ((gt,(40,210,70)),(est,(235,60,60))):
        for i in range(1,len(traj)):
            a0,b0 = top(traj[i-1][0],traj[i-1][1]); a1,b1 = top(traj[i][0],traj[i][1])
            line(a0,b0,a1,b1,col)
    a,b = top(veh[0],veh[1]); px(a,b,(255,255,255),3)

    # ---------- SIDE panel (X-Z) ----------
    side = make_panel(-12,12, -2,5, 30, 40 + (W-60) + 40, W-30, H-30)
    x0,y0 = side(-HX,ZLO); x1,y1 = side(HX,ZHI)
    line(x0,y0,x1,y0,(120,130,150)); line(x1,y0,x1,y1,(120,130,150))
    line(x1,y1,x0,y1,(120,130,150)); line(x0,y1,x0,y0,(120,130,150))   # room cross-section
    rsc = (x1-x0)/(2*HX)
    for cx,cy,cz,r in SPHERES:
        a,b = side(cx,cz); circle(a,b,int(r*rsc),(230,200,60))
    for p in scan:
        a,b = side(p[0],p[2]); px(a,b,(40,200,230))
    for traj,col in ((gt,(40,210,70)),(est,(235,60,60))):
        for i in range(1,len(traj)):
            a0,b0 = side(traj[i-1][0],traj[i-1][2]); a1,b1 = side(traj[i][0],traj[i][2])
            line(a0,b0,a1,b1,col)
    a,b = side(veh[0],veh[2]); px(a,b,(255,255,255),3)

    # PNG encode
    raw = bytearray()
    for y in range(H):
        raw.append(0); raw += img[y*W*3:(y+1)*W*3]
    def ch(tag,data):
        return struct.pack(">I",len(data))+tag+data+struct.pack(">I",zlib.crc32(tag+data)&0xffffffff)
    png = (b"\x89PNG\r\n\x1a\n" + ch(b"IHDR",struct.pack(">IIBBBBB",W,H,8,2,0,0,0)) +
           ch(b"IDAT",zlib.compress(bytes(raw),9)) + ch(b"IEND",b""))
    open(out,"wb").write(png)
    print(f"wrote {out}  (TOP = X-Y plan, SIDE = X-Z elevation; "
          f"scan={len(scan)}, gt={len(gt)}, est={len(est)})")

if __name__ == "__main__":
    main()
