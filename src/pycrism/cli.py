"""Command-line interface: ``pycrism list | params | browse``."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from ._version import __version__
from .browse import RECIPES, RECIPE_TITLES, browse_rgb
from .parameters import PARAM_INFO, PARAMS, compute, default_names, stack_maps
from .readers import read_cube, write_envi


def _summary(name: str, arr: np.ndarray) -> str:
    finite = np.isfinite(arr)
    if not finite.any():
        return f"{name:<14} all NaN"
    a = arr[finite]
    return (f"{name:<14} valid {finite.mean() * 100:5.1f} %  min {a.min():+.4f}  "
            f"median {np.median(a):+.4f}  max {a.max():+.4f}")


def cmd_list(args) -> int:
    rows = [(n, i) for n, i in PARAM_INFO.items()
            if args.detector is None or i["detector"] == args.detector]
    print(f"{'name':<14} {'det':<3} {'kind':<17} {'bands [kernels] nm':<52} {'status':<10} description")
    for n, i in rows:
        print(f"{n:<14} {i['detector']:<3} {i['kind']:<17} {i['bands']:<52} {i['status']:<10} {i['description']}")
    print("\nbrowse recipes (R, G, B):")
    for r, ch in RECIPES.items():
        if args.detector is None or all(PARAM_INFO[c]["detector"] == args.detector for c in ch):
            print(f"  {r}  {RECIPE_TITLES.get(r, ''):<22} {', '.join(ch)}")
    return 0


def cmd_params(args) -> int:
    cube, wl, meta = read_cube(args.cube, wa=args.wa)
    if wl is None:
        print("error: no wavelengths found; pass the WA CDR with --wa", file=sys.stderr)
        return 2
    names = args.params or default_names(wl)
    print(f"{Path(args.cube).name}: {cube.shape[0]} lines x {cube.shape[1]} samples x "
          f"{cube.shape[2]} bands, {np.nanmin(wl):.1f}-{np.nanmax(wl):.1f} nm")
    maps = compute(cube, wl, names)
    for n in names:
        print("  " + _summary(n, maps[n]))
    if args.out:
        stack, order = stack_maps(maps, names)
        write_envi(args.out, stack, band_names=order,
                   description=f"pycrism {__version__} summary parameters of {Path(args.cube).name}")
        print(f"wrote {args.out} ({len(order)} bands) and its .hdr")
    return 0


def cmd_browse(args) -> int:
    cube, wl, meta = read_cube(args.cube, wa=args.wa)
    if wl is None:
        print("error: no wavelengths found; pass the WA CDR with --wa", file=sys.stderr)
        return 2
    rgb = browse_rgb(cube, wl, args.recipe, percentiles=(args.low, args.high))
    out = Path(args.out) if args.out else Path(args.cube).with_name(
        f"{Path(args.cube).stem}_{args.recipe.upper()}.png")
    if out.suffix.lower() == ".png":
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("error: writing PNG needs matplotlib (pip install pycrism[plot]); "
                  "use an .img output name for ENVI", file=sys.stderr)
            return 2
        plt.imsave(out, rgb)
    else:
        write_envi(out, rgb, band_names=list(RECIPES[args.recipe.upper()]),
                   description=f"pycrism {__version__} {args.recipe.upper()} browse of {Path(args.cube).name}")
    print(f"wrote {out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="pycrism",
                                description="CRISM summary parameters and browse products in pure NumPy")
    p.add_argument("--version", action="version", version=f"pycrism {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("list", help="list the parameters and browse recipes")
    s.add_argument("--detector", choices=("S", "L"), help="VNIR (S) or IR (L) only")
    s.set_defaults(func=cmd_list)

    s = sub.add_parser("params", help="compute parameters for a cube")
    s.add_argument("cube", help="ENVI .img (with .hdr) or PDS3 .img (with .lbl)")
    s.add_argument("-p", "--params", nargs="+", metavar="NAME", choices=list(PARAMS),
                   help="parameter names (default: all the wavelength axis covers)")
    s.add_argument("-o", "--out", help="write the maps as a float32 BSQ ENVI file")
    s.add_argument("--wa", help="WA wavelength CDR for a PDS3 TRR3 cube")
    s.set_defaults(func=cmd_params)

    s = sub.add_parser("browse", help="write a browse-product RGB image")
    s.add_argument("cube")
    s.add_argument("recipe", type=str.upper, choices=list(RECIPES))
    s.add_argument("-o", "--out", help=".png (needs matplotlib) or .img (ENVI); "
                                        "default <cube>_<RECIPE>.png")
    s.add_argument("--wa", help="WA wavelength CDR for a PDS3 TRR3 cube")
    s.add_argument("--low", type=float, default=2.0, help="lower stretch percentile")
    s.add_argument("--high", type=float, default=98.0, help="upper stretch percentile")
    s.set_defaults(func=cmd_browse)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
