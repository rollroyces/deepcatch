#!/usr/bin/env python3
"""Export a self-contained diagram HTML to a PNG of its <svg> alone.

Why this script exists
----------------------
The two schematic figures under ``docs/figures/`` (the encoder
architecture and the validation protocol) are authored as
self-contained HTML files with an inline SVG, because that is what
makes them editable, reviewable and accessible. The manuscript, the
README and the docs site need a raster file as well, and a raster file
is an export -- so the export has to be reproducible, not a one-off
command somebody ran once.

What it does
-------------
It renders the **original HTML** (not an extracted SVG), so the fonts
load the way the page intends, and screenshots only the ``<svg>``
element's bounding box. The editorial header and footer are dropped by
design: the deliverable is the diagram.

The PNG keeps a transparent background, so it can sit on any page or
slide colour without a white halo. The pixel size is the SVG's
``viewBox`` multiplied by the device scale factor -- pass 2 for docs and
slides, 3 for print.

Usage
-----
    python3 scripts/export_diagram_png.py <figure.html> <figure.png> [scale]

Requires Playwright with Chromium:

    pip install playwright && playwright install chromium
"""
from __future__ import annotations

import pathlib
import sys


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    src, out = sys.argv[1], sys.argv[2]
    scale = int(sys.argv[3]) if len(sys.argv) > 3 else 2

    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print("Playwright isn't installed. To enable PNG export, run:\n"
              "    pip install playwright\n"
              "    playwright install chromium\n"
              "then export again.", file=sys.stderr)
        return 1

    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(device_scale_factor=scale)
        page.goto(f"file://{pathlib.Path(src).resolve()}")
        page.wait_for_load_state("networkidle")
        page.locator("svg").first.screenshot(path=out, omit_background=True)
        browser.close()
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
