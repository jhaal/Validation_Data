#!/usr/bin/env python3
"""Run installed operating decks and require explicit completion, not exit code alone."""
import argparse
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "non_upload/operating_coefficients"
DRAGON = Path("/Users/alokjha/DRAGON/5.1/Dragon")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--burnups", nargs="+", type=float, default=[0., 7.5, 20.])
    p.add_argument("--prefix", default="OP_CONV_J2")
    p.add_argument("--states-per-deck", type=int, default=21)
    args = p.parse_args()
    out = WORK / "dragon_results"; out.mkdir(exist_ok=True)
    for bu in args.burnups:
        name = f"{args.prefix}_BU{bu:04.1f}".replace(".", "p")
        listing = DRAGON / "Darwin_arm64" / (name+".result")
        old = listing.read_text(errors="replace") if listing.is_file() else ""
        if not listing.is_file() or ">|VRF_COMPLETED" not in old or "CONVERGENCE NOT REACHED" in old:
            print(f"Running {name}", flush=True)
            with (WORK / (name+".driver.log")).open("w") as log:
                subprocess.run([str(DRAGON / "rdragon"), name+".x2m"], cwd=DRAGON, stdout=log, stderr=subprocess.STDOUT)
        text = listing.read_text(errors="replace")
        count = sum(line.startswith(">|OP_RESULT ") for line in text.splitlines())
        if ">|VRF_COMPLETED" not in text or count != args.states_per_deck:
            raise SystemExit(f"Incomplete DRAGON listing: {name}; {count}/{args.states_per_deck} states")
        if "CONVERGENCE NOT REACHED" in text:
            raise SystemExit(f"Unconverged DRAGON listing: {name}; results must not be used")
        shutil.copy2(listing, out / listing.name)
        print(f"Complete {name}: {count} states", flush=True)


if __name__ == "__main__": main()
