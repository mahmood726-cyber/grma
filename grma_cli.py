"""
GRMA command-line interface
===========================

Run Grey Relational Meta-Analysis on your own data from a CSV file, without
writing any Python. This is a thin, dependency-light wrapper around the
estimator in ``grey_meta_v8.py``; it changes no scientific behaviour and simply
exposes the existing ``GRMA`` / ``compare_methods`` API on the command line.

Input CSV
---------
A comma-separated file with a header row containing an effect-size column and a
variance column. Column names are auto-detected (case-insensitive) from the
common aliases below, or you can name them explicitly with ``--effect-col`` /
``--variance-col``.

    effect  aliases: effect, yi, y, es, estimate, theta, logrr, logor, smd
    variance aliases: variance, vi, v, var, se2, sampling_variance

If you have standard errors instead of variances, pass ``--se-col`` and the
column will be squared to obtain the variance.

Examples
--------
    # Point estimate + bootstrap CI on your own data
    python grma_cli.py mydata.csv

    # Also compare against IV/DL/REML/HK/Huber/t and GRMA-noguard
    python grma_cli.py mydata.csv --compare

    # Run on the embedded BCG example (no input file needed)
    python grma_cli.py --example bcg --compare

    # Reproducible bootstrap, JSON output, write a results CSV
    python grma_cli.py mydata.csv --seed 20260213 --json --out results.csv

Exit codes
----------
    0  success
    2  usage / input error (bad file, missing columns, invalid data)
"""

import argparse
import csv
import json
import sys

import numpy as np

from grey_meta_v8 import (
    GRMA,
    compare_methods,
    get_bcg_data,
    get_morris_data,
)

EFFECT_ALIASES = ("effect", "yi", "y", "es", "estimate", "theta", "logrr", "logor", "smd")
VARIANCE_ALIASES = ("variance", "vi", "v", "var", "se2", "sampling_variance")
SE_ALIASES = ("se", "std_err", "stderr", "standard_error")


class CLIError(Exception):
    """Raised for user-facing usage/input errors (mapped to exit code 2)."""


def _pick_column(header, requested, aliases, kind):
    """Resolve a column name from an explicit request or a list of aliases.

    Returns the resolved header string. Raises CLIError if it cannot be found
    or is ambiguous.
    """
    lower = {h.lower(): h for h in header}
    if requested is not None:
        if requested in header:
            return requested
        if requested.lower() in lower:
            return lower[requested.lower()]
        raise CLIError(
            f"{kind} column {requested!r} not found. Available columns: {list(header)}"
        )
    matches = [lower[a] for a in aliases if a in lower]
    if not matches:
        raise CLIError(
            f"could not auto-detect a {kind} column. Tried aliases {list(aliases)}; "
            f"available columns: {list(header)}. "
            f"Specify one explicitly with the corresponding --*-col flag."
        )
    return matches[0]


def load_csv(path, effect_col=None, variance_col=None, se_col=None):
    """Load effect/variance arrays from a CSV file.

    Returns ``(yi, vi, label)`` where ``label`` records the resolved columns.
    Raises CLIError on any structural or numeric problem.
    """
    try:
        with open(path, "r", newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames is None:
                raise CLIError(f"{path}: file is empty or has no header row")
            header = list(reader.fieldnames)
            rows = list(reader)
    except FileNotFoundError:
        raise CLIError(f"input file not found: {path}")
    except OSError as exc:
        raise CLIError(f"could not read {path}: {exc}")

    if not rows:
        raise CLIError(f"{path}: no data rows found (only a header?)")

    eff_name = _pick_column(header, effect_col, EFFECT_ALIASES, "effect")

    use_se = se_col is not None
    if use_se:
        var_name = _pick_column(header, se_col, SE_ALIASES, "standard-error")
    else:
        # Prefer an explicit variance column; only fall back to SE auto-detect
        # if the user asked for it. Here we auto-detect variance.
        var_name = _pick_column(header, variance_col, VARIANCE_ALIASES, "variance")

    yi, raw = [], []
    for i, row in enumerate(rows, start=2):  # start=2: header is line 1
        try:
            yi.append(float(row[eff_name]))
            raw.append(float(row[var_name]))
        except (TypeError, ValueError):
            raise CLIError(
                f"{path} line {i}: could not parse numeric value "
                f"from columns {eff_name!r}/{var_name!r} (got "
                f"{row.get(eff_name)!r}, {row.get(var_name)!r})"
            )

    yi = np.asarray(yi, dtype=np.float64)
    raw = np.asarray(raw, dtype=np.float64)
    vi = raw**2 if use_se else raw
    label = f"{path} (effect={eff_name}, {'se' if use_se else 'variance'}={var_name})"
    return yi, vi, label


def _fmt(x, nd=6):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{x:.{nd}f}"


def run(yi, vi, args):
    """Execute the requested analysis and return a JSON-serialisable dict."""
    g = GRMA(
        zeta=args.zeta,
        effect_guard=not args.no_guard,
        tukey_c=args.tukey_c,
    )
    # fit() and bootstrap_ci() both validate inputs internally.
    fit = g.fit(yi, vi)
    ci = g.bootstrap_ci(yi, vi, B=args.boot, conf_level=args.conf, seed=args.seed)

    out = {
        "k": int(len(yi)),
        "method": fit["method"],
        "estimate": ci["estimate"],
        "se": ci.get("se"),
        "ci_lo_pct": ci.get("ci_lo_pct"),
        "ci_hi_pct": ci.get("ci_hi_pct"),
        "ci_lo_bca": ci.get("ci_lo_bca"),
        "ci_hi_bca": ci.get("ci_hi_bca"),
        "pvalue": ci.get("pvalue"),
        "conf_level": args.conf,
        "n_boot_ok": ci.get("n_boot_ok"),
        "w_max": fit.get("w_max"),
        "n_eff": fit.get("n_eff"),
    }

    comparison = None
    if args.compare:
        comparison = []
        for r in compare_methods(yi, vi, seed=args.seed):
            comparison.append(
                {
                    "method": r.get("method"),
                    "estimate": r.get("estimate"),
                    "se": r.get("se"),
                    "ci_lo": r.get("ci_lo"),
                    "ci_hi": r.get("ci_hi"),
                }
            )
        out["comparison"] = comparison
    return out


def _print_human(out):
    conf_pct = int(round(out["conf_level"] * 100))
    print("=" * 60)
    print(f"GRMA result  (k={out['k']}, method={out['method']})")
    print("=" * 60)
    print(f"  estimate      : {_fmt(out['estimate'])}")
    print(f"  se            : {_fmt(out['se'])}")
    print(
        f"  {conf_pct}% CI (pct)  : "
        f"[{_fmt(out['ci_lo_pct'], 4)}, {_fmt(out['ci_hi_pct'], 4)}]"
    )
    print(
        f"  {conf_pct}% CI (BCa)  : "
        f"[{_fmt(out['ci_lo_bca'], 4)}, {_fmt(out['ci_hi_bca'], 4)}]"
    )
    print(f"  p-value       : {_fmt(out['pvalue'], 4)}")
    if out.get("w_max") is not None:
        print(f"  w_max         : {_fmt(out['w_max'])}")
    if out.get("n_eff") is not None:
        print(f"  n_eff         : {_fmt(out['n_eff'], 3)}")
    if out.get("comparison"):
        print("\n  Method comparison")
        print("  " + "-" * 56)
        print(f"  {'method':<14}{'estimate':>12}{'ci_lo':>12}{'ci_hi':>12}")
        for r in out["comparison"]:
            print(
                f"  {str(r['method']):<14}"
                f"{_fmt(r['estimate'], 4):>12}"
                f"{_fmt(r['ci_lo'], 4):>12}"
                f"{_fmt(r['ci_hi'], 4):>12}"
            )


def _write_csv(path, out):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["quantity", "value"])
        for key in (
            "k", "method", "estimate", "se", "ci_lo_pct", "ci_hi_pct",
            "ci_lo_bca", "ci_hi_bca", "pvalue", "conf_level",
            "n_boot_ok", "w_max", "n_eff",
        ):
            w.writerow([key, out.get(key)])
        if out.get("comparison"):
            w.writerow([])
            w.writerow(["method", "estimate", "se", "ci_lo", "ci_hi"])
            for r in out["comparison"]:
                w.writerow([r["method"], r["estimate"], r["se"], r["ci_lo"], r["ci_hi"]])


def build_parser():
    p = argparse.ArgumentParser(
        prog="grma_cli.py",
        description="Run Grey Relational Meta-Analysis (GRMA) on a CSV of effect/variance pairs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example: python grma_cli.py --example bcg --compare",
    )
    src = p.add_argument_group("input")
    src.add_argument("input", nargs="?", help="path to input CSV (omit if using --example)")
    src.add_argument(
        "--example",
        choices=("bcg", "morris"),
        help="use a built-in example dataset instead of an input file",
    )
    src.add_argument("--effect-col", help="name of the effect-size column")
    src.add_argument("--variance-col", help="name of the variance column")
    src.add_argument(
        "--se-col",
        help="name of a standard-error column (squared to variance); overrides --variance-col",
    )

    est = p.add_argument_group("estimator")
    est.add_argument("--zeta", type=float, default=0.5, help="grey distinguishing coefficient in (0,1] (default 0.5)")
    est.add_argument("--tukey-c", type=float, default=4.685, help="Tukey bisquare tuning constant (default 4.685)")
    est.add_argument("--no-guard", action="store_true", help="disable the redescending effect guard")

    inf = p.add_argument_group("inference")
    inf.add_argument("--boot", type=int, default=999, help="bootstrap replicates (default 999)")
    inf.add_argument("--conf", type=float, default=0.95, help="confidence level (default 0.95)")
    inf.add_argument("--seed", type=int, default=None, help="RNG seed for reproducible bootstrap")
    inf.add_argument("--compare", action="store_true", help="also fit comparator methods")

    outg = p.add_argument_group("output")
    outg.add_argument("--json", action="store_true", help="print results as JSON")
    outg.add_argument("--out", help="write results to this CSV path")
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if not (0.0 < args.conf < 1.0):
            raise CLIError("--conf must be in (0, 1)")
        if args.boot < 1:
            raise CLIError("--boot must be a positive integer")

        if args.example:
            if args.input:
                raise CLIError("provide either an input file or --example, not both")
            if args.example == "bcg":
                yi, vi, _ = get_bcg_data()
            else:
                yi, vi, _ = get_morris_data()
        elif args.input:
            yi, vi, label = load_csv(
                args.input,
                effect_col=args.effect_col,
                variance_col=args.variance_col,
                se_col=args.se_col,
            )
            if not args.json:
                print(f"Loaded {label}: k={len(yi)}")
        else:
            parser.error("provide an input CSV or use --example")

        out = run(yi, vi, args)
    except CLIError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except ValueError as exc:
        # Raised by GRMA input validation (e.g. non-positive variance).
        print(f"error: invalid data: {exc}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(out, indent=2, default=lambda o: None if (isinstance(o, float) and not np.isfinite(o)) else o))
    else:
        _print_human(out)

    if args.out:
        try:
            _write_csv(args.out, out)
            if not args.json:
                print(f"\nWrote {args.out}")
        except OSError as exc:
            print(f"error: could not write {args.out}: {exc}", file=sys.stderr)
            return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
