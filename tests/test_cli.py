"""Tests for the GRMA command-line interface (grma_cli.py)."""

import json
import os
import subprocess
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import grma_cli  # noqa: E402
from grey_meta_v8 import get_bcg_data  # noqa: E402


def _write_csv(tmp_path, header, rows, name="data.csv"):
    p = tmp_path / name
    lines = [",".join(header)] + [",".join(str(c) for c in r) for r in rows]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(p)


# --------------------------------------------------------------------------- #
# load_csv
# --------------------------------------------------------------------------- #
def test_load_csv_autodetect(tmp_path):
    yi, vi, _ = get_bcg_data()
    path = _write_csv(tmp_path, ["study", "yi", "vi"],
                      [(i + 1, y, v) for i, (y, v) in enumerate(zip(yi, vi))])
    loaded_y, loaded_v, label = grma_cli.load_csv(path)
    assert np.allclose(loaded_y, yi)
    assert np.allclose(loaded_v, vi)
    assert "yi" in label and "vi" in label


def test_load_csv_se_column_squared(tmp_path):
    path = _write_csv(tmp_path, ["effect", "se"], [(0.2, 0.3), (0.5, 0.25), (-0.1, 0.4)])
    _, v, label = grma_cli.load_csv(path, se_col="se")
    assert np.allclose(v, np.array([0.3, 0.25, 0.4]) ** 2)
    assert "se=se" in label


def test_load_csv_explicit_columns(tmp_path):
    path = _write_csv(tmp_path, ["my_es", "my_var"], [(0.2, 0.3), (0.5, 0.25), (-0.1, 0.4)])
    y, v, _ = grma_cli.load_csv(path, effect_col="my_es", variance_col="my_var")
    assert np.allclose(y, [0.2, 0.5, -0.1])


def test_load_csv_missing_file():
    with pytest.raises(grma_cli.CLIError):
        grma_cli.load_csv("definitely_not_a_real_file_12345.csv")


def test_load_csv_no_effect_column(tmp_path):
    path = _write_csv(tmp_path, ["a", "b"], [(1, 2), (3, 4)])
    with pytest.raises(grma_cli.CLIError):
        grma_cli.load_csv(path)


def test_load_csv_non_numeric(tmp_path):
    path = _write_csv(tmp_path, ["yi", "vi"], [("abc", 0.1), (0.2, 0.2), (0.3, 0.3)])
    with pytest.raises(grma_cli.CLIError):
        grma_cli.load_csv(path)


def test_load_csv_header_only(tmp_path):
    path = _write_csv(tmp_path, ["yi", "vi"], [])
    with pytest.raises(grma_cli.CLIError):
        grma_cli.load_csv(path)


# --------------------------------------------------------------------------- #
# main() end-to-end via in-process argv
# --------------------------------------------------------------------------- #
def test_main_example_bcg_ok(capsys):
    rc = grma_cli.main(["--example", "bcg", "--seed", "20260213"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "GRMA result" in out
    assert "-0.598953" in out  # regression anchor for the estimate line


def test_main_json_output_matches_estimate(capsys):
    rc = grma_cli.main(["--example", "bcg", "--seed", "20260213", "--json"])
    assert rc == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["estimate"] == pytest.approx(-0.5989525560283318, abs=1e-12)
    assert payload["k"] == 13
    assert payload["method"] == "GRMA"


def test_main_compare_lists_all_methods(capsys):
    rc = grma_cli.main(["--example", "morris", "--seed", "1", "--compare", "--json"])
    assert rc == 0
    payload = json.loads(capsys.readouterr().out)
    methods = {r["method"] for r in payload["comparison"]}
    assert {"IV_FE", "DL_RE", "GRMA", "GRMA_noguard"} <= methods


def test_main_bad_variance_returns_2(tmp_path, capsys):
    path = _write_csv(tmp_path, ["yi", "vi"], [(0.2, -0.3), (0.5, 0.25), (0.1, 0.4)])
    rc = grma_cli.main([path])
    assert rc == 2
    assert "error" in capsys.readouterr().err.lower()


def test_main_conflicting_input_and_example_returns_2(tmp_path, capsys):
    path = _write_csv(tmp_path, ["yi", "vi"], [(0.2, 0.3), (0.5, 0.25), (0.1, 0.4)])
    rc = grma_cli.main([path, "--example", "bcg"])
    assert rc == 2


def test_main_bad_conf_returns_2():
    rc = grma_cli.main(["--example", "bcg", "--conf", "1.5"])
    assert rc == 2


def test_main_writes_out_csv(tmp_path):
    out = tmp_path / "res.csv"
    rc = grma_cli.main(["--example", "bcg", "--seed", "1", "--out", str(out), "--json"])
    assert rc == 0
    text = out.read_text(encoding="utf-8")
    assert "estimate" in text
    assert "method" in text


def test_cli_runs_as_subprocess():
    # Exercises the __main__ entry point and exit code plumbing.
    result = subprocess.run(
        [sys.executable, os.path.join(ROOT, "grma_cli.py"), "--example", "bcg", "--seed", "1"],
        cwd=ROOT, capture_output=True, text=True, timeout=120, check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "GRMA result" in result.stdout
