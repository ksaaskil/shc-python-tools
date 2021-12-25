from pathlib import Path
import subprocess
import tempfile

import pytest

from sdhc.utils import create_compact_vels_file

RESOURCES_PATH = Path("tests") / "resources"

VELS_DAT_PATH = RESOURCES_PATH / "small.simu.vels.dat"


@pytest.mark.cpp
def test_compactify():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpd = Path(tmpdir)
        out_path = tmpd / "vels.dat.compact"
        res = subprocess.run(
            ["./scripts/compactify_vels", str(VELS_DAT_PATH), str(out_path)]
        )
        assert res.returncode == 0

        assert out_path.is_file(), f"Expected to exist: {out_path}"

        content = out_path.read_text()
        expected_content = (RESOURCES_PATH / "small.simu.vels.dat.compact").read_text()

        assert content == expected_content


@pytest.mark.cpp
def test_compactify_v2():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpd = Path(tmpdir)
        out_path = tmpd / "vels.dat.compact"
        create_compact_vels_file(velocities_file=VELS_DAT_PATH, output_file=out_path)

        assert out_path.is_file(), f"Expected to exist: {out_path}"

        content = out_path.read_text()
        expected_content = (RESOURCES_PATH / "small.simu.vels.dat.compact").read_text()

        assert content == expected_content
