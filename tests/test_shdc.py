from pathlib import Path
from sdhc.SHCPostProc import SHCPostProc

RESOURCES_PATH = Path("tests") / "resources"


def compute_sdhc_e2e(folder: Path, restart_file: Path):

    output_compact_velocities_file = folder.joinpath("small.simu.vels.dat.compact")
    atomic_velocities_file = folder.joinpath("small.simu.vels.dat")
    frequency_window_width = 0.5e12

    backup_prefix = folder.joinpath("backup")
    output_force_constant_file_prefix = folder.joinpath("force_constants")
    unit_scaling_factor = 1.602e-19 / 1e-20 * 1e4
    md_timestep = 2.5e-15

    postprocessor = SHCPostProc(
        compactVelocityFile=str(output_compact_velocities_file),
        KijFilePrefix=str(output_force_constant_file_prefix),
        dt_md=md_timestep,
        scaleFactor=unit_scaling_factor,
        LAMMPSDumpFile=str(atomic_velocities_file),
        widthWin=frequency_window_width,
        NChunks=20,
        chunkSize=50000,
        backupPrefix=str(backup_prefix),
        LAMMPSRestartFile=str(restart_file),
        reCalcVels=False,
        reCalcFC=False,
    )

    postprocessor.postProcess()
    return postprocessor


def test_compute_sdhc():

    folder = RESOURCES_PATH
    restart_file = RESOURCES_PATH.joinpath("quenched.restart")

    postprocessor = compute_sdhc_e2e(folder=folder, restart_file=restart_file)
    assert postprocessor.oms_fft is not None
    assert postprocessor.SHC_smooth is not None
