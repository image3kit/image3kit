from __future__ import annotations

import math
import tempfile
from pathlib import Path

import numpy as np

import image3kit as ik


def test_voxcylinder_tif():
    img = ik.VxlImgU8((20, 20, 20), 1)
    img.paint(ik.cylinder((0, 10, 10), (20, 10, 10), 5, 0))
    img.spacing = ik.dbl3(1e-6, 1e-6, 1e-6)

    expected_porosity = math.pi * 5 * 5 / (20 * 20)
    actual_porosity = np.mean(img.data == 0)
    assert abs(actual_porosity - expected_porosity) < 0.05 * expected_porosity

    with tempfile.TemporaryDirectory() as tmp_dir:
        mhd_path = str(Path(tmp_dir) / "voxcylinder.mhd")
        tif_path = str(Path(tmp_dir) / "voxcylinder.tif")
        img.write(mhd_path)
        img.write(tif_path)

        img_read = ik.read_image(tif_path)
        assert img_read.shape == (20, 20, 20)
        assert np.array_equal(img_read.data, img.data)


def test_voxcylinder_z():
    img = ik.VxlImgU8((20, 20, 20), 1)
    img.paint(ik.cylinder((10, 10, 0), (10, 10, 20), 5, 0))
    img.spacing = ik.dbl3(1e-6, 1e-6, 1e-6)

    expected_porosity = math.pi * 5 * 5 / (20 * 20)
    actual_porosity = np.mean(img.data == 0)
    assert abs(actual_porosity - expected_porosity) < 0.1 * expected_porosity

    with tempfile.TemporaryDirectory() as tmp_dir:
        slice_path = str(Path(tmp_dir) / "Zcylinder.png")
        img.plot_slice(filename=slice_path, normal_axis="z", slice_index=10)
        png_files = list(Path(tmp_dir).glob("*.png"))
        assert len(png_files) > 0


def test_voxcylinder_noisy():
    img = ik.VxlImgU8((20, 20, 20), 1)
    img.paint(ik.cylinder((10, 10, 0), (10, 10, 20), 5, 0))
    img.spacing = ik.dbl3(1e-6, 1e-6, 1e-6)

    # In VxlPro, addSurfNoise 1 1 13 used mask 1<<2 = 4
    img.add_surf_noise(4, 4, 13)
    img.add_surf_noise(4, 4, 3)

    expected_porosity = math.pi * 5 * 5 / (20 * 20)
    actual_porosity = np.mean(img.data == 0)
    assert abs(actual_porosity - expected_porosity) < 0.1 * expected_porosity

    with tempfile.TemporaryDirectory() as tmp_dir:
        slice_path = str(Path(tmp_dir) / "ZcylRough13x2.png")
        img.plot_slice(filename=slice_path, normal_axis="z", slice_index=10)
        png_files = list(Path(tmp_dir).glob("*.png"))
        assert len(png_files) > 0


if __name__ == "__main__":
    test_voxcylinder_tif()
    test_voxcylinder_z()
    test_voxcylinder_noisy()
