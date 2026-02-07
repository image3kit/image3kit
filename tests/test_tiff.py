from __future__ import annotations

from pathlib import Path

import image3kit as ik

cwd = Path(__file__).parent


def test_readPng():
    img = ik.VxlImgU8()
    assert img.data().shape == (0, 0, 0)
    img = ik.VxlImgU8(cwd / "piskelapp.png")
    Path("fig").mkdir(exist_ok=True)
    img.write("fig/piskelapp.tif")
    assert Path("fig/piskelapp.tif").exists()
    img2 = ik.VxlImgU8("fig/piskelapp.tif")
    assert img2.data().shape == (72, 54, 1)
    print(img2)
    img2.plotSlice(filename="fig/piskelapp2", normal_axis="z")


if __name__ == "__main__":
    # test_version()
    # test_voxlibI()
    test_readPng()
