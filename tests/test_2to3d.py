from __future__ import annotations

from pathlib import Path

import voxlib as vx

cwd = Path(__file__).parent
def test_readPng():
    img = vx.VxlImgU8()
    assert img.data().shape == (0,0,0)
    img = vx.readImage(cwd / "piskelapp.png")
    print(img, img.data().shape)

if __name__ == "__main__":
    # test_version()
    # test_voxlibI()
    test_readPng()