from __future__ import annotations

from pathlib import Path

import image3kit as ik

cwd = Path(__file__).parent
def test_readPng():
    img = ik.VxlImgU8()
    assert img.data().shape == (0,0,0)
    img = ik.readImage(cwd / "piskelapp.png")
    assert img.shape() == img.data().shape
    assert img.shape() == img.data().shape
    img.plotAll()
    print(img, img.data().shape)

if __name__ == "__main__":
    # test_version()
    # test_voxlibI()
    test_readPng()