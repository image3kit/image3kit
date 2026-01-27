from __future__ import annotations

import numpy as np

import image3kit as ik


def test_version():
    assert ik.__version__ == "0.0.1"

def test_voxlibI():
    img = ik.VxlImgU8((20, 20, 1), 22)

    assert np.array(img, copy=False)[0,0,0] == 22
    assert ik.VxlImgU8((20, 20, 1), 22).data()[0,0,0] == 22

    numpy_array = np.ones((3, 3, 3))
    print(numpy_array)
    img = ik.VxlImgU8((20, 20, 1), 22)
    print(img.data())

    assert np.array(img, copy=False)[0, 0, 0] == 22
    assert ik.VxlImgU8((20, 20, 1), 22).data()[0, 0, 0] == 22


def test_writePng():
    lnt, rad = 128, 16
    img = ik.VxlImgU8((lnt, rad * 2, rad * 2), 1)
    img.Paint(ik.cylinder((0, rad, rad), (lnt, rad, rad), (rad * 3) // 4, 0))
    img.printInfo()
    img.sliceToPng("x", "cross_x.png", lnt // 2, 0, 2, "RGB")
    img.sliceToPng("y", "cross_y.png", rad, 0, 2, "RGB")
    img.sliceToPng("z", "cross_z.png", rad, 0, 2, "RGB")
    indices = np.where(img.data() == 0)
    print(indices)
    print(lnt * rad * rad * 4)


if __name__ == "__main__":
    # test_version()
    # test_voxlibI()
    test_writePng()
