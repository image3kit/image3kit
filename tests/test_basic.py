from __future__ import annotations

import numpy as np

import image3kit as ik


def test_version():
    assert ik.__version__ == "0.0.1"

def test_voxlibI():
    img = ik.VxlImgU8((20, 20, 1), 22)

    assert np.array(img, copy=False)[0,0,0] == 22
    assert ik.VxlImgU8((20, 20, 1), 22).data()[0,0,0] == 22

def test_writePng():
    lnt, rad = 128, 16
    img = ik.VxlImgU8((lnt, rad*2, rad*2), 1)
    img.paint(ik.cylinder((0,rad,rad), (lnt, rad, rad), (rad*3)//4, 0))
    img.printInfo()
    img.plotSlice(filename="cross_x.png", normal_axis='x', slice_index=lnt//2, min_val=0, max_val=2, color_map="RGB")
    img.plotSlice(filename="cross_y.png", normal_axis='y', slice_index=rad, min_val=0, max_val=2, color_map="RGB")
    img.plotSlice(filename="cross_z.png", normal_axis='z', slice_index=rad, min_val=0, max_val=2, color_map="RGB")
    indices = np.where(img.data()==0)
    print(indices)
    print(lnt * rad * rad * 4)


if __name__ == "__main__":
    # test_version()
    # test_voxlibI()
    test_writePng()
