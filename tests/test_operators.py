from __future__ import annotations

import image3kit as ik


def test_bitwise_not():
    # Test for U8 image (0-255)
    img = ik.VxlImgU8((10, 10, 10), 1)

    # Set some values
    img2 = img.copy()
    data = img2.data()
    data[0, 0, 0] = 0
    data[1, 1, 1] = 255
    data[2, 2, 2] = 100

    assert img2.data()[0, 0, 0] == 0
    assert img2.data()[1, 1, 1] == 255
    assert img2.data()[2, 2, 2] == 100
    assert img2.data()[3, 3, 3] == 1

    img.NOT(img2)

    # Check logical inversion
    assert img.data()[0, 0, 0] == 1
    assert img.data()[1, 1, 1] == 0
    assert img.data()[2, 2, 2] == 0
    assert img.data()[3, 3, 3] == 0

    # Verify original is unchanged
    assert data[0, 0, 0] == 0
