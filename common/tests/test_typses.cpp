#include "utest.h"
#include "typses.h"
#include "profilers.h"
#include <array>
#include <cmath>

UTEST(Typses, PieceAndDbls) {
    std::array<int, 6> oneto6 = {1, 2, 3, 4, 5, 6};
    EXPECT_EQ(21, piece<int>(&oneto6[0], 6).sum());
    EXPECT_NEAR(21.0, dbls({1., 2., 3., 4., 5., 6.}).sum(), 1e-5);
    EXPECT_NEAR(0.0, mag(dbl3s(5, dbl3(1, 1, 1)).sum() - dbl3(5, 5, 5)), 1e-5);
}

UTEST(Timing, DataCollection) {
    Timing tim;
    tim("Testing InputFile");
    tim("Testing vars<T>");
    tim("Testing Timing");
    EXPECT_EQ(2u, tim.times.size());
}
