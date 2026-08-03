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

UTEST(Typses, PathHelpersWindowsSeparators) {
    // Windows absolute paths use only backslashes; these helpers must not
    // treat such paths as bare filenames (regression for the CI failure
    // where slice PNGs were not written on windows-latest runners).
    std::string winPath = "C:\\Users\\RUNNER~1\\AppData\\Local\\Temp\\tmpXXXX\\Zcylinder.png";

    std::string prepended = prepend("grey_", winPath);
    std::string nameOfPrepended = nameOf(prepended);
    std::string nameOfWinPath = nameOf(winPath);
    std::string basePathOfWinPath = basePath(winPath);

    EXPECT_STREQ("grey_Zcylinder", nameOfPrepended.c_str());
    EXPECT_STREQ("C:\\Users\\RUNNER~1\\AppData\\Local\\Temp\\tmpXXXX\\grey_Zcylinder.png",
                 prepended.c_str());
    EXPECT_STREQ("Zcylinder", nameOfWinPath.c_str());
    EXPECT_STREQ("C:\\Users\\RUNNER~1\\AppData\\Local\\Temp\\tmpXXXX\\Zcylinder",
                 basePathOfWinPath.c_str()); // strips ".png" suffix, keeps the backslash directory intact
}
