#include "utest.h"
#include "InputFile.h"
#include "typses.h"
#include <string>

UTEST(InputFile, BasicLookupAndGet) {
    InputFile inp;
    inp.add("one", "1");
    int oneInInputFile = 0;
    EXPECT_TRUE(inp.lookup("one", oneInInputFile));
    EXPECT_EQ(1, oneInInputFile);
    EXPECT_EQ(1, inp.getOr("one", 2));
    EXPECT_EQ(2, inp.getOr("notthere", 2));
}

UTEST(InputFile, Dbl3AndRawString) {
    InputFile inp;
    inp.add("dbl3_235", " 2 3 5.");
    dbl3 val = inp.getOr("dbl3_235", dbl3());
    EXPECT_NEAR(2.0, val.x, 1e-5);
    EXPECT_NEAR(3.0, val.y, 1e-5);
    EXPECT_NEAR(5.0, val.z, 1e-5);
    EXPECT_STREQ("2", inp.getOr("dbl3_235", std::string()).c_str());
    EXPECT_STREQ(" 2 3 5.", inp["dbl3_235"].c_str());
}
