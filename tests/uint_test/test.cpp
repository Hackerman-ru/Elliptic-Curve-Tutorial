#include "pch.h"
#include "uint.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/fwd.hpp>

using boost::multiprecision::uint512_t;

#include <random>

using namespace ECG;

static constexpr size_t ConvertCorrectnessN = 10000000;
static constexpr size_t StringCorrectnessN = 10000;
static constexpr size_t ShiftCorrectnessN = 1000;
static constexpr size_t ArithmeticCorrectnessN = 100000;

static constexpr size_t StringTimingN = 10000;
static constexpr size_t ShiftTimingN = 10000000;
static constexpr size_t ArithmeticTimingN = 1000000;
static constexpr size_t DivisionTimingN = 100000;

static uint_t<512> conv(uint512_t value) {
    uint_t<512> result;

    for (size_t i = 16; i > 0; --i) {
        result <<= 32;
        result |= (value >> ((i - 1) * 32)).convert_to<uint32_t>();
    }

    return result;
}

static void comp(uint512_t lhs, uint_t<512> rhs) {
    for (size_t i = 0; i < 16; ++i) {
        uint32_t a = lhs.convert_to<uint32_t>();
        uint32_t b = rhs.convert_to<uint32_t>();
        EXPECT_EQ(a, b);
        lhs >>= 32;
        rhs >>= 32;

        if (lhs == 0 && rhs == 0) {
            break;
        }
    }
}

TEST(StringConversionCorrectness, DecimalStringConversion) {
    uint512_t boost_value = 1;

    for (size_t i = 1; i < StringCorrectnessN; ++i) {
        boost_value *= i;

        if (boost_value == 0) {
            boost_value = 1;
        }

        auto boost_str = boost_value.convert_to<std::string>();
        uint_t<512> my_value(boost_str.c_str());
        auto my_str = my_value.convert_to<std::string>();
        EXPECT_EQ(boost_str, my_str);
    }
}

TEST(StringConversionTiming, DecimalStringConversionBoost) {
    uint512_t a = 1;

    for (size_t i = 1; i < StringTimingN; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }

        std::string a_str = a.convert_to<std::string>();
    }
}

TEST(StringConversionTiming, DecimalStringConversion) {
    uint_t<512> a = 1;

    for (size_t i = 1; i < StringTimingN; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }

        std::string a_str = a.convert_to<std::string>();
    }
}

// Convert_to

TEST(IntConversionCorrectness, uint32_t) {
    std::mt19937 gen;

    for (size_t i = 0; i < ConvertCorrectnessN; ++i) {
        uint32_t a = gen();
        uint_t<512> b(a);
        EXPECT_EQ(a, b.convert_to<uint32_t>());
    }
}

TEST(IntConversionCorrectness, size_t) {
    std::mt19937_64 gen;

    for (size_t i = 0; i < ConvertCorrectnessN; ++i) {
        size_t a = gen();
        uint_t<512> b(a);
        EXPECT_EQ(a, b.convert_to<size_t>());
    }
}

// Shift

TEST(ShiftCorrectness, LeftShift) {
    uint512_t a("9999999999999999999999999999999999");

    for (size_t i = 1; i < ShiftCorrectnessN; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }

        uint_t<512> b(a.convert_to<std::string>().c_str());

        for (size_t j = 0; j < 512; ++j) {
            comp(a << j, b << j);

            if (a << j == 0) {
                break;
            }
        }
    }
}

TEST(ShiftCorrectness, RightShift) {
    uint512_t a("9999999999999999999999999999999999");

    for (size_t i = 1; i < ShiftCorrectnessN; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }

        uint_t<512> b(a.convert_to<std::string>().c_str());

        for (size_t j = 0; j < 512; ++j) {
            comp(a >> j, b >> j);

            if (a >> j == 0) {
                break;
            }
        }
    }
}

TEST(ShiftTiming, LeftShift) {
    uint_t<512> a("9999999999999999999999999999999999");

    for (size_t i = 0; i < ShiftTimingN; ++i) {
        auto b = a << i;
    }
}

TEST(ShiftTiming, LeftShiftBoost) {
    uint512_t a("9999999999999999999999999999999999");

    for (size_t i = 0; i < ShiftTimingN; ++i) {
        auto b = a << i;
    }
}

TEST(ShiftTiming, RightShift) {
    uint_t<512> a("9999999999999999999999999999999999");

    for (size_t i = 0; i < ShiftTimingN; ++i) {
        auto b = a >> i;
    }
}

TEST(ShiftTiming, RightShiftBoost) {
    uint512_t a("9999999999999999999999999999999999");

    for (size_t i = 0; i < ShiftTimingN; ++i) {
        auto b = a >> i;
    }
}

// Arithmetic
TEST(ArithmeticCorrectness, Addition) {
    uint512_t a = 1;
    uint_t<512> b = 1;

    for (size_t i = 0; i < ArithmeticCorrectnessN; ++i) {
        a += i;
        b += i;

        comp(a, b);
    }
}

TEST(ArithmeticCorrectness, Subtraction) {
    uint512_t a = 1;
    uint_t<512> b = 1;

    for (size_t i = 0; i < ArithmeticCorrectnessN; ++i) {
        a -= i;
        b -= i;

        comp(a, b);
    }
}

TEST(ArithmeticCorrectness, Multiplication) {
    uint512_t a = 1;
    uint_t<512> b = 1;

    for (size_t i = 1; i < ArithmeticCorrectnessN; ++i) {
        a *= i;
        b *= i;

        if (a == 0) {
            a = 1;
        }

        if (b == 0) {
            b = 1;
        }

        comp(a, b);
    }
}

TEST(ArithmeticCorrectness, Division) {
    uint512_t a = 1;
    uint_t<512> b = 1;

    for (size_t i = 1; i < ArithmeticCorrectnessN; ++i) {
        a *= i;
        b *= i;

        if (a == 0) {
            a = 1;
        }

        if (b == 0) {
            b = 1;
        }

        uint512_t c = 1;
        uint_t<512> d = 1;

        for (size_t j = 93; j < ArithmeticCorrectnessN >> 10; j += 3) {
            c *= j;
            d *= j;

            if (c == 0) {
                c = 1;
            }

            if (d == 0) {
                d = 1;
            }

            comp(a / c, b / d);
        }
    }
}

TEST(ArithmeticTiming, Addition) {
    uint_t<512> b = 1;

    for (size_t i = 0; i < ArithmeticTimingN; ++i) {
        b += i;
    }
}

TEST(ArithmeticTiming, AdditionBoost) {
    uint512_t a = 1;
    uint_t<512> b = 1;

    for (size_t i = 0; i < ArithmeticTimingN; ++i) {
        a += i;
    }
}

TEST(ArithmeticTiming, Subtraction) {
    uint_t<512> b = 1;

    for (size_t i = 0; i < ArithmeticTimingN; ++i) {
        b -= i;
    }
}

TEST(ArithmeticTiming, SubtractionBoost) {
    uint512_t a = 1;
    uint_t<512> b = 1;

    for (size_t i = 0; i < ArithmeticTimingN; ++i) {
        a -= i;
    }
}

TEST(ArithmeticTiming, Multiplication) {
    uint_t<512> a = 1;

    for (size_t i = 1; i < ArithmeticTimingN; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }
    }
}

TEST(ArithmeticTiming, MultiplicationBoost) {
    uint512_t a = 1;

    for (size_t i = 1; i < ArithmeticTimingN; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }
    }
}

TEST(ArithmeticTiming, Division) {
    uint_t<512> a = 1;

    for (size_t i = 1; i < DivisionTimingN; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }

        uint_t<512> b = 1;

        for (size_t j = 93; j < DivisionTimingN >> 10; j += 3) {
            b *= j;

            if (b == 0) {
                b = 1;
            }

            a / b;
        }
    }
}

TEST(ArithmeticTiming, DivisionBoost) {
    uint512_t a = 1;

    for (size_t i = 1; i < DivisionTimingN; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }

        uint512_t b = 1;

        for (size_t j = 93; j < DivisionTimingN >> 10; j += 3) {
            b *= j;

            if (b == 0) {
                b = 1;
            }

            a / b;
        }
    }
}
