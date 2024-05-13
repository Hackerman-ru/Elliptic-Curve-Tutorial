// clang-format off
#include "pch.h"
// clang-format on
#include "long-arithmetic.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/fwd.hpp>

using boost::multiprecision::uint512_t;

#include <random>

using namespace elliptic_curve_guide;

static constexpr size_t convert_correctness_n = 10000000;
static constexpr size_t string_correctness_n = 10000;
static constexpr size_t shift_correctness_n = 1000;
static constexpr size_t arithmetic_correctness_n = 100000;
static constexpr size_t division_correctness_n = 10000;

static constexpr size_t string_timing_n = 10000;
static constexpr size_t shift_timin_n = 10000000;
static constexpr size_t arithmetic_timing_n = 1000000;
static constexpr size_t division_timing_n = 100000;

template<typename From, typename To>
static To conv(const From& value) {
    To result = 0;

    for (size_t i = 16; i > 0; --i) {
        result <<= 32;
        result |= (value >> ((i - 1) * 32)).convert_to<uint32_t>();
    }

    return result;
}

static void comp(uint512_t boost_value, uint_t<512> my_value_) {
    uint512_t my_value = conv<uint_t<512>, uint512_t>(my_value_);
    ASSERT_EQ(my_value, boost_value);
}

TEST(StringConversionCorrectness, SimpleConversions) {
    uint_t<512> my_p("0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff");
    uint512_t boost_p("0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff");
    comp(boost_p, my_p);
    uint_t<512> my_a = "0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc";
    uint512_t boost_a("0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc");
    comp(boost_a, my_a);
    uint_t<512> my_b = "0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b";
    uint512_t boost_b("0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b");
    comp(boost_b, my_b);
    uint_t<512> my_x = "0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296";
    uint512_t boost_x("0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296");
    comp(boost_x, my_x);
    uint_t<512> my_y = "0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5";
    uint512_t boost_y("0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5");
    comp(boost_y, my_y);
    uint_t<512> my_n = "0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551";
    uint512_t boost_n("0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551");
    comp(boost_n, my_n);

    uint_t<512> my_lhs = (my_y * my_y) % my_p;
    uint512_t boost_lhs = (boost_y * boost_y) % boost_p;
    comp(boost_lhs, my_lhs);

    uint_t<512> my_rhs = (((my_x * my_x) % my_p) * my_x) % my_p + (my_a * my_x) % my_p + my_b - my_p;
    uint512_t boost_rhs = (((boost_x * boost_x) % boost_p) * boost_x) % boost_p
                        + (boost_a * boost_x) % boost_p + boost_b - boost_p;
    comp(boost_rhs, my_rhs);
    ASSERT_EQ(my_lhs, my_rhs);
    ASSERT_EQ(boost_lhs, boost_rhs);
}

TEST(StringConversionCorrectness, DecimalStringConversion) {
    uint512_t boost_value = 1;

    for (size_t i = 1; i < string_correctness_n; ++i) {
        boost_value *= i;

        if (boost_value == 0) {
            boost_value = 1;
        }

        auto boost_str = boost_value.convert_to<std::string>();
        uint_t<512> my_value(boost_str.c_str());
        auto my_str = my_value.convert_to<std::string>();
        ASSERT_EQ(boost_str, my_str);
    }
}

TEST(StringConversionTiming, DecimalStringConversionBoost) {
    uint512_t a = 1;

    for (size_t i = 1; i < string_timing_n; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }

        std::string a_str = a.convert_to<std::string>();
    }
}

TEST(StringConversionTiming, DecimalStringConversion) {
    uint_t<512> a = 1;

    for (size_t i = 1; i < string_timing_n; ++i) {
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

    for (size_t i = 0; i < convert_correctness_n; ++i) {
        uint32_t a = gen();
        uint_t<512> b(a);
        ASSERT_EQ(a, b.convert_to<uint32_t>());
    }
}

TEST(IntConversionCorrectness, size_t) {
    std::mt19937_64 gen;

    for (size_t i = 0; i < convert_correctness_n; ++i) {
        size_t a = gen();
        uint_t<512> b(a);
        ASSERT_EQ(a, b.convert_to<size_t>());
    }
}

// Shift

TEST(ShiftCorrectness, LeftShift) {
    uint512_t a("9999999999999999999999999999999999");

    for (size_t i = 1; i < shift_correctness_n; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }

        uint_t<512> b = conv<uint512_t, uint_t<512>>(a);

        for (size_t j = 0; j < 512; ++j) {
            uint512_t boost_value = a << j;
            uint512_t my_value = conv<uint_t<512>, uint512_t>(b << j);
            ASSERT_EQ(my_value, boost_value);

            if (a << j == 0) {
                break;
            }
        }
    }
}

TEST(ShiftCorrectness, RightShift) {
    uint512_t a("9999999999999999999999999999999999");

    for (size_t i = 1; i < shift_correctness_n; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }

        uint_t<512> b = conv<uint512_t, uint_t<512>>(a);

        for (size_t j = 0; j < 512; ++j) {
            uint512_t boost_value = a >> j;
            uint512_t my_value = conv<uint_t<512>, uint512_t>(b >> j);
            ASSERT_EQ(my_value, boost_value);

            if (a >> j == 0) {
                break;
            }
        }
    }
}

TEST(ShiftTiming, LeftShift) {
    uint_t<512> a("9999999999999999999999999999999999");

    for (size_t i = 0; i < shift_timin_n; ++i) {
        auto b = a << i;
    }
}

TEST(ShiftTiming, LeftShiftBoost) {
    uint512_t a("9999999999999999999999999999999999");

    for (size_t i = 0; i < shift_timin_n; ++i) {
        auto b = a << i;
    }
}

TEST(ShiftTiming, RightShift) {
    uint_t<512> a("9999999999999999999999999999999999");

    for (size_t i = 0; i < shift_timin_n; ++i) {
        auto b = a >> i;
    }
}

TEST(ShiftTiming, RightShiftBoost) {
    uint512_t a("9999999999999999999999999999999999");

    for (size_t i = 0; i < shift_timin_n; ++i) {
        auto b = a >> i;
    }
}

// Arithmetic
TEST(ArithmeticCorrectness, Addition) {
    uint512_t a = 1;
    uint_t<512> b = 1;

    for (size_t i = 0; i < arithmetic_correctness_n; ++i) {
        a += i;
        b += i;

        uint512_t boost_value = a;
        uint512_t my_value = conv<uint_t<512>, uint512_t>(b);
        ASSERT_EQ(my_value, boost_value);
    }
}

TEST(ArithmeticCorrectness, Subtraction) {
    uint512_t a = 1;
    uint_t<512> b = 1;

    for (size_t i = 0; i < arithmetic_correctness_n; ++i) {
        a -= i;
        b -= i;

        uint512_t boost_value = a;
        uint512_t my_value = conv<uint_t<512>, uint512_t>(b);
        ASSERT_EQ(my_value, boost_value);
    }
}

TEST(ArithmeticCorrectness, Multiplication) {
    uint512_t a = 1;
    uint_t<512> b = 1;

    for (size_t i = 1; i < arithmetic_correctness_n; ++i) {
        a *= i;
        b *= i;

        if (a == 0 || b == 0) {
            ASSERT_EQ(a, 0);
            ASSERT_EQ(b, 0);
            a = 1;
            b = 1;
        }

        uint512_t boost_value = a;
        uint512_t my_value = conv<uint_t<512>, uint512_t>(b);
        ASSERT_EQ(my_value, boost_value);
    }
}

TEST(ArithmeticCorrectness, Division) {
    uint512_t a = 1;
    uint_t<512> b = 1;

    for (size_t i = 1; i < division_correctness_n; ++i) {
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

        for (size_t j = 0; j < division_correctness_n >> 8; j += 3) {
            c *= j;
            d *= j;

            if (c == 0) {
                c = 1;
            }

            if (d == 0) {
                d = 1;
            }

            uint512_t boost_value = a / c;
            uint512_t my_value = conv<uint_t<512>, uint512_t>(b / d);
            ASSERT_EQ(my_value, boost_value);
        }
    }
}

TEST(ArithmeticTiming, Addition) {
    uint_t<512> b = 1;

    for (size_t i = 0; i < arithmetic_timing_n; ++i) {
        b += i;
    }
}

TEST(ArithmeticTiming, AdditionBoost) {
    uint512_t a = 1;
    uint_t<512> b = 1;

    for (size_t i = 0; i < arithmetic_timing_n; ++i) {
        a += i;
    }
}

TEST(ArithmeticTiming, Subtraction) {
    uint_t<512> b = 1;

    for (size_t i = 0; i < arithmetic_timing_n; ++i) {
        b -= i;
    }
}

TEST(ArithmeticTiming, SubtractionBoost) {
    uint512_t a = 1;
    uint_t<512> b = 1;

    for (size_t i = 0; i < arithmetic_timing_n; ++i) {
        a -= i;
    }
}

TEST(ArithmeticTiming, Multiplication) {
    uint_t<512> a = 1;

    for (size_t i = 1; i < arithmetic_timing_n; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }
    }
}

TEST(ArithmeticTiming, MultiplicationBoost) {
    uint512_t a = 1;

    for (size_t i = 1; i < arithmetic_timing_n; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }
    }
}

TEST(ArithmeticTiming, Division) {
    uint_t<512> a = 1;

    for (size_t i = 1; i < division_timing_n; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }

        uint_t<512> b = 1;

        for (size_t j = 93; j < division_timing_n >> 10; j += 3) {
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

    for (size_t i = 1; i < division_timing_n; ++i) {
        a *= i;

        if (a == 0) {
            a = 1;
        }

        uint512_t b = 1;

        for (size_t j = 93; j < division_timing_n >> 10; j += 3) {
            b *= j;

            if (b == 0) {
                b = 1;
            }

            a / b;
        }
    }
}
