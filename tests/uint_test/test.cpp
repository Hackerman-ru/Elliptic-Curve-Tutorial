// clang-format off
#include "pch.h"
// clang-format on
#include "long-arithmetic.h"
#include "utils/csprng/csprng.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/fwd.hpp>

using boost::multiprecision::uint512_t;

#include <bitset>
#include <random>

using namespace elliptic_curve_guide;

static constexpr size_t int_conversion_correctness_n = 1000000;
static constexpr size_t string_correctness_n = 10000;
static constexpr size_t shift_correctness_n = 1000;
static constexpr size_t arithmetic_correctness_n = 10000;
static constexpr size_t division_correctness_n = 10000;

static constexpr size_t string_timing_n = 10000;
static constexpr size_t shift_timin_n = 10000;
static constexpr size_t arithmetic_timing_n = 10000;

static uint512_t generate_random_boost_uint() {
    duthomhas::csprng rng;
    std::seed_seq sseq {228};
    rng.seed(sseq);
    uint512_t result = 0;
    uint32_t x = rng(uint32_t());
    std::vector<uint32_t> values = rng(std::vector<uint32_t>(16));

    for (size_t i = 0; i < 16; ++i) {
        result <<= 32;
        result += values[i];
    }

    return result;
}

static uint_t<512> generate_random_my_uint() {
    duthomhas::csprng rng;
    std::seed_seq sseq {228};
    rng.seed(sseq);
    uint_t<512> result = 0;
    uint32_t x = rng(uint32_t());
    std::vector<uint32_t> values = rng(std::vector<uint32_t>(16));

    for (size_t i = 0; i < 16; ++i) {
        result <<= 32;
        result += values[i];
    }

    return result;
}

static std::string convert_to_binary(uint512_t value) {
    std::string result;

    do {
        std::bitset<64> bits = value.convert_to<uint64_t>();
        std::string bits_string = bits.to_string();
        std::reverse(bits_string.begin(), bits_string.end());
        result += bits_string;
        value >>= 64;
    } while (value > 0);

    result += "b0";
    std::reverse(result.begin(), result.end());
    return result;
}

template<typename From, typename To>
static To conv(const From& value) {
    To result = 0;

    for (size_t i = 16; i > 0; --i) {
        result <<= 32;
        result |= (value >> ((i - 1) * 32)).convert_to<uint32_t>();
    }

    return result;
}

static bool is_equal(uint512_t lhs, uint_t<512> rhs) {
    uint512_t rhs_value = conv<uint_t<512>, uint512_t>(rhs);
    return lhs == rhs_value;
}

#define comp(a, b)                                           \
    if (!is_equal(a, b)) {                                   \
        std::string boost_str = a.convert_to<std::string>(); \
        std::string my_str = b.convert_to<std::string>();    \
        ASSERT_EQ(boost_str, my_str);                        \
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
    ASSERT_EQ(boost_lhs, boost_rhs);
    ASSERT_EQ(my_lhs, my_rhs);
}

TEST(StringConversionCorrectness, DecimalStringConversion) {
    for (size_t i = 0; i < string_correctness_n; ++i) {
        uint512_t boost_value = generate_random_boost_uint();
        uint_t<512> my_value = conv<uint512_t, uint_t<512>>(boost_value);
        comp(boost_value, my_value);
    }
}

TEST(StringConversionCorrectness, HexadecimalStringConversion) {
    for (size_t i = 0; i < string_correctness_n; ++i) {
        uint512_t boost_value = generate_random_boost_uint();
        std::stringstream ss;
        ss << std::hex << std::showbase << boost_value;
        std::string hex_str = ss.str();
        uint_t<512> my_value(hex_str.c_str());
        comp(boost_value, my_value);
    }
}

TEST(StringConversionCorrectness, OctalStringConversion) {
    for (size_t i = 0; i < string_correctness_n; ++i) {
        uint512_t boost_value = generate_random_boost_uint();
        std::stringstream ss;
        ss << std::oct << std::showbase << boost_value;
        std::string oct_str = ss.str();
        uint_t<512> my_value(oct_str.c_str());
        comp(boost_value, my_value);
    }
}

TEST(StringConversionCorrectness, BinaryStringConversion) {
    for (size_t i = 0; i < string_correctness_n; ++i) {
        uint512_t boost_value = generate_random_boost_uint();
        std::string binary_str = convert_to_binary(boost_value);
        uint_t<512> my_value(binary_str.c_str());
        comp(boost_value, my_value);
    }
}

TEST(StringConversionTiming, DecimalStringConversionBoost) {
    for (size_t i = 0; i < string_timing_n; ++i) {
        uint512_t a = generate_random_boost_uint();
        std::string a_str = a.convert_to<std::string>();
    }
}

TEST(StringConversionTiming, DecimalStringConversion) {
    for (size_t i = 0; i < string_timing_n; ++i) {
        uint_t<512> a = generate_random_my_uint();
        std::string a_str = a.convert_to<std::string>();
    }
}

// Convert_to

TEST(IntConversionCorrectness, uint32_t) {
    std::mt19937 gen(42);

    for (size_t i = 0; i < int_conversion_correctness_n; ++i) {
        uint32_t value = gen();
        uint_t<512> my_value(value);
        ASSERT_EQ(value, my_value.convert_to<uint32_t>());
        uint512_t boost_value = value;
        comp(boost_value, my_value);
    }
}

TEST(IntConversionCorrectness, size_t) {
    std::mt19937_64 gen(42);

    for (size_t i = 0; i < int_conversion_correctness_n; ++i) {
        size_t value = gen();
        uint_t<512> my_value(value);
        ASSERT_EQ(value, my_value.convert_to<size_t>());
        uint512_t boost_value = value;
        comp(boost_value, my_value);
    }
}

// Shift

TEST(ShiftCorrectness, LeftShift) {
    for (size_t i = 0; i < shift_correctness_n; ++i) {
        uint512_t a = generate_random_boost_uint();
        uint_t<512> b = conv<uint512_t, uint_t<512>>(a);

        for (size_t j = 0; j < 512; ++j) {
            uint512_t boost_value = a << j;
            uint_t<512> my_value = b << j;
            comp(boost_value, my_value);

            if (boost_value == 0) {
                break;
            }
        }
    }
}

TEST(ShiftCorrectness, RightShift) {
    for (size_t i = 0; i < shift_correctness_n; ++i) {
        uint512_t a = generate_random_boost_uint();
        uint_t<512> b = conv<uint512_t, uint_t<512>>(a);

        for (size_t j = 0; j < 512; ++j) {
            uint512_t boost_value = a >> j;
            uint_t<512> my_value = b >> j;
            comp(boost_value, my_value);

            if (boost_value == 0) {
                break;
            }
        }
    }
}

TEST(ShiftTiming, LeftShift) {
    for (size_t i = 0; i < shift_timin_n; ++i) {
        uint_t<512> a = generate_random_my_uint();

        for (size_t shift = 0; shift < 512; ++shift) {
            auto b = a << shift;
        }
    }
}

TEST(ShiftTiming, LeftShiftBoost) {
    for (size_t i = 0; i < shift_timin_n; ++i) {
        uint512_t a = generate_random_boost_uint();

        for (size_t shift = 0; shift < 512; ++shift) {
            auto b = a << shift;
        }
    }
}

TEST(ShiftTiming, RightShift) {
    for (size_t i = 0; i < shift_timin_n; ++i) {
        uint_t<512> a = generate_random_my_uint();

        for (size_t shift = 0; shift < 512; ++shift) {
            auto b = a >> shift;
        }
    }
}

TEST(ShiftTiming, RightShiftBoost) {
    for (size_t i = 0; i < shift_timin_n; ++i) {
        uint512_t a = generate_random_boost_uint();

        for (size_t shift = 0; shift < 512; ++shift) {
            auto b = a >> shift;
        }
    }
}

// Arithmetic
TEST(ArithmeticCorrectness, Addition) {
    for (size_t i = 0; i < arithmetic_correctness_n; ++i) {
        uint_t<512> left = generate_random_my_uint();
        uint_t<512> right = generate_random_my_uint();
        uint_t<512> my_value = left + right;
        uint512_t boost_value = conv<uint_t<512>, uint512_t>(left) + conv<uint_t<512>, uint512_t>(right);
        comp(boost_value, my_value);
    }
}

TEST(ArithmeticCorrectness, Subtraction) {
    for (size_t i = 0; i < arithmetic_correctness_n; ++i) {
        uint_t<512> left = generate_random_my_uint();
        uint_t<512> right = generate_random_my_uint();
        uint_t<512> my_value = left - right;
        uint512_t boost_value = conv<uint_t<512>, uint512_t>(left) - conv<uint_t<512>, uint512_t>(right);
        comp(boost_value, my_value);
    }
}

TEST(ArithmeticCorrectness, Multiplication) {
    for (size_t i = 0; i < arithmetic_correctness_n; ++i) {
        uint_t<512> left = generate_random_my_uint();
        uint_t<512> right = generate_random_my_uint();
        uint_t<512> my_value = left * right;
        uint512_t boost_value = conv<uint_t<512>, uint512_t>(left) * conv<uint_t<512>, uint512_t>(right);
        comp(boost_value, my_value);
    }
}

TEST(ArithmeticCorrectness, Division) {
    for (size_t i = 0; i < arithmetic_correctness_n; ++i) {
        uint_t<512> left = generate_random_my_uint();
        uint_t<512> right = generate_random_my_uint();
        uint_t<512> my_value = left / right;
        uint512_t boost_value = conv<uint_t<512>, uint512_t>(left) / conv<uint_t<512>, uint512_t>(right);
        comp(boost_value, my_value);
    }
}

TEST(ArithmeticTiming, Addition) {
    for (size_t i = 0; i < arithmetic_timing_n; ++i) {
        uint_t<512> a = generate_random_my_uint() + generate_random_my_uint();
    }
}

TEST(ArithmeticTiming, AdditionBoost) {
    for (size_t i = 0; i < arithmetic_timing_n; ++i) {
        uint512_t a = generate_random_boost_uint() + generate_random_boost_uint();
    }
}

TEST(ArithmeticTiming, Subtraction) {
    for (size_t i = 0; i < arithmetic_timing_n; ++i) {
        uint_t<512> a = generate_random_my_uint() - generate_random_my_uint();
    }
}

TEST(ArithmeticTiming, SubtractionBoost) {
    for (size_t i = 0; i < arithmetic_timing_n; ++i) {
        uint512_t a = generate_random_boost_uint() - generate_random_boost_uint();
    }
}

TEST(ArithmeticTiming, Multiplication) {
    for (size_t i = 0; i < arithmetic_timing_n; ++i) {
        uint_t<512> a = generate_random_my_uint() * generate_random_my_uint();
    }
}

TEST(ArithmeticTiming, MultiplicationBoost) {
    for (size_t i = 0; i < arithmetic_timing_n; ++i) {
        uint512_t a = generate_random_boost_uint() * generate_random_boost_uint();
    }
}

TEST(ArithmeticTiming, Division) {
    for (size_t i = 1; i < arithmetic_timing_n; ++i) {
        uint_t<512> a = generate_random_my_uint() / generate_random_my_uint();
    }
}

TEST(ArithmeticTiming, DivisionBoost) {
    for (size_t i = 1; i < arithmetic_timing_n; ++i) {
        uint512_t a = generate_random_boost_uint() / generate_random_boost_uint();
    }
}
