#include "Field/Field.h"
#include "Uint/uint.h"

#include <iostream>

int main() {
    return 0;
}

// Tests that are in another project
//#include "../Elliptic-Curve-Tutorial/Uint/uint.h"
//#include "../Elliptic-Curve-Tutorial/util.h"
//#include "pch.h"
//
//#include <boost/multiprecision/cpp_int.hpp>
//
//using boost::multiprecision::uint512_t;
//using namespace ECG;
//static constexpr size_t N = 100000;
//
//static constexpr size_t StringCorrectnessN = 10000;
//static constexpr size_t StringTimingN = 10000;
//
//TEST(StringConversionCorrectness, DecimalStringConversion) {
//    uint512_t a = 1;
//
//    for (size_t i = 1; i < StringCorrectnessN; ++i) {
//        a *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//
//        auto a_str = a.convert_to<std::string>();
//        uint_t<512> b(a_str);
//        auto b_str = b.into_string();
//        EXPECT_EQ(a_str, b_str);
//    }
//}
//
//TEST(StringConversionCorrectness, BinaryStringConversion) {
//    uint512_t a = 1;
//
//    for (size_t i = 1; i < StringCorrectnessN; ++i) {
//        a *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//
//        uint512_t a_clone = a;
//        std::string a_str;
//
//        while (a_clone != 0) {
//            a_str.push_back(((a_clone & 1) != 0) + '0');
//            a_clone >>= 1;
//        }
//
//        std::reverse(a_str.begin(), a_str.end());
//
//        uint_t<512> b(a_str, StringType::BINARY);
//        auto b_str = b.into_string(StringType::BINARY);
//        EXPECT_EQ(a_str, b_str);
//    }
//}
//
//TEST(StringConversionCorrectness, HexadecimalStringConversion) {
//    uint512_t a = 1;
//
//    for (size_t i = 1; i < StringCorrectnessN; ++i) {
//        a *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//
//        uint512_t a_clone = a;
//        std::string a_str;
//
//        while (a_clone != 0) {
//            auto n = a_clone.convert_to<uint32_t>() & 0xF;
//
//            if (n < 10) {
//                a_str.push_back(n + '0');
//            } else {
//                n -= 10;
//                a_str.push_back(n + 'a');
//            }
//
//            a_clone >>= 4;
//        }
//
//        std::reverse(a_str.begin(), a_str.end());
//
//        uint_t<512> b(a_str, StringType::HEXADECIMAL);
//        auto b_str = b.into_string(StringType::HEXADECIMAL);
//        EXPECT_EQ(a_str, b_str);
//    }
//}
//
//TEST(StringConversionTiming, DecimalStringConversionBoost) {
//    uint512_t a = 1;
//
//    for (size_t i = 1; i < StringTimingN; ++i) {
//        a *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//
//        std::string a_str = a.convert_to<std::string>();
//    }
//}
//
//TEST(StringConversionTiming, DecimalStringConversion) {
//    uint_t<512> a = 1;
//
//    for (size_t i = 1; i < StringTimingN; ++i) {
//        a *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//
//        std::string a_str = a.into_string();
//    }
//}
//
//TEST(StringConversionTiming, BinaryStringConversionBoost) {
//    uint512_t a = 1;
//
//    for (size_t i = 1; i < StringCorrectnessN; ++i) {
//        a *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//
//        std::string a_str;
//
//        while (a != 0) {
//            a_str += ((a & 1) != 0) + '0';
//            a >>= 1;
//        }
//
//        std::reverse(a_str.begin(), a_str.end());
//    }
//}
//
//TEST(StringConversionTiming, BinaryStringConversion) {
//    uint_t<512> a = 1;
//
//    for (size_t i = 1; i < StringCorrectnessN; ++i) {
//        a *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//
//        std::string a_str = a.into_string(StringType::BINARY);
//    }
//}
//
//TEST(StringConversionTiming, HexadecimalStringConversionBoost) {
//    uint512_t a = 1;
//
//    for (size_t i = 1; i < StringCorrectnessN; ++i) {
//        a *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//
//        std::string a_str;
//
//        while (a != 0) {
//            auto n = a.convert_to<uint32_t>() & 0xF;
//
//            if (n < 10) {
//                a_str.push_back(n + '0');
//            } else {
//                n -= 10;
//                a_str.push_back(n + 'a');
//            }
//
//            a >>= 4;
//        }
//
//        std::reverse(a_str.begin(), a_str.end());
//    }
//}
//
//TEST(StringConversionTiming, HexadecimalStringConversion) {
//    uint_t<512> a = 1;
//
//    for (size_t i = 1; i < StringCorrectnessN; ++i) {
//        a *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//
//        std::string a_str = a.into_string(StringType::HEXADECIMAL);
//    }
//}
//
//TEST(ArithmeticCorrectness, Addition) {
//    uint512_t a = 1;
//    uint_t<512> b = 1;
//
//    for (size_t i = 0; i < N; ++i) {
//        a += i;
//        b += i;
//
//        std::string a_str = a.convert_to<std::string>();
//        std::string b_str = b.into_string();
//
//        EXPECT_EQ(a_str, b_str);
//    }
//}
//
//TEST(ArithmeticCorrectness, Subtraction) {
//    uint512_t a = 1;
//    uint_t<512> b = 1;
//
//    for (size_t i = 0; i < N; ++i) {
//        a -= i;
//        b -= i;
//
//        std::string a_str = a.convert_to<std::string>();
//        std::string b_str = b.into_string();
//
//        EXPECT_EQ(a_str, b_str);
//    }
//}
//
//TEST(ArithmeticCorrectness, Multiplication) {
//    uint512_t a = 1;
//    uint_t<512> b = 1;
//
//    for (size_t i = 1; i < N; ++i) {
//        a *= i;
//        b *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//
//        if (b == 0) {
//            b = 1;
//        }
//
//        std::string a_str = a.convert_to<std::string>();
//        std::string b_str = b.into_string();
//
//        EXPECT_EQ(a_str, b_str);
//    }
//}
//
//TEST(ArithmeticCorrectness, Division) {
//    uint512_t a = 1;
//    uint_t<512> b = 1;
//
//    for (size_t i = 1; i < N; ++i) {
//        a *= i;
//        b *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//
//        if (b == 0) {
//            b = 1;
//        }
//
//        uint512_t c = 1;
//        uint_t<512> d = 1;
//
//        for (size_t j = 93; j < N >> 10; j += 3) {
//            c *= j;
//            d *= j;
//
//            if (c == 0) {
//                c = 1;
//            }
//
//            if (d == 0) {
//                d = 1;
//            }
//
//            std::string boost_str = (a / c).convert_to<std::string>();
//            std::string my_str = (b / d).into_string();
//
//            EXPECT_EQ(boost_str, my_str);
//        }
//    }
//}
//
//static constexpr size_t ArithmeticTimingN = 1000000;
//
//TEST(ArithmeticTiming, Addition) {
//    uint_t<512> b = 1;
//
//    for (size_t i = 0; i < ArithmeticTimingN; ++i) {
//        b += i;
//    }
//}
//
//TEST(ArithmeticTiming, AdditionBoost) {
//    uint512_t a = 1;
//    uint_t<512> b = 1;
//
//    for (size_t i = 0; i < ArithmeticTimingN; ++i) {
//        a += i;
//    }
//}
//
//TEST(ArithmeticTiming, Subtraction) {
//    uint_t<512> b = 1;
//
//    for (size_t i = 0; i < ArithmeticTimingN; ++i) {
//        b -= i;
//    }
//}
//
//TEST(ArithmeticTiming, SubtractionBoost) {
//    uint512_t a = 1;
//    uint_t<512> b = 1;
//
//    for (size_t i = 0; i < ArithmeticTimingN; ++i) {
//        a -= i;
//    }
//}
//
//TEST(ArithmeticTiming, Multiplication) {
//    uint_t<512> a = 1;
//
//    for (size_t i = 1; i < ArithmeticTimingN; ++i) {
//        a *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//    }
//}
//
//TEST(ArithmeticTiming, MultiplicationBoost) {
//    uint512_t a = 1;
//
//    for (size_t i = 1; i < ArithmeticTimingN; ++i) {
//        a *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//    }
//}
//
//static constexpr size_t DivisionTimingN = 100000;
//
//TEST(ArithmeticTiming, Division) {
//    uint_t<512> a = 1;
//
//    for (size_t i = 1; i < DivisionTimingN; ++i) {
//        a *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//
//        uint_t<512> b = 1;
//
//        for (size_t j = 93; j < DivisionTimingN >> 10; j += 3) {
//            b *= j;
//
//            if (b == 0) {
//                b = 1;
//            }
//
//            a / b;
//        }
//    }
//}
//
//TEST(ArithmeticTiming, DivisionBoost) {
//    uint512_t a = 1;
//
//    for (size_t i = 1; i < DivisionTimingN; ++i) {
//        a *= i;
//
//        if (a == 0) {
//            a = 1;
//        }
//
//        uint512_t b = 1;
//
//        for (size_t j = 93; j < DivisionTimingN >> 10; j += 3) {
//            b *= j;
//
//            if (b == 0) {
//                b = 1;
//            }
//
//            a / b;
//        }
//    }
//}
