#define BOOST_TEST_MODULE uint_test
#include "LongInt/long-int.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/test/included/unit_test.hpp>
#include <format>

using namespace boost::multiprecision;

static constexpr size_t N = 100000;

BOOST_AUTO_TEST_CASE(OneStringConvertion) {
    std::string str = "6227020800";
    uint_t<512> my(str);
    std::string my_str = my.into_string();
    BOOST_TEST(str == my_str);
}

BOOST_AUTO_TEST_CASE(StringConvertion) {
    uint512_t boost(1);
    for (size_t i = 1; i < N; ++i) {
        boost *= uint512_t(i);
        std::string boost_str = boost.convert_to<std::string>();
        std::string my_str = (uint_t<512>(boost_str)).into_string();
        BOOST_TEST_REQUIRE(boost_str == my_str);
    }
}

BOOST_AUTO_TEST_CASE(Addition) {
    uint512_t boost(1);
    uint_t<512> my(1);
    for (size_t i = 1; i < N; ++i) {
        boost += uint512_t(i);
        my += uint_t<512>(i);
        std::string boost_str = boost.convert_to<std::string>();
        std::string my_str = my.into_string();
        BOOST_TEST_REQUIRE(boost_str == my_str);
    }
}

BOOST_AUTO_TEST_CASE(Subtraction) {
    uint512_t boost(1);
    uint_t<512> my(1);
    for (size_t i = 1; i < N; ++i) {
        boost -= uint512_t(i);
        my -= uint_t<512>(i);
        std::string boost_str = boost.convert_to<std::string>();
        std::string my_str = my.into_string();
        BOOST_TEST_REQUIRE(boost_str == my_str);
    }
}

BOOST_AUTO_TEST_CASE(Multiplication) {
    uint512_t boost(1);
    uint_t<512> my(1);
    for (size_t i = 1; i < N; ++i) {
        boost *= uint512_t(i);
        my *= uint_t<512>(i);
        std::string boost_str = boost.convert_to<std::string>();
        std::string my_str = my.into_string();
        BOOST_TEST_REQUIRE(boost_str == my_str);
    }
}

void divison_check(const std::string& dividend, const std::string& divisor) {
    uint512_t boost_lhs(dividend);
    uint512_t boost_rhs(divisor);
    uint_t<512> my_lhs(dividend);
    uint_t<512> my_rhs(divisor);

    uint512_t boost_result = boost_lhs / boost_rhs;
    uint_t<512> my_result = my_lhs / my_rhs;

    std::string boost_str = boost_result.convert_to<std::string>();
    std::string my_str = my_result.into_string();
    std::string msg = std::format("\n{} divided by {} equals {},\nbut uint {} divided by {} equals {}\n",
                                  boost_lhs.convert_to<std::string>(),
                                  boost_rhs.convert_to<std::string>(),
                                  boost_result.convert_to<std::string>(),
                                  my_lhs.into_string(),
                                  my_rhs.into_string(),
                                  my_result.into_string());

    BOOST_TEST_REQUIRE(boost_str == my_str, msg);
}

BOOST_AUTO_TEST_CASE(OneDivision) {
    std::string dividend = "1307674368000";
    std::string divisor = "1307674368000";
    divison_check(dividend, divisor);
    dividend = "51090942171709440000";
    divisor = "37822859361";
    divison_check(dividend, divisor);
}

void check_division(const uint512_t& boost_lhs,
                    const uint512_t& boost_rhs,
                    const uint_t<512>& my_lhs,
                    const uint_t<512>& my_rhs) {
    BOOST_TEST_REQUIRE(boost_lhs.convert_to<std::string>() == my_lhs.into_string());
    BOOST_TEST_REQUIRE(boost_rhs.convert_to<std::string>() == my_rhs.into_string());

    uint512_t boost_result = boost_lhs / boost_rhs;
    uint_t<512> my_result = my_lhs / my_rhs;
    std::string boost_str = boost_result.convert_to<std::string>();
    std::string my_str = my_result.into_string();
    std::string msg = std::format("\n{} mod {} equals {},\nbut uint {} mod {} equals {}\n",
                                  boost_lhs.convert_to<std::string>(),
                                  boost_rhs.convert_to<std::string>(),
                                  boost_result.convert_to<std::string>(),
                                  my_lhs.into_string(),
                                  my_rhs.into_string(),
                                  my_result.into_string());

    BOOST_TEST_REQUIRE(boost_str == my_str, msg);
}

void check_mod(const uint512_t& boost_lhs,
               const uint512_t& boost_rhs,
               const uint_t<512>& my_lhs,
               const uint_t<512>& my_rhs) {
    BOOST_TEST_REQUIRE(boost_lhs.convert_to<std::string>() == my_lhs.into_string());
    BOOST_TEST_REQUIRE(boost_rhs.convert_to<std::string>() == my_rhs.into_string());

    uint512_t boost_result = boost_lhs % boost_rhs;
    uint_t<512> my_result = my_lhs % my_rhs;
    std::string boost_str = boost_result.convert_to<std::string>();
    std::string my_str = my_result.into_string();
    std::string msg = std::format("\n{} mod {} equals {},\nbut uint {} mod {} equals {}\n",
                                  boost_lhs.convert_to<std::string>(),
                                  boost_rhs.convert_to<std::string>(),
                                  boost_result.convert_to<std::string>(),
                                  my_lhs.into_string(),
                                  my_rhs.into_string(),
                                  my_result.into_string());

    BOOST_TEST_REQUIRE(boost_str == my_str, msg);
}

BOOST_AUTO_TEST_CASE(Division) {
    uint512_t boost_lhs = 1;
    uint_t<512> my_lhs = uint_t<512>(1);
    for (size_t i = 1; i < N >> 7; i += 3) {
        boost_lhs *= uint512_t(i);
        my_lhs *= uint_t<512>(i);

        if (boost_lhs == 0) {
            BOOST_TEST_REQUIRE(boost_lhs.convert_to<std::string>() == my_lhs.into_string());
            boost_lhs = 1;
            my_lhs = uint_t<512>(1);
        }

        uint512_t boost_rhs(1);
        uint_t<512> my_rhs(1);

        for (size_t j = i >> 2; j <= i; j += 2) {
            boost_rhs *= uint512_t(j);
            my_rhs *= uint_t<512>(j);

            if (boost_rhs == 0) {
                BOOST_TEST_REQUIRE(boost_rhs.convert_to<std::string>() == my_rhs.into_string());
                boost_rhs = 1;
                my_rhs = uint_t<512>(1);
            }

            check_division(boost_lhs, boost_rhs, my_lhs, my_rhs);
            check_mod(boost_lhs, boost_rhs, my_lhs, my_rhs);
        }
    }
}
