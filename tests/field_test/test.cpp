// clang-format off
#include "pch.h"
// clang-format on
#include "field.h"
#include "utils/primes.h"
#include "utils/random.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/fwd.hpp>

using namespace elliptic_curve_guide;
using namespace field;

using namespace algorithm::random;

static constexpr size_t c_primes_n = 150;
static constexpr size_t c_correctness_test_arithmetic_n = 50;
static constexpr size_t c_correctness_test_inversion_n = 1600000;
static constexpr size_t c_correctness_test_inversion_prime = 1599827;
static constexpr size_t c_correctness_test_power_n = 100;
static constexpr size_t c_correctness_test_shift_n = 100;

static constexpr size_t c_timing_test_arithmetic_n = 200;
static constexpr size_t c_timing_test_power_n = 200;
static constexpr size_t c_timing_test_shift_n = 200;

#define UCMP(my_value, correct_value)                                      \
    if (my_value != correct_value) {                                       \
        std::string my_str = my_value.convert_to<std::string>();           \
        std::string correct_str = correct_value.convert_to<std::string>(); \
        ASSERT_EQ(my_str, correct_str);                                    \
    }

#define FCMP(my_value, correct_value)                                              \
    if (my_value != correct_value) {                                               \
        std::string my_str = my_value.value().convert_to<std::string>();           \
        std::string correct_str = correct_value.value().convert_to<std::string>(); \
        ASSERT_EQ(my_str, correct_str);                                            \
    }

TEST(SimpleTest, Creating) {
    Field f("7");
    FieldElement a = f.element(10);
    ASSERT_EQ(a.value(), 3);
    a = f.element(99);
    ASSERT_EQ(a.value(), 1);
}

TEST(SimpleTest, Multiplication) {
    Field f("7");
    FieldElement a = f.element(3);
    FieldElement b = a;
    uint c = (a * b).value();
    ASSERT_EQ(c, 2);
}

TEST(SimpleTest, Negotiation) {
    Field f("7");
    FieldElement a = -f.element(3);
    uint c = FieldElement::pow(a, 10).value();
    ASSERT_EQ(c, 4);
}

TEST(SimpleTest, Inversion) {
    Field f("1000000007");
    FieldElement a = f.element("999999999");
    FieldElement b = f.element(2);
    FieldElement inverse = a * b;
    inverse.inverse();
    uint result = inverse.value();
    uint correct_result("437500003");
    UCMP(result, correct_result);
}

TEST(SimpleTest, Comparison) {
    Field f("1000000007");
    FieldElement a = f.element("999999999");
    FieldElement b = f.element(2);
    ASSERT_EQ(a < b, false);
}

TEST(SimpleTest, Shift) {
    Field f("1000000007");
    FieldElement a = f.element("999999999");
    FieldElement b = a << uint("30");
    uint result = b.value();
    uint correct_result("410065471");
    UCMP(result, correct_result);
}

TEST(CorrectnessTest, Creating) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);
        UCMP(f.modulus(), p);
    }
}

TEST(CorrectnessTest, Comparison) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_correctness_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement b = generate_random_field_element(f);
            uint a_value = a.value();
            uint b_value = b.value();
            ASSERT_EQ(a != b, a_value != b_value);
            ASSERT_EQ(a == b, a_value == b_value);
            ASSERT_EQ(a < b, a_value < b_value);
            ASSERT_EQ(a > b, a_value > b_value);
        }
    }
}

TEST(CorrectnessTest, Addition) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_correctness_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement b = generate_random_field_element(f);
            uint my_value = (a + b).value();
            uint correct_value = (a.value() + b.value()) % p;
            UCMP(my_value, correct_value);
        }
    }
}

TEST(CorrectnessTest, Negotiation) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_correctness_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement b = -a;
            uint my_value = (a + b).value();
            uint correct_value = 0;
            UCMP(my_value, correct_value);
        }
    }
}

TEST(CorrectnessTest, Subtraction) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_correctness_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement b = generate_random_field_element(f);
            uint my_value = (a - b).value();
            uint correct_value = (a.value() + p - b.value()) % p;
            UCMP(my_value, correct_value);
        }
    }
}

TEST(CorrectnessTest, Multiplication) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_correctness_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement b = generate_random_field_element(f);
            uint my_value = (a * b).value();
            uint correct_value = (a.value() * b.value()) % p;
            UCMP(my_value, correct_value);
        }
    }
}

TEST(CorrectnessTest, Division) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_correctness_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement b = generate_random_non_zero_field_element(f);
            FieldElement q = a / b;
            uint my_value = (q * b).value();
            uint correct_value = a.value();
            UCMP(my_value, correct_value);
        }
    }
}

TEST(CorrectnessTest, Inversion) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);
        FieldElement one = f.element(1);

        for (size_t j = 0; j < c_correctness_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_non_zero_field_element(f);
            FieldElement inv_a = FieldElement::inverse(a);
            FieldElement result = a * inv_a;
            FCMP(result, one);
        }
    }
}

TEST(CorrectnessTest, Power) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_correctness_test_power_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            uint my_value = FieldElement::pow(a, j).value();

            uint correct_value = 1;
            size_t current_power = 0;

            while (current_power < j) {
                correct_value *= a.value();
                correct_value %= p;
                current_power++;
            }

            UCMP(my_value, correct_value);
        }
    }
}

TEST(CorrectnessTest, Shift) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_correctness_test_shift_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            uint my_value = (a << j).value();

            uint correct_value = a.value();
            size_t current_shift = 0;

            while (current_shift < j) {
                correct_value <<= 1;
                correct_value %= p;
                current_shift++;
            }

            UCMP(my_value, correct_value);
        }
    }
}

// Timing measurements
TEST(TimingTest, Addition) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_timing_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement b = generate_random_field_element(f);
            FieldElement result = a + b;
        }
    }
}

TEST(TimingTest, Negotiation) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_timing_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement result = -a;
        }
    }
}

TEST(TimingTest, Subtraction) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_timing_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement b = generate_random_field_element(f);
            FieldElement result = a - b;
        }
    }
}

TEST(TimingTest, Multiplication) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_timing_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement b = generate_random_field_element(f);
            FieldElement result = a * b;
        }
    }
}

TEST(TimingTest, Division) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_timing_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement b = generate_random_non_zero_field_element(f);
            FieldElement result = a / b;
        }
    }
}

TEST(TimingTest, Inversion) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_timing_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_non_zero_field_element(f);
            FieldElement result = FieldElement::inverse(a);
        }
    }
}

TEST(TimingTest, Power) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_timing_test_power_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            uint power = generate_random_uint_modulo(p);
            FieldElement result = FieldElement::pow(a, power);
        }
    }
}

TEST(TimingTest, Shift) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = primes::prime_number_list[i];
        Field f(p);

        for (size_t j = 0; j < c_timing_test_shift_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement result = a << j;
        }
    }
}
