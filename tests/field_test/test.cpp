// clang-format off
#include "pch.h"
// clang-format on
#include "field.h"
#include "utils/primes.h"
#include "utils/random.h"

#include <random>

using namespace elliptic_curve_guide;
using namespace field;
using namespace algorithm::random;

static constexpr size_t c_primes_n = 150;
static constexpr size_t c_correctness_test_n = 100;

static constexpr size_t c_stress_test_arithmetic_n = 200;
static constexpr size_t c_stress_test_inversion_n = 400;
static constexpr size_t c_stress_test_negotiation_n = 400;
static constexpr size_t c_stress_test_power_n = 200;
static constexpr size_t c_stress_test_shift_n = 200;

#define UINT_EQ(my_value, correct_value)                                   \
    if (my_value != correct_value) {                                       \
        std::string my_str = my_value.convert_to<std::string>();           \
        std::string correct_str = correct_value.convert_to<std::string>(); \
        ASSERT_EQ(my_str, correct_str);                                    \
    }

#define FIELD_EQ(my_value, correct_value)                                          \
    if (my_value != correct_value) {                                               \
        std::string my_str = my_value.value().convert_to<std::string>();           \
        std::string correct_str = correct_value.value().convert_to<std::string>(); \
        ASSERT_EQ(my_str, correct_str);                                            \
    }

std::mt19937_64 gen(42);

static uint get_random_prime() {
    size_t pos = gen() % (c_primes_n - 2) + 2;
    return primes::prime_number_list[pos];
}

// Simple tests
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
    UINT_EQ(result, correct_result);
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
    UINT_EQ(result, correct_result);
}

// Correctness tests

TEST(CorrectnessTest, Creating) {
    for (size_t i = 0; i < c_correctness_test_n; ++i) {
        uint p = get_random_prime();
        Field f(p);
        UINT_EQ(f.modulus(), p);
    }
}

TEST(CorrectnessTest, Comparison) {
    for (size_t i = 0; i < c_correctness_test_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

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

TEST(CorrectnessTest, Addition) {
    for (size_t i = 0; i < c_correctness_test_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        FieldElement a = generate_random_field_element(f);
        FieldElement b = generate_random_field_element(f);
        uint my_value = (a + b).value();
        uint correct_value = (a.value() + b.value()) % p;
        UINT_EQ(my_value, correct_value);
    }
}

TEST(CorrectnessTest, Negotiation) {
    for (size_t i = 0; i < c_correctness_test_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        FieldElement a = generate_random_field_element(f);
        FieldElement b = -a;
        uint my_value = (a + b).value();
        uint correct_value = 0;
        UINT_EQ(my_value, correct_value);
    }
}

TEST(CorrectnessTest, Subtraction) {
    for (size_t i = 0; i < c_correctness_test_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        FieldElement a = generate_random_field_element(f);
        FieldElement b = generate_random_field_element(f);
        uint my_value = (a - b).value();
        uint correct_value = (a.value() + p - b.value()) % p;
        UINT_EQ(my_value, correct_value);
    }
}

TEST(CorrectnessTest, Multiplication) {
    for (size_t i = 0; i < c_correctness_test_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        FieldElement a = generate_random_field_element(f);
        FieldElement b = generate_random_field_element(f);
        uint my_value = (a * b).value();
        uint correct_value = (a.value() * b.value()) % p;
        UINT_EQ(my_value, correct_value);
    }
}

TEST(CorrectnessTest, Division) {
    for (size_t i = 0; i < c_correctness_test_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        FieldElement a = generate_random_field_element(f);
        FieldElement b = generate_random_non_zero_field_element(f);
        FieldElement q = a / b;
        uint my_value = (q * b).value();
        uint correct_value = a.value();
        UINT_EQ(my_value, correct_value);
    }
}

TEST(CorrectnessTest, Inversion) {
    for (size_t i = 0; i < c_correctness_test_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        FieldElement a = generate_random_non_zero_field_element(f);
        FieldElement inv_a = FieldElement::inverse(a);
        FieldElement result = a * inv_a;
        FIELD_EQ(result, f.element(1));
    }
}

TEST(CorrectnessTest, Power) {
    for (size_t i = 0; i < c_correctness_test_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        for (uint j = 1; j < p; j <<= 1) {
            FieldElement a = generate_random_field_element(f);
            uint my_value = FieldElement::pow(a, j).value();

            uint correct_value = 1;
            size_t current_power = 0;

            while (current_power < j) {
                correct_value *= a.value();
                correct_value %= p;
                current_power++;
            }

            UINT_EQ(my_value, correct_value);
        }
    }
}

TEST(CorrectnessTest, Shift) {
    for (size_t i = 0; i < c_correctness_test_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        for (size_t j = 0; j < p; ++j) {
            FieldElement a = generate_random_field_element(f);
            uint my_value = (a << j).value();

            uint correct_value = a.value();
            size_t current_shift = 0;

            while (current_shift < j) {
                correct_value <<= 1;
                correct_value %= p;
                current_shift++;
            }

            UINT_EQ(my_value, correct_value);
        }
    }
}

// Stress tests
TEST(StressTest, Addition) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        for (size_t j = 0; j < c_stress_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement b = generate_random_field_element(f);
            FieldElement result = a + b;
        }
    }
}

TEST(StressTest, Negotiation) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        for (size_t j = 0; j < c_stress_test_negotiation_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement result = -a;
        }
    }
}

TEST(StressTest, Subtraction) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        for (size_t j = 0; j < c_stress_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement b = generate_random_field_element(f);
            FieldElement result = a - b;
        }
    }
}

TEST(StressTest, Multiplication) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        for (size_t j = 0; j < c_stress_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement b = generate_random_field_element(f);
            FieldElement result = a * b;
        }
    }
}

TEST(StressTest, Division) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        for (size_t j = 0; j < c_stress_test_arithmetic_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement b = generate_random_non_zero_field_element(f);
            FieldElement result = a / b;
        }
    }
}

TEST(StressTest, Inversion) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        for (size_t j = 0; j < c_stress_test_inversion_n; ++j) {
            FieldElement a = generate_random_non_zero_field_element(f);
            FieldElement result = FieldElement::inverse(a);
        }
    }
}

TEST(StressTest, Power) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        for (size_t j = 0; j < c_stress_test_power_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            uint power = generate_random_uint_modulo(p);
            FieldElement result = FieldElement::pow(a, power);
        }
    }
}

TEST(StressTest, Shift) {
    for (size_t i = 0; i < c_primes_n; ++i) {
        uint p = get_random_prime();
        Field f(p);

        for (size_t j = 0; j < c_stress_test_shift_n; ++j) {
            FieldElement a = generate_random_field_element(f);
            FieldElement result = a << j;
        }
    }
}
