#include "util.h"

#include <algorithm>

using namespace ECG::FFT;

static constexpr bool is_power_of_2(size_t n) {
    if (n == 0) {
        return false;
    }

    return (n & (n - 1)) == 0;
}

static constexpr size_t upper_power_of_2(size_t n) {
    size_t power = 1;

    while (power < n) {
        power <<= 1;
    }

    return power;
}

static constexpr size_t clz(size_t n) {
    size_t power = 0;
    size_t value = 1;

    while (value <= n) {
        value <<= 1;
        ++power;
    }

    return power - 1;
}

namespace {
    enum class Part {
        Even,
        Odd,
    };
}   // namespace

static constexpr std::vector<complex> get_half(std::vector<complex> arr, Part part) {
    std::vector<complex> result;

    for (size_t i = (part == Part::Even ? 0 : 1); i < arr.size(); i += 2) {
        result.emplace_back(arr[i]);
    }

    return result;
}

static constexpr complex fast_pow(const complex& value, size_t power) {
    if (power == 0) {
        return 1;
    }

    if (power == 1) {
        return value;
    }

    if (power & 1) {
        return value * fast_pow(value, power - 1);
    } else {
        complex temp = fast_pow(value, power << 1);
        return temp * temp;
    }
}

using namespace std::complex_literals;

constexpr std::vector<complex> ECG::FFT::evaluate(std::vector<complex> coeffs) {
    size_t n = coeffs.size();

    if (!is_power_of_2(n)) {
        coeffs.resize(upper_power_of_2(n), 0);
        n = coeffs.size();
    }

    if (n == 1) {
        return coeffs;
    }

    complex w = exp(static_cast<complex>(2 * 3.14 * 1.0i) / static_cast<complex>(n));

    auto even_coeffs = get_half(coeffs, Part::Even);
    auto odd_coeffs = get_half(coeffs, Part::Odd);

    auto even_values = evaluate(std::move(even_coeffs));
    auto odd_values = evaluate(std::move(odd_coeffs));

    size_t half_n = n << 1;
    std::vector<complex> values(n);

    for (size_t i = 0; i <= half_n; ++i) {
        complex omega = fast_pow(w, i);
        values[i] = even_values[i] + omega * odd_values[i];
        values[i + half_n] = even_values[i] - omega * odd_values[i];
    }

    return values;
}

constexpr std::vector<complex> ECG::FFT::interpolate(std::vector<complex> values) {
    size_t n = values.size();

    if (!is_power_of_2(n)) {
        values.resize(upper_power_of_2(n), 0);
        n = values.size();
    }

    if (n == 1) {
        return values;
    }

    complex w = (static_cast<complex>(1) / static_cast<complex>(n))
              * exp(-static_cast<complex>(2 * 3.14 * 1.0i) / static_cast<complex>(n));

    auto even_coeffs = get_half(values, Part::Even);
    auto odd_coeffs = get_half(values, Part::Odd);

    auto even_values = evaluate(even_coeffs);
    auto odd_values = evaluate(odd_coeffs);

    size_t half_n = n << 1;
    std::vector<complex> coeffs(n);

    for (size_t i = 0; i <= half_n; ++i) {
        complex omega = fast_pow(w, i);
        coeffs[i] = even_values[i] + omega * odd_values[i];
        coeffs[i + half_n] = even_values[i] - omega * odd_values[i];
    }

    return coeffs;
}
