#ifndef ECG_UTIL_H
#define ECG_UTIL_H

#include <complex>
#include <concepts>
#include <vector>

namespace ECG {
    template<typename From, typename To>
    concept is_convertible_to = requires(From f) {
        { static_cast<To>(f) } noexcept;
    };

    template<typename From, typename To>
    concept is_upcastable_to = sizeof(From) <= sizeof(To) && is_convertible_to<From, To>;

    template<typename From, typename To>
    concept is_downcastable_to = sizeof(From) > sizeof(To) && is_convertible_to<From, To>;

    template<typename Container, typename T>
    concept is_convertible_container = requires(Container t, size_t i) {
        { t[i] } -> is_convertible_to<T>;
        { t.size() } -> std::same_as<size_t>;
    };

    template<typename T>
    concept is_integral = std::is_integral_v<T> || requires(T t, T* p, void (*f)(T)) {
        f(0);
        p + t;
    };

    namespace FFT {
        using complex = std::complex<long double>;

        static constexpr std::vector<complex> evaluate(std::vector<complex> coeffs);
        static constexpr std::vector<complex> interpolate(std::vector<complex> values);

        template<typename T>
        requires std::is_integral_v<T>
        constexpr std::vector<T> multiply(std::vector<T> lhs, std::vector<T> rhs) {
            size_t degree = lhs.size() + rhs.size();

            lhs.resize(degree);
            rhs.resize(degree);

            std::vector<complex> lhs_values = evaluate(std::move(lhs));
            std::vector<complex> rhs_values = evaluate(std::move(rhs));

            std::vector<complex> result_values;

            for (size_t i = 0; i < lhs_values.size(); ++i) {
                result_values[i] = lhs_values[i] * rhs_values[i];
            }

            std::vector<T> result = interpolate(std::move(result_values));
            return result;
        }
    }   // namespace FFT
}   // namespace ECG
#endif
