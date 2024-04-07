#ifndef ECG_UTIL_H
#define ECG_UTIL_H

#include <algorithm>
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

        constexpr std::vector<complex> evaluate(std::vector<complex> coeffs);
        constexpr std::vector<complex> interpolate(std::vector<complex> values);

        template<typename To>
        static constexpr std::vector<To> convert_vector(std::vector<complex>&& vec) {
            std::vector<To> result(vec.size());
            std::transform(vec.begin(), vec.end(), result.begin(), [](complex value) {
                return static_cast<To>(std::abs(value));
            });
            return result;
        }

        template<typename From, typename To>
        requires is_convertible_to<From, To>
        static constexpr std::vector<To> convert_vector(std::vector<From>&& vec) {
            std::vector<To> result(vec.size());
            std::transform(
                vec.begin(), vec.end(), result.begin(), [](From value) { return static_cast<To>(value); });
            return result;
        }

        template<typename T>
        requires std::is_integral_v<T>
        constexpr std::vector<T> multiply(std::vector<T> lhs, std::vector<T> rhs) {
            size_t degree = lhs.size() + rhs.size();

            lhs.resize(degree);
            rhs.resize(degree);

            std::vector<complex> lhs_values = evaluate(convert_vector<T, complex>(std::move(lhs)));
            std::vector<complex> rhs_values = evaluate(convert_vector<T, complex>(std::move(rhs)));

            std::vector<complex> result_values(degree);

            for (size_t i = 0; i < lhs_values.size(); ++i) {
                result_values[i] = lhs_values[i] * rhs_values[i];
            }

            std::vector<T> result = convert_vector<T>(interpolate(std::move(result_values)));
            return result;
        }

        template<typename T, size_t S>
        constexpr std::array<T, S> convert_to_array(const std::vector<T>& arr) {
            std::array<T, S> result;
            std::copy(std::begin(arr), std::end(arr), std::begin(result));
            return result;
        }

        template<typename T, size_t S>
        constexpr std::vector<T> convert_to_vector(const std::array<T, S>& arr) {
            return std::vector<T>(std::begin(arr), std::end(arr));
        }

        template<typename T, size_t S>
        requires std::is_integral_v<T>
        constexpr std::array<T, S> multiply(const std::array<T, S>& lhs, const std::array<T, S>& rhs) {
            std::vector<T> result = multiply<T>(convert_to_vector<T, S>(lhs), convert_to_vector<T, S>(rhs));
            return convert_to_array<T, S>(result);
        }
    }   // namespace FFT
}   // namespace ECG
#endif
