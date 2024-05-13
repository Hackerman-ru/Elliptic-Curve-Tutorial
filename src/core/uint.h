#ifndef ECG_UINT_H
#define ECG_UINT_H

//#define ECG_USE_BOOST

namespace elliptic_curve_guide {
    namespace uint_info {
        constexpr size_t uint_bits_number = 512;
        constexpr size_t uint_bytes_number = uint_bits_number >> 3;
    }   // namespace uint_info
}   // namespace elliptic_curve_guide

#ifdef ECG_USE_BOOST
    #include "boost/multiprecision/cpp_int.hpp"
    #include "boost/multiprecision/fwd.hpp"

namespace elliptic_curve_guide {
    using uint = boost::multiprecision::uint512_t;
}   // namespace elliptic_curve_guide
#else
    #include "long-arithmetic.h"

namespace elliptic_curve_guide {
    using uint = uint_t<512>;
}   // namespace elliptic_curve_guide
#endif
#endif
