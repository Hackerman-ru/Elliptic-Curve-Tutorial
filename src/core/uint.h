#ifndef ECG_UINT_H
#define ECG_UINT_H

#define ECG_USE_BOOST

#ifdef ECG_USE_BOOST
    #include "boost/multiprecision/cpp_int.hpp"
    #include "boost/multiprecision/fwd.hpp"

namespace ECG {
    using uint = boost::multiprecision::uint512_t;
}   // namespace ECG
#else
    #include "long-arithmetic.h"

namespace ECG {
    using uint = uint_t<512>;
}   // namespace ECG
#endif
#endif
