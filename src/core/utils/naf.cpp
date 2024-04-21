#include "naf.h"

namespace ECG {
    auto get_naf(const uint& value) {
        uint half_value = value >> 1;
        uint one_and_a_half_value = value + half_value;
        uint c = half_value ^ one_and_a_half_value;
        uint positive_bits = one_and_a_half_value & c;
        uint negative_bits = half_value & c;
    }

}   // namespace ECG
