#ifndef ECG_FIELD_ROOT_H
#define ECG_FIELD_ROOT_H

#include "field.h"

#include <optional>

namespace elliptic_curve_guide {
    namespace algorithm {
        std::optional<field::FieldElement> find_root(const field::FieldElement& value,
                                                     const field::Field& field);
    }   // namespace algorithm
}   // namespace elliptic_curve_guide
#endif
