#include "Field/Field.h"
#include "Uint/uint.h"

#include <iostream>

int main() {
    ECG::uint_t<512> a = 1;
    std::cout << a.into_string();
}
