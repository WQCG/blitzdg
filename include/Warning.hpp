// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file Warning.hpp
 * @brief File containing standard GNU-style copyright and warranty notices.
 */
#include <iostream>

namespace blitzdg {
    /**
     * Prints the product name, version information and standard GPL warning/disclaimer.
     */
    void printDisclaimer() {
        std::cout << "blitzdg, version 1.0b" << std::endl;
        std::cout << "Copyright (C) 2017-2022 Waterloo Quantitative Consulting Group, Inc." << std::endl;
        std::cout << "This is free software; see the source code for copying conditions." << std::endl;
        std::cout << "There is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or" << std::endl;
        std::cout << "FITNESS FOR A PARTICULAR PURPOSE." << std::endl << std::endl;
    }
}