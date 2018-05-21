// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file CsvOutputter.hpp
 * @brief Defines the CsvOutputter class for writing data to output files.
 */
#pragma once
#include "Types.hpp"
#include <string>

namespace blitzdg {
    class CsvOutputter {

    public:
        void writeFieldToFile(const std::string fileName, const matrix_type & field, const char delimeter);
        std::string generateFileName(const std::string fieldName, const index_type fileNumber);
    };
}