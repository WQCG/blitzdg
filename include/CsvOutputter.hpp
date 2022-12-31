// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file CsvOutputter.hpp
 * @brief Defines the CsvOutputter class for writing data to output files.
 */
#pragma once
#include "Types.hpp"
#include "OutputterBase.hpp"
#include <string>

namespace blitzdg {
    class CsvOutputter : OutputterBase {

    public:
        void writeFieldToFile(const std::string & fileName, const real_matrix_type & field, const char delimeter);
        std::string generateFileName(const std::string & fieldName, const index_type fileNumber);
        void writeFieldsToFiles(std::map<std::string, real_matrix_type>& fields, index_type tstep);
    };
}