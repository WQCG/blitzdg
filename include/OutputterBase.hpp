// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#pragma once
#include "Types.hpp"
#include <string>
#include <map>

namespace blitzdg {
    class OutputterBase {
    public:
        virtual void writeFieldsToFiles(std::map<std::string, real_matrix_type>& fields, index_type tstep) = 0;
        virtual ~OutputterBase() = default;
    };
}
