// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "CsvOutputter.hpp"
#include "Types.hpp"
#include <blitz/array.h>
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>

using std::ofstream;
using std::stringstream;
using std::setfill;
using std::setw;
using std::string;
using std::endl;

namespace blitzdg {

    /**
     * Writes a blitz array to plain-text file.
     * @param[in] fileName Name of the file (e.g., field0000010.dat).
     * @param[in] field Two-dimensional blitz array to be written to the file. Usually a 'field' of the PDE (system) being solved.
     * @param[in] delimeter Character that will be used to separate columns. Rows are always separated by line-endings.
     */
    void CsvOutputter::writeFieldToFile(const string & fileName, const real_matrix_type & field, const char delimeter) {
        ofstream outFile;
        outFile.open(fileName);
        for(index_type i=0; i < field.rows(); i++) {
            for(index_type k=0; k < field.cols(); k++) {
                real_type num = field(i, k);
                outFile << num << delimeter;
            }
            outFile << endl;
        }
        outFile.close();
    }

    /**
     * Generates a file name for storing a blitzdg object in plain-text.
     * @param[in] fieldName Name of field from the PDE (system), e.g., "u".
     * @param[in] fileNumber An integral intex indicating a logical ordering on the output files. It is usually related to time-level.
     */
    string CsvOutputter::generateFileName(const string & fieldName, const index_type fileNumber) {
        stringstream fileNameStrm;
        fileNameStrm << "u" << setfill('0') << setw(7) << fileNumber << ".dat";
        return fileNameStrm.str();
    }
}