// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "CSVFileReader.hpp"

using std::getline;
using std::runtime_error;
using std::string;
using std::vector;

namespace blitzdg {
	void CSVFileReader::countCols() {
		ncols_ = 0;
		string line;
		if (getNonemptyLine(line)) {
			vector<string> splitVec;
			tokenizeLine(line, splitVec); // parse line into string tokens stored in splitVec
			ncols_ = splitVec.size();
		}
		setToStart();
	}

	index_type CSVFileReader::getNumRows() {
		index_type nrows = 0;
		string line;
		while (getNonemptyLine(line))
			++nrows;
		return nrows;
	}

	CSVFileReader::CSVFileReader(const string& filename, index_type nskip, const string& delims)
		: filename_{ filename }, delims_{ delims }, strm_{ filename }, lineno_{ 0 }, 
        nskip_{ nskip }, ncols_{ -1 }
	{
        if (!checkDelimiters())
            throw runtime_error("CSVFileReader: invalid delimiter");
		if (!strm_.is_open())
			throw runtime_error("CSVFileReader: unable to open file " + filename_);
		if (!skipLines(nskip))
			throw runtime_error("CSVFileReader: number of lines to skip exceeds number of lines in file");
		countCols();
	}

	void CSVFileReader::openFile(const string& filename, index_type nskip, const string& delims) {
		filename_ = filename;
        delims_ = delims;
		strm_.close(); strm_.clear();
		strm_.open(filename);
		lineno_ = 0;
		nskip_ = nskip;
		ncols_ = -1;
        if (!checkDelimiters())
            throw runtime_error("CSVFileReader: invalid delimiter");
		if (!strm_.is_open())
			throw runtime_error("CSVFileReader: unable to open file " + filename_);
		if (!skipLines(nskip))
			throw runtime_error("CSVFileReader: number of lines to skip exceeds number of lines in file");
		countCols();
	}

	void CSVFileReader::setNumCols(index_type cols) {
        ncols_ = cols;
    }

	bool CSVFileReader::readLine(string& line) {
		if (getline(strm_, line))
			++lineno_;
		if (strm_.bad())
			throw runtime_error("CSVFileReader: an error occurred while reading file " + filename_);
		return static_cast<bool>(strm_);
	}
} // namespace blitzdg