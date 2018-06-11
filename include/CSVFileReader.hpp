// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file CSVFileReader.hpp
 * @brief Defines the class CSVFileReader and functions csvread for reading
 * csv files.
 * @note The CSVFileReader class does not support multiple line records with
 * quoted text.
 */
#pragma once
#include "Traits.hpp"
#include "Types.hpp"
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace blitzdg {
    namespace details {
        /**
         * Defines a class with static member convert
         * for converting from strings to another type
         * T. Typically T is either a string, char, or
         * numeric type.
         */
        template <typename T>
        struct StrCast {
            static bool convert(const std::string& s, T& ret) {
                ss_.clear(); ss_.str(s);
                if (!(ss_ >> ret))
                    return false;
                return true;
            }
        private:
            static std::istringstream ss_;
        };
        template <typename T> std::istringstream StrCast<T>::ss_;

        /**
        * Specialization of StrCast for type std::string.
        */
        template <>
        struct StrCast<std::string> {
            static bool convert(const std::string& s, std::string& ret) {
                ret = s;
                return true;
            }
        };

        /**
         * Specialization of StrCast for type char.
         */
        template <>
        struct StrCast<char> {
            static bool convert(const std::string& s, char& ret) {
                if (s.size() == 1) {
                    ret = s[0];
                    return true;
                }
                return false;
            }
        };

        /**
         * Specialization of StrCast for type double.
         */
        template <>
        struct StrCast<double> {
            static bool convert(const std::string& s, double& ret) {
                try {
                    size_t lastChar;
                    ret = std::stod(s, &lastChar);
                    return lastChar == s.size();
                }
                catch (std::exception&) {
                    return false;
                }
            }
        };

        /**
         * Specialization of StrCast for type float.
         */
        template <>
        struct StrCast<float> {
            static bool convert(const std::string& s, float& ret) {
                try {
                    size_t lastChar;
                    ret = std::stof(s, &lastChar);
                    return lastChar == s.size();
                }
                catch (std::exception&) {
                    return false;
                }
            }
        };

        /**
         * Specialization of StrCast for type int.
         */
        template <>
        struct StrCast<int> {
            static bool convert(const std::string& s, int& ret) {
                try {
                    size_t lastChar;
                    ret = std::stoi(s, &lastChar);
                    return lastChar == s.size();
                }
                catch(std::exception&) {
                    return false;
                }
            }
        };

        /**
         * Specialization of StrCast for type long long.
         */
        template <>
        struct StrCast<long long> {
            static bool convert(const std::string& s, long long& ret) {
                try {
                    size_t lastChar;
                    ret = std::stoll(s, &lastChar);
                    return lastChar == s.size();
                }
                catch (std::exception&) {
                    return false;
                }
            }
        };

        /**
         * Specialization of StrCast for type unsigned long.
         */
        template <>
        struct StrCast<unsigned long> {
            static bool convert(const std::string& s, unsigned long& ret) {
                try {
                    size_t lastChar;
                    ret = std::stoul(s, &lastChar);
                    return lastChar == s.size();
                }
                catch(std::exception&) {
                    return false;
                }
            }
        };

        /**
         * Specialization of StrCast for type unsigned long long.
         */
        template <>
        struct StrCast<unsigned long long> {
            static bool convert(const std::string& s, unsigned long long& ret) {
                try {
                    size_t lastChar;
                    ret = std::stoull(s, &lastChar);
                    return lastChar == s.size();
                }
                catch(std::exception&) {
                    return false;
                }
            }
        };
    } // namespace details
    
    /**
     * Implements a class for reading csv files.
     */
	class CSVFileReader {
	public:
        /**
         * Constructor.
         * @param[in] filename The path to the csv file to read.
         * @param[in] nskip The number of lines to skip at the start of the file. Defaults to zero.
         * @param[in] delims A string that represents the delimiter sequence. Defaults to "\t " (tab and space).
         * @note The delimiter string may contain any of the following characters: '\t', ' ',
         * ',', ';', '|', '^', and no others.
         */
		explicit CSVFileReader(const std::string& filename, index_type nskip = 0, const std::string& delims = "\t ");
		
		CSVFileReader(const CSVFileReader&) = delete;
		CSVFileReader& operator=(const CSVFileReader&) = delete;
		
        /**
         * Closes any currently open file and opens the csv file with path filename.
         * @see CSVFileReader::CSVFileReader().
         */
		void openFile(const std::string& filename, index_type nskip = 0, const std::string& delims = "\t ");
		
        /**
         * Returns the filename.
         */
		std::string getFilename() const {
			return filename_;
		}
		
        /**
         * Returns the current line number to be read (indexed from one).
         */
		index_type getLineNum() const {
			return lineno_;
		}
		
        /**
         * Returns the number of fields per line in a nonblank and
         * nonempty line.
         */
		index_type getNumCols() const {
			return ncols_;
		}

        /**
         * Returns the number of nonblank and nonempty lines starting
         * from the line nskip.
         * @note This function iterates over the entire file and resets
         * the file stream to the start of line nskip.
         */
        index_type getNumRows();
		
        /**
         * Returns the file stream to the start of line nskip.
         */
		void setToStart() {
			strm_.clear();
			strm_.seekg(0, std::ios::beg);
			lineno_ = 0;
			skipLines(nskip_);
		}
		
        /**
         * Skips the file stream ahead by n lines or until
         * the end-of-file is reached.
         * @param[in] n The number of lines to skip. Should be > 0.
         * @return The state of the file stream. Returns true if 
         * no error flags were set and the EOF has not
         * been reached. Returns false if EOF has been reached.
         */
		bool skipLines(index_type n) {
			if (n > 0) {
                std::string dummy;
                while(n > 0 && readLine(dummy))
                    --n;
            }
			return static_cast<bool>(strm_);
		}
		
        /**
         * Reads the next line of the file if possible.
         * @param[out] line The line of the file that was read.
         * @return The state of the file stream. Returns true if 
         * no error flags were set and the EOF has not
         * been reached. Returns false if EOF has been reached.
         * @note Does not skip blank or empty lines.
         */
		bool readLine(std::string& line);
		
        /**
         * Parses the next line of the file and writes the tokens
         * to the output iterator if possible.
         * @param[out] itr An output iterator.
         * @return The state of the file stream. Returns true if 
         * no error flags were set and the EOF has not
         * been reached. Returns false if EOF has been reached.
         * @note Skips blank or empty lines.
         */
		template <typename OutputItr>
		bool parseRowIterator(OutputItr& itr);
		
        /**
         * Parses the next line of the file and writes the tokens
         * to the output arguments if possible.
         * @param[out] args A list of output arguments.
         * @return The state of the file stream. Returns true if 
         * no error flags were set and the EOF has not
         * been reached. Returns false if EOF has been reached.
         * @note Skips blank or empty lines.
         */
		template <typename... Args>
		bool parseRowValues(Args&... args);
	private:
		std::string filename_; // name of current file
        std::string delims_;   // csv delimiters
		std::ifstream strm_;   // file stream associated with current file
		index_type lineno_;    // current number of lines read (including blank or empty lines)
		index_type nskip_;     // the number of initial lines to skip
		index_type ncols_;     // number of fields per row

        /**
         * Counts the number of fields in the first nonblank and nonempty
         * line starting from line nskip.
         */
        void countCols();

        /**
         * Reads lines until either a nonblank and nonempty line is
         * found or until EOF is reached. Returns true if a line is
         * found and false in the case of EOF.
         */
        bool getNonemptyLine(std::string& line) {
            while (readLine(line)) {
                boost::algorithm::trim(line);
                if (!line.empty())
                    return true;
            }
            return false;
        }
		
        /**
         * Returns true if the delimiter sequence is valid.
         */
        bool checkDelimiters() const {
            return !delims_.empty() && std::all_of(delims_.begin(), delims_.end(),
                [](char c) { return (c == ' ' || c == '\t' || c == ',' 
                || c == ';' || c == '^' || c == '|'); });
        }

        /**
         * Parses a line into string tokens based on the delimiter sequence.
         * @param[in] line The line to be parsed.
         * @param[out] splitVec A vector container of the tokens.
         */
		void tokenizeLine(const std::string& line, std::vector<std::string>& splitVec) const {
			using boost::token_compress_on;
			using boost::algorithm::is_any_of;
			using boost::algorithm::split;
			split(splitVec, line, is_any_of(delims_), token_compress_on);
		}

        /**
         * Converts a string to the templated type T, which
         * should be either and integral or floating point
         * numeric type.
         * @param[in] s The string to convert.
         */
        template <typename T>
        T strCast(const std::string& s) {
            T ret;
            if (!details::StrCast<T>::convert(s, ret))
                throw std::runtime_error("CSVFileReader: conversion failed for '"
                    + s + "' on line " + std::to_string(lineno_) 
                    + " of file " + filename_);
            return ret;
        }
		
        /**
         * Helper function for parsing a line into a sequence of arguments.
         */
		template <typename OutputItr, typename T, typename... Args>
		void parseRowRec(OutputItr strItr, T& first, Args&... rest) {
            first = strCast<T>(*strItr);
		    parseRowRec(++strItr, rest...);
        }
		
        /**
         * Helper function for parsing a line into a sequence of arguments.
         */
		template <typename OutputItr, typename T>
		void parseRowRec(OutputItr strItr, T& val) {
            val = strCast<T>(*strItr);
        }
	};

	template <typename OutputItr>
	bool CSVFileReader::parseRowIterator(OutputItr& itr) {
		using value_type = typename IteratorValueType<OutputItr>::type;
		std::string line;
		if (getNonemptyLine(line)) {
			std::vector<std::string> splitVec;
			tokenizeLine(line, splitVec); // parse line into string tokens stored in splitVec
			if (splitVec.size() != static_cast<size_t>(ncols_)) {
				throw std::runtime_error("CSVFileReader: invalid number of fields on line "
					+ std::to_string(lineno_) + " of file " + filename_);
			}
			for (const auto& s : splitVec) { // write to the output arguments
				*itr = strCast<value_type>(s);
				++itr;
			}
		}
		return static_cast<bool>(strm_);
	}

	template <typename... Args>
	bool CSVFileReader::parseRowValues(Args&... args) {
		if (sizeof...(Args) != static_cast<size_t>(ncols_))
			throw std::runtime_error("CSVFileReader: number of output arguments does not match number of fields");
		std::string line;
		if (getNonemptyLine(line)) {
			std::vector<std::string> splitVec;
			tokenizeLine(line, splitVec); // parse line into string tokens stored in splitVec
			if (splitVec.size() != static_cast<size_t>(ncols_)) {
				throw std::runtime_error("CSVFileReader: invalid number of fields on line "
					+ std::to_string(lineno_) + " of file " + filename_);
			}
			parseRowRec(splitVec.begin(), args...); // write to the output arguments
		}
		return static_cast<bool>(strm_);
	}

    /**
     * Reads a numeric csv file and returns its contents as a unique_ptr<matrix_type<T>>, 
     * where each row of the matrix corresponds to a nonblank and nonempty line in the file.
     * @param[in] filename The path of the file to read.
     * @param[in] nskip The number of lines to skip at the start of the file. Defaults to zero.
     * @param[in] delims A string that represents the delimiter sequence. Defaults to "\t " (tab and space).
     * @note The delimiter string may contain any of the following characters: '\t', ' ',
     * ',', ';', '|', '^', and no others.
     */
	template <typename T>
	std::unique_ptr<matrix_type<T>> 
    csvread(const std::string& filename, index_type nskip = 0, const std::string& delims = "\t ") {
		static_assert(isArithmetic<T>(), "csvread template type T is not a numeric type");
        CSVFileReader reader(filename, nskip, delims);
		std::vector<T> data;
        auto itr = std::back_inserter(data);
        index_type nrows = 0;
		while (reader.parseRowIterator(itr))
            ++nrows;
        std::unique_ptr<matrix_type<T>> ret(new matrix_type<T>(nrows, reader.getNumCols()));
		std::copy(data.begin(), data.end(), ret->begin());
		return ret;
	}

    /**
     * Reads a numeric csv file and returns its contents as a unique_ptr<vector_type<T>>, 
     * in which nonblank and nonempty lines in the file are stored contiguously.
     * @param[in] filename The path of the file to read.
     * @param[out] nrows The number of nonblank and nonempty lines in the file.
     * @param[out] ncols The number of fields in a nonblank and nonempty line.
     * @param[in] nskip The number of lines to skip at the start of the file. Defaults to zero.
     * @param[in] delims A string that represents the delimiter sequence. Defaults to "\t " (tab and space).
     * @note The delimiter string may contain any of the following characters: '\t', ' ',
     * ',', ';', '|', '^', and no others.
     */
    template <typename T>
	std::unique_ptr<vector_type<T>> 
    csvread(const std::string& filename, index_type& nrows, index_type& ncols, index_type nskip = 0, const std::string& delims = "\t ") {
		static_assert(isArithmetic<T>(), "csvread template type T is not a numeric type");
        CSVFileReader reader(filename, nskip, delims);
		nrows = 0;
        ncols = reader.getNumCols();
        std::vector<T> data;
        auto itr = std::back_inserter(data);
		while (reader.parseRowIterator(itr))
            ++nrows;
        std::unique_ptr<vector_type<T>> ret(new vector_type<T>(nrows * ncols));
		std::copy(data.begin(), data.end(), ret->begin());
		return ret;
	}
} // namespace blitzdg
