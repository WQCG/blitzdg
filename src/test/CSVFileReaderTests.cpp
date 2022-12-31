// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "CSVFileReader.hpp"
#include "LinAlgHelpers.hpp"
#include "Types.hpp"
#include <igloo/igloo_alt.h>
#include <boost/algorithm/string.hpp>
#include <whereami.h>
#include <memory>
#include <string>
#include <vector>

using boost::algorithm::find_all;
using boost::algorithm::join;
using boost::algorithm::replace_last;
using boost::algorithm::trim_right;
using boost::iterator_range;
using std::cout;
using std::endl;
using std::string;
using std::unique_ptr;
using std::vector;

namespace blitzdg {
    namespace CSVFileReaderTests {
        using namespace igloo;
        Describe(CSVFileReader_Object) {
            using find_vector_type = vector<iterator_range<string::iterator>>;

            string ExePath;

            void SetUp() {
                if (ExePath.empty()) {
                    // Deal with paths to the test input files.
                    index_type cap = 1024;
                    char* pathBuffer = new char[cap];
                    index_type length = wai_getExecutablePath(pathBuffer, cap, NULL);

                    for(index_type i = 0; i < length; ++i) {
                        ExePath += pathBuffer[i]; 
                    }
                    trim_right(ExePath);
                    delete[] pathBuffer;
                }
            }

            string getFilePath(const string& fname) {    
                find_vector_type FindVec;
                string PathDelimeter = "/";
                string path(ExePath);
                replace_last(path, ".exe", "");
                replace_last(path, "/bin/test", "");
                find_all(FindVec, path, "\\");
                if (FindVec.size() > 0) {
                    PathDelimeter = "\\";
                    replace_last(path, "\\bin\\test", "");
                }

                vector<string> pathVec;
                pathVec.push_back(path);
                pathVec.push_back("input");
                pathVec.push_back(fname);

                return join(pathVec, PathDelimeter);
            }

            It(Reads_A_CSV_File) {
                cout << "CSVFileReaderTests: Read a CSV file" << endl;
                string fname = getFilePath("csvtest1.csv");
                CSVFileReader reader(fname);
                index_type nrows = reader.getNumRows();
                index_type ncols = reader.getNumCols();
                
                string line;
                while (reader.readLine(line));
                
                Assert::That(nrows, Equals(6));
                Assert::That(ncols, Equals(2));
                Assert::That(reader.getFilename(), Equals(fname));
                Assert::That(reader.getLineNum(), Equals(10));
            }

            It(Reads_A_CSV_File_To_A_Matrix) {
                cout << "CSVFileReaderTests: Read a CSV file to a matrix" << endl;
                string fname = getFilePath("csvtest1.csv");
                unique_ptr<matrix_type<double>> out = csvread<double>(fname);
                matrix_type<double> ret(6,2);
                ret = 0.0,0.0,
                      0.5,0.0,
                      1.0,0.0,
                      1.0,1.0,
                      0.5,1.0,
                      0.0,1.0;
                ret -= *out;
                Assert::That(normMax(ret), Equals(0));
            }

            It(Reads_A_CSV_File_To_A_Vector) {
                cout << "CSVFileReaderTests: Read a CSV file to a vector" << endl;
                string fname = getFilePath("csvtest1.csv");
                index_type nrows, ncols;
                unique_ptr<vector_type<double>> out = csvread<double>(fname, nrows, ncols);
                vector_type<double> ret(12);
                ret = 0.0,0.0,
                      0.5,0.0,
                      1.0,0.0,
                      1.0,1.0,
                      0.5,1.0,
                      0.0,1.0;
                ret -= *out;
                Assert::That(normInf(ret), Equals(0));
            }

            It(Skips_The_First_Few_Lines) {
                cout << "CSVFileReaderTests: Read a CSV file skipping the first few lines" << endl;
                string fname = getFilePath("csvtest1.csv");
                unique_ptr<matrix_type<double>> out = csvread<double>(fname, 4);
                matrix_type<double> ret(4,2);
                ret = 1.0,0.0,
                      1.0,1.0,
                      0.5,1.0,
                      0.0,1.0;
                ret -= *out;
                Assert::That(normMax(ret), Equals(0.0));
            }

            It(Reads_Different_Field_Types) {
                cout << "CSVFileReaderTests: Read a CSV file with different field types" << endl;
                string fname = getFilePath("csvtest2.csv");
                CSVFileReader reader(fname, 1);
                vector<string> name;
                vector<float> price;
                vector<int> number;
                string nm;
                float prc;
                int num;
                while (reader.parseRowValues(nm, prc, num)) {
                    name.push_back(nm);
                    price.push_back(prc);
                    number.push_back(num);
                }
                Assert::That(name.size(), Equals(5));
                Assert::That(name[0], Equals("pencil"));
                Assert::That(price[0], Equals(1.0));
                Assert::That(number[0], Equals(20));
                Assert::That(name[1], Equals("scissors"));
                Assert::That(price[1], Equals(2.0));
                Assert::That(number[1], Equals(10));
                Assert::That(name[2], Equals("booklets"));
                Assert::That(price[2], Equals(1.5));
                Assert::That(number[2], Equals(20));
                Assert::That(name[3], Equals("erasers"));
                Assert::That(price[3], Equals(0.75));
                Assert::That(number[3], Equals(20));
                Assert::That(name[4], Equals("crayons"));
                Assert::That(price[4], Equals(3.75));
                Assert::That(number[4], Equals(15));
            }

            It(Opens_A_New_File_To_Read) {
                cout << "CSVFileReaderTests: Open new CSV file to read" << endl;
                string fname1 = getFilePath("csvtest1.csv");
                string fname2 = getFilePath("csvtest2.csv");
                CSVFileReader reader(fname2, 1);
                reader.openFile(fname1);
                vector_type<double> ret(12);
                ret = 0.0,0.0,
                      0.5,0.0,
                      1.0,0.0,
                      1.0,1.0,
                      0.5,1.0,
                      0.0,1.0;
                vector_type<double> out(12);
                vector_type<double>::iterator itr = out.begin();
                while (reader.parseRowIterator(itr));
                ret -= out;
                Assert::That(normInf(ret), Equals(0.0));
            }

            It(Detects_A_Missing_Field) {
                cout << "CSVFileReaderTests: Detects a missing field" << endl;
                string fname = getFilePath("csvtest3.csv");
                bool caught = false;
                try {
                    unique_ptr<matrix_type<double>> out = csvread<double>(fname);
                }
                catch (std::runtime_error& e) {
                    cout << e.what() << endl;
                    caught = true;
                }
                Assert::That(caught, Equals(true));
            }

            It(Detects_An_Invalid_Field) {
                cout << "CSVFileReaderTests: Detects an invalid field" << endl;
                string fname = getFilePath("csvtest4.csv");
                bool caught = false;
                try {
                    unique_ptr<matrix_type<double>> out = csvread<double>(fname);
                }
                catch (std::runtime_error& e) {
                    cout << e.what() << endl;
                    caught = true;
                }
                Assert::That(caught, Equals(true));
            }

            It(Detects_Too_Many_Skipped_Lines) {
                cout << "CSVFileReaderTests: Detects too many skipped lines" << endl;
                string fname = getFilePath("csvtest1.csv");
                bool caught = false;
                try {
                    CSVFileReader reader(fname, 11);
                }
                catch (std::runtime_error& e) {
                    cout << e.what() << endl;
                    caught = true;
                }
                Assert::That(caught, Equals(true));
                caught = false;
                 try {
                    CSVFileReader reader(fname);
                    reader.openFile(fname, 11);
                }
                catch (std::runtime_error& e) {
                    cout << e.what() << endl;
                    caught = true;
                }
                Assert::That(caught, Equals(true));
            }

            It(Detects_Invalid_Delimiter_Sequence) {
                cout << "CSVFileReaderTests: Detects invalid delimiter sequence" << endl;
                string fname = getFilePath("csvtest1.csv");
                bool caught = false;
                try {
                    CSVFileReader reader(fname, 0, "~$#@*|^8x,");
                }
                catch (std::runtime_error& e) {
                    cout << e.what() << endl;
                    caught = true;
                }
                Assert::That(caught, Equals(true));
            }

             It(Detects_Invalid_File) {
                cout << "CSVFileReaderTests: Detects the file could not be opened" << endl;
                bool caught = false;
                try {
                    CSVFileReader reader("not_a_file.csv");
                }
                catch (std::runtime_error& e) {
                    cout << e.what() << endl;
                    caught = true;
                }
                Assert::That(caught, Equals(true));
            }
        };
    } // namespace CSVFileReaderTests
} // namespace blitzdg