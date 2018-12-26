// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "PathResolver.hpp"
#include "Types.hpp"

using boost::algorithm::find_all;
using boost::algorithm::join;
using boost::algorithm::replace_last;
using boost::algorithm::replace_first;
using boost::algorithm::trim_right;
using boost::iterator_range;
using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::numeric_limits;
using std::abs;
using std::shared_ptr;
using std::unique_ptr;

namespace blitzdg{
    PathResolver::PathResolver() {
        index_type cap = 1024;
        char pathBuffer[1024];

        // Get path to executable with 'whereami'
        index_type length = wai_getExecutablePath(pathBuffer, cap, NULL);
        ExePath = "";

        for(index_type i=0; i < length; i++) {
            ExePath += pathBuffer[i];
        }
        trim_right(ExePath);

        resolveDelimeter();
    }

    void PathResolver::resolveDelimeter() {
        PathDelimeter = "/";
        using find_vector_type = vector<iterator_range<string::iterator>>;
        
        string path(ExePath);
        find_vector_type FindVec;

        find_all( FindVec, path, "\\" );
        if (FindVec.size() > 0) {
            PathDelimeter = "\\";
        }
    }

    string PathResolver::get_RootPath() {
        string path(ExePath);
        vector<string> pathVec;
        replace_last(path, ".exe", "");
        replace_last(path, PathDelimeter + "bin" + PathDelimeter + "test", "");

        return path;
    }

    string PathResolver::joinPaths(string path1, string path2) {

        char delim = *PathDelimeter.cbegin();
        if (*path2.cbegin() == delim)
            replace_last(path2, PathDelimeter, "");
        
        if (*path1.cend() == delim)
            replace_first(path1, PathDelimeter, "");

        return path1 + PathDelimeter + path2;
    }
}