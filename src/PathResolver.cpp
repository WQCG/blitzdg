// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "PathResolver.hpp"
#include "Types.hpp"

using boost::algorithm::find_all;
using boost::algorithm::join;
using boost::algorithm::replace_last;
using boost::algorithm::replace_first;
using boost::algorithm::trim_right;
using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::numeric_limits;
using std::abs;

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

        resolveDelimiter();
    }

    void PathResolver::resolveDelimiter() {
        PathDelimiter = "/";

        if (ExePath.find('\\') != std::string::npos)
            PathDelimiter = "\\";
    }

    string PathResolver::get_RootPath() {
        string path(ExePath);
        vector<string> pathVec;
        replace_last(path, ".exe", "");
        replace_last(path, PathDelimiter + "bin" + PathDelimiter + "test", "");

        return path;
    }

    string PathResolver::joinPaths(string path1, string path2) {
        char d = PathDelimiter.at(0);
        if (path1.back() == d && path2.front() == d)
            path1.pop_back();
        else if (path1.back() != d && path2.front() != d)
            path1.push_back(d);
        path1 += path2;

        return path1;
    }
}
