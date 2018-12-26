#include "PathResolver.hpp"
#include "Types.hpp"

namespace blitzdg{
    PathResolver::PathResolver() {
        index_type cap = 1024;
        char pathBuffer[1024];

        // Get path to executable with 'whereami'
        index_type length = wai_getExecutablePath(pathBuffer, cap, NULL);
        ExePath = shared_ptr<string>(new string());

        for(index_type i=0; i < length; i++) {
            *ExePath += pathBuffer[i];
        }
        trim_right(*ExePath);

        resolveDelimeter();
    }

    void PathResolver::resolveDelimeter() {
        PathDelimeter = "/";
        using find_vector_type = vector<iterator_range<string::iterator>>;
        
        string path(*ExePath);
        find_vector_type FindVec;

        find_all( FindVec, path, "\\" );
        if (FindVec.size() > 0) {
            PathDelimeter = "\\";
        }
    }

    string PathResolver::get_RootPath() {
        string path(*ExePath);
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