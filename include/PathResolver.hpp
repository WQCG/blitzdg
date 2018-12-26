// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file PathResolver.hpp
 * @brief Defines the PathResolver class that allows for the resolution of paths
 * relative to the project root path.
 */

#pragma once
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <whereami.h>
#include <boost/algorithm/string.hpp>

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

namespace blitzdg {
  class PathResolver {
  
    string PathDelimeter = "/";
    shared_ptr<string> ExePath = nullptr;
    unique_ptr<vector<string>> InputPathVec = nullptr;

    void resolveDelimeter();

public:
    /**
     * Constructor
     * @note Assumes uniform elements.
     */
    PathResolver();

    /**
     * Copy constructor (deleted).
     */
    PathResolver(const PathResolver&) = delete;

    /**
     * Copy assignment operator (deleted).
     */
    PathResolver& operator=(const PathResolver&) = delete;

    /**
     * Move constructor.
     */
    PathResolver(PathResolver&&) = default;

    /**
     * Move assignment operator.
     */
    PathResolver& operator=(PathResolver&&) = default;

    /**
     * Returns absolute project root path.
     */
    string get_RootPath();

    /**
     * Joins to paths in a way that handles possible trailing/leading delimeters.
     */
    string joinPaths(string path1, string path2);
  };
}
