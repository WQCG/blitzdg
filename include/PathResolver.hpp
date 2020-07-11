// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc. 
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

namespace blitzdg {
  class PathResolver {
  
    std::string PathDelimiter;
    std::string ExePath;

    void resolveDelimiter();

public:
    /**
     * Constructor
     */
    PathResolver();

    /**
     * Copy constructor (deleted).
     */
    PathResolver(const PathResolver&) = delete;

    /**
     * Copy assignment operator
     */
    PathResolver& operator=(const PathResolver&) = default;

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
    std::string get_RootPath();

    /**
     * Joins to paths in a way that handles possible trailing/leading delimeters.
     */
    std::string joinPaths(std::string path1, std::string path2);
  };
}
