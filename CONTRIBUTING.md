# Guidelines for Contributing

Below are some C++ programming guidelines for developers who wish to contribute to blitzdg. Some of these are C++ best practices, others are project-specific.

In general, we aim for a good, readable, consistent programming style, so our code is easy to use and easy to understand.

## Header Files 

All header files must begin with a pair of include guards, i.e., `#ifndef ... #define ...`, or with the compiler directive `#pragma once`. When including header files defined in this project put their filenames in double quotes, e.g., `Types.hpp`. Filenames of third party header files and C++ standard library header files are put in angle brackets, e.g., `<blitz/array.h>` or `<vector>`. Include header files defined in this project at the top of the file, followed by third party libraries, and the C++ standard libraries last. Preferably header files should be organized alphabetically within each group, however, this is not mandatory. Only include what you need and put the rest in the source files. If possible, use forward declarations of classes to help reduce the number of included files.

## Namespaces

All code defined in this project is put under the `blitzdg` namespace.

In header files, use fully qualified names of any types not under the `blitzdg` namespaces, e.g.,
`std::vector<double> x;`. Do not include statements such as `using namespace std;` or `using namespace boost;` in a header file. Doing so imports all the names in these namespaces into the file, increasing the likelihood of name clashes.

In source files, use declarations such as `using std::vector;` at global scope as opposed to `using namespace std;`. Doing so allows you to control exactly which names are imported from the included files. Obviously, any names that are already part of the blitzdg namespace do not require any `using blitzdg::` declarations. In certain cases such as with the main function, which must be at global scope, you may add `using namespace blitzdg;` inside the main function, so as to limit its scope. The same practice can be applied to BDD testing with the Igloo testing framework.

## Types

Commonly used types such as `int`, `double`, `Array<double, 1>`, `Array<double, 2>`, `Array<int, 1>`, `Array<int, 2>` have project-specific aliases defined the header `Types.hpp`. These aliases should be used exclusively. The only exception to this caveat is for argument types of double and int in the extern declaration of BLAS, UMFPACK, and LAPACK functions. Since these functions are named based on the real argument type, it is best to be more explicit about their argument types.

## Const vs Non-Const

Use "const correctness" judiciously. For example, if a member function does not modify any member variables of a class then it should be declared `const`. If a member function returns a member variable it should do so by value, by `const` reference, or by pointer to `const`, unless you actually want to be able to modify that member variable. Be careful with returning non-`const` pointers or references since it provides the opportunity to modify the internals of an object, possibly invalidating any class invariants.
