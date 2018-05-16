# Guidelines for Contributing

Below are some C++ programming guidelines for developers who wish to contribute to blitzdg. Some of these are C++ best practices, others are project-specific.

In general, we aim for a good, readable, consistent programming style, so our code is easy to use and easy to understand.

## Header files 

When including headers it is common practice to put headers that we wrote in double quotes, e.g., `"Types.hpp"`. Names of third party headers and the C++ standard library names can be put in angle brackets, e.g., `<blitz/array.h>`. Put headers in our project at the top of the file, followed by third party libraries, and the standard libraries at the end.

## Namespaces

In source files, use declarations such as `using std::vector` at global scope as opposed to `using namespace std`. Doing so allows you to see exactly what names are being imported from the included libraries, which is useful. Obviously, if the source code is under the blitzdg namespace, then we don't need any `using blitzdg::` declarations. In certain cases such as with main(), which must be at global scope, add "using namespace blitzdg" inside the main function, so as to limit its scope to that function. This same practice can be applied for igloo tests.

It is best practice to use fully qualified names in header files, i.e., don't use `using...`.

## Types

Instead of explicitly referencing types, use the project-specific aliases. That is, `int`, `double`, `Array<double, 1>`, `Array<double, 2>`, `Array<int, 1>`, `Array<int, 2>` with their corresponding aliases defined in `Types.hpp`. These aliases should be used exclusively. The only exception to this caveat is for double and int was in the extern declaration of BLAS, UMFPACK, and LAPACK functions. Since these functions provide are overloaded based on argument type, it is best to be more explicit about the argument types.

## Const vs Non-Const

If a member function does not modify any internal members of a class then it should be declared `const`. Additionally, if a member function returns a member variable it should do so by value, by `const` reference, or by pointer to `const`, unless you actually want to be able to modify that member variable. Be careful with returning non-`const` pointers or references since it provides the opportunity to mess with the internals of an object, potentially destroying any class invariants.
