#ifndef FMT_SUPPORT_HPP
#define FMT_SUPPORT_HPP

// GCC < 13 doesn't have std::format. We use fmtlib as a workaround.
#if __has_include(<format>)
#include <format>
#endif

#if defined(__cpp_lib_format) && __cpp_lib_format >= 201907L
using std::format;
#else
#include <fmt/format.h>
using fmt::format;
#endif

#endif
