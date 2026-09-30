# FindGMP.cmake : GMP et son interface C++ (gmpxx).
#
# Cibles importées : GMP::gmp et GMP::gmpxx (qui dépend de GMP::gmp).
# Variables de recherche : GMP_ROOT, ou les chemins habituels (CMAKE_PREFIX_PATH).

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(PC_GMP QUIET gmp)
    pkg_check_modules(PC_GMPXX QUIET gmpxx)
endif()

find_path(GMP_INCLUDE_DIR gmp.h HINTS ${PC_GMP_INCLUDE_DIRS} ${GMP_ROOT}/include)
find_path(GMPXX_INCLUDE_DIR gmpxx.h HINTS ${PC_GMPXX_INCLUDE_DIRS} ${GMP_ROOT}/include)
find_library(GMP_LIBRARY NAMES gmp HINTS ${PC_GMP_LIBRARY_DIRS} ${GMP_ROOT}/lib)
find_library(GMPXX_LIBRARY NAMES gmpxx HINTS ${PC_GMPXX_LIBRARY_DIRS} ${GMP_ROOT}/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
    REQUIRED_VARS GMP_LIBRARY GMP_INCLUDE_DIR GMPXX_LIBRARY GMPXX_INCLUDE_DIR)

if(GMP_FOUND)
    if(NOT TARGET GMP::gmp)
        add_library(GMP::gmp UNKNOWN IMPORTED)
        set_target_properties(GMP::gmp PROPERTIES
            IMPORTED_LOCATION "${GMP_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}")
    endif()
    if(NOT TARGET GMP::gmpxx)
        add_library(GMP::gmpxx UNKNOWN IMPORTED)
        set_target_properties(GMP::gmpxx PROPERTIES
            IMPORTED_LOCATION "${GMPXX_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMPXX_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES GMP::gmp)
    endif()
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMPXX_INCLUDE_DIR GMP_LIBRARY GMPXX_LIBRARY)
