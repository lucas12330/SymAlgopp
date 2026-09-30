# ==========================================
# SymAlgo++ : bibliothèque de calcul symbolique C++17
# ==========================================
#
# Construire, tester, installer :
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build
#   ctest --test-dir build --output-on-failure
#   cmake --install build --prefix /usr/local
#
# Utiliser depuis un autre projet CMake :
#   find_package(symalgopp REQUIRED)
#   target_link_libraries(mon_programme PRIVATE symalgopp::symalgopp)
#
# (le Makefile reste disponible pour le développement : make, make check, make bench)

cmake_minimum_required(VERSION 3.16)

project(symalgopp
    VERSION 0.1.0
    DESCRIPTION "Calcul symbolique C++17 : dérivées, intégrales, équations, fonctions de plusieurs variables"
    HOMEPAGE_URL "https://github.com/lucas12330/SymAlgopp"
    LANGUAGES CXX)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# --- Options -----------------------------------------------------------------

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(SYMALGOPP_EST_RACINE ON)
else()
    set(SYMALGOPP_EST_RACINE OFF) # utilisé par add_subdirectory / FetchContent
endif()

option(BUILD_SHARED_LIBS "Construire une bibliothèque partagée plutôt que statique" OFF)
option(SYMALGOPP_BUILD_TESTS "Construire les tests (ctest)" ${SYMALGOPP_EST_RACINE})
option(SYMALGOPP_BUILD_EXAMPLES "Construire la démo et les cas d'usage" ${SYMALGOPP_EST_RACINE})
option(SYMALGOPP_BUILD_BENCHMARKS "Construire les benchmarks (Google Benchmark et GiNaC requis)" OFF)
option(SYMALGOPP_SANITIZE "Compiler avec AddressSanitizer et UBSan" OFF)
option(SYMALGOPP_INSTALL "Générer les règles d'installation" ${SYMALGOPP_EST_RACINE})

if(SYMALGOPP_EST_RACINE AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Type de construction" FORCE)
endif()

# --- Dépendances -------------------------------------------------------------

# GMP (rationnels exacts de taille arbitraire) : privée à la bibliothèque
find_package(GMP REQUIRED)

# Eigen (matrice compagnon des équations différentielles) : fait partie de l'API publique.
# On préfère la version du système ; à défaut, la copie fournie dans vendor/eigen sert à
# construire (un paquet installé exigera alors Eigen3 chez l'utilisateur).
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(TARGET Eigen3::Eigen)
    set(SYMALGOPP_EIGEN_CONSTRUCTION Eigen3::Eigen)
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/vendor/eigen/Eigen/Dense")
    message(STATUS "Eigen3 introuvable : utilisation de vendor/eigen pour la construction "
                   "(installez libeigen3-dev pour un paquet relogeable)")
    add_library(symalgopp_eigen_interne INTERFACE)
    target_include_directories(symalgopp_eigen_interne SYSTEM INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vendor/eigen>")
    set(SYMALGOPP_EIGEN_CONSTRUCTION symalgopp_eigen_interne)
else()
    message(FATAL_ERROR "Eigen 3.3 ou plus est requis (paquet libeigen3-dev, ou dossier vendor/eigen)")
endif()

# --- Bibliothèque ------------------------------------------------------------

file(GLOB SYMALGOPP_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_library(symalgopp ${SYMALGOPP_SOURCES})
add_library(symalgopp::symalgopp ALIAS symalgopp)

target_compile_features(symalgopp PUBLIC cxx_std_17)
target_include_directories(symalgopp PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/symalgopp>")
target_link_libraries(symalgopp
    PUBLIC
        "$<BUILD_INTERFACE:${SYMALGOPP_EIGEN_CONSTRUCTION}>"
        "$<INSTALL_INTERFACE:Eigen3::Eigen>"
    PRIVATE
        GMP::gmpxx GMP::gmp)
set_target_properties(symalgopp PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR} # ABI instable avant la 1.0
    CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(symalgopp PRIVATE -Wall -Wextra)
endif()

if(SYMALGOPP_SANITIZE)
    set(SYMALGOPP_OPTIONS_SANITIZE -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all)
    target_compile_options(symalgopp PUBLIC ${SYMALGOPP_OPTIONS_SANITIZE} -g)
    target_link_options(symalgopp PUBLIC ${SYMALGOPP_OPTIONS_SANITIZE})
endif()

# --- Tests, démo, cas d'usage, benchmarks ------------------------------------

if(SYMALGOPP_BUILD_TESTS)
    enable_testing()
    file(GLOB SYMALGOPP_TESTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.cpp")
    foreach(source IN LISTS SYMALGOPP_TESTS)
        get_filename_component(nom "${source}" NAME_WE)
        add_executable(${nom} "${source}")
        target_link_libraries(${nom} PRIVATE symalgopp)
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(${nom} PRIVATE -Wall -Wextra)
        endif()
        add_test(NAME ${nom} COMMAND ${nom})
    endforeach()
endif()

if(SYMALGOPP_BUILD_EXAMPLES)
    add_executable(demo demo.cpp)
    target_link_libraries(demo PRIVATE symalgopp)
    file(GLOB SYMALGOPP_EXEMPLES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/examples/*.cpp")
    foreach(source IN LISTS SYMALGOPP_EXEMPLES)
        get_filename_component(nom "${source}" NAME_WE)
        add_executable(exemple_${nom} "${source}")
        target_link_libraries(exemple_${nom} PRIVATE symalgopp)
        # chaque cas d'usage vérifie son propre résultat (code de retour non nul sinon)
        if(SYMALGOPP_BUILD_TESTS)
            add_test(NAME exemple_${nom} COMMAND exemple_${nom})
        endif()
    endforeach()
endif()

if(SYMALGOPP_BUILD_BENCHMARKS)
    find_package(benchmark REQUIRED)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GINAC REQUIRED IMPORTED_TARGET ginac)
    add_executable(bench_suite benchmarks/bench_suite.cpp)
    target_link_libraries(bench_suite PRIVATE symalgopp benchmark::benchmark PkgConfig::GINAC)
    target_compile_options(bench_suite PRIVATE -O3 -DNDEBUG)
endif()

# --- Installation ------------------------------------------------------------

if(SYMALGOPP_INSTALL)
    set(SYMALGOPP_DOSSIER_CMAKE "${CMAKE_INSTALL_LIBDIR}/cmake/symalgopp")

    install(TARGETS symalgopp
        EXPORT symalgoppTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

    # Les en-têtes vont dans <prefix>/include/symalgopp/ : #include <symalgopp> fonctionne
    # (le dossier est dans le chemin d'inclusion de la cible). Canonique.hpp est interne.
    install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/symalgopp
        FILES_MATCHING PATTERN "*.hpp" PATTERN "symalgopp"
        PATTERN "Canonique.hpp" EXCLUDE)

    install(EXPORT symalgoppTargets
        NAMESPACE symalgopp::
        DESTINATION ${SYMALGOPP_DOSSIER_CMAKE})

    configure_package_config_file(cmake/symalgoppConfig.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/symalgoppConfig.cmake"
        INSTALL_DESTINATION ${SYMALGOPP_DOSSIER_CMAKE})
    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/symalgoppConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMinorVersion) # avant la 1.0, chaque version mineure peut casser l'API
    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/symalgoppConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/symalgoppConfigVersion.cmake"
        cmake/FindGMP.cmake
        DESTINATION ${SYMALGOPP_DOSSIER_CMAKE})

    # pkg-config : pour les projets qui n'utilisent pas CMake
    find_package(PkgConfig QUIET)
    set(SYMALGOPP_PC_EIGEN "")
    if(PkgConfig_FOUND)
        pkg_check_modules(SYMALGOPP_EIGEN_PC QUIET eigen3)
        if(SYMALGOPP_EIGEN_PC_FOUND)
            set(SYMALGOPP_PC_EIGEN "eigen3")
        endif()
    endif()
    # GMP n'est requis par les utilisateurs que pour la version statique
    if(BUILD_SHARED_LIBS)
        set(SYMALGOPP_PC_REQUIRES "${SYMALGOPP_PC_EIGEN}")
        set(SYMALGOPP_PC_REQUIRES_PRIVATE "gmpxx gmp")
    else()
        set(SYMALGOPP_PC_REQUIRES "${SYMALGOPP_PC_EIGEN} gmpxx gmp")
        set(SYMALGOPP_PC_REQUIRES_PRIVATE "")
    endif()
    configure_file(cmake/symalgopp.pc.in "${CMAKE_CURRENT_BINARY_DIR}/symalgopp.pc" @ONLY)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/symalgopp.pc"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

    install(FILES README.md CHANGELOG.md DESTINATION ${CMAKE_INSTALL_DOCDIR})

    # cpack : archive .tar.gz (cpack -G TGZ) ou paquet Debian (cpack -G DEB)
    set(CPACK_PACKAGE_NAME symalgopp)
    set(CPACK_PACKAGE_VENDOR "Lucas Bezanilla")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
    set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
    set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
    set(CPACK_GENERATOR TGZ)
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Lucas Bezanilla")
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "libgmp-dev, libgmpxx4ldbl, libeigen3-dev")
    set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
    include(CPack)
endif()
