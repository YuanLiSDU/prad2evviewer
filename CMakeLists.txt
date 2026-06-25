cmake_minimum_required(VERSION 3.14)
project(prad2evviewer VERSION 1.0.0 LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Layout
#   cmake/           helper modules (PradHelpers, WebDeps)
#   prad2dec/        evio/ET decoder library
#   prad2det/        HyCal + GEM reconstruction
#   python/          prad2py pybind11 module          (BUILD_PYTHON)
#   analysis/        ROOT replay tools + data-source  (BUILD_ANALYSIS)
#   test/            command-line dump/diagnostic exes
#   src/             web server + Qt GUI applications (BUILD_GUI)
# ---------------------------------------------------------------------------

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(PradHelpers)        # prad_option, prad_optional_probe, add_prad_tool

# ---------------------------------------------------------------------------
# Global compile settings
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -fPIC for all static libs so they can be pulled into shared objects
# (ACLiC-compiled ROOT macros and the pybind11 module).
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# Compile-time paths (auto-discovered, overridable)
# ---------------------------------------------------------------------------
if(NOT DEFINED DATABASE_DIR)
    set(DATABASE_DIR "${CMAKE_SOURCE_DIR}/database")
endif()
get_filename_component(DATABASE_DIR "${DATABASE_DIR}" ABSOLUTE)

if(NOT DEFINED RESOURCE_DIR)
    set(RESOURCE_DIR "${CMAKE_SOURCE_DIR}/resources")
endif()
get_filename_component(RESOURCE_DIR "${RESOURCE_DIR}" ABSOLUTE)

# ---------------------------------------------------------------------------
# Install-dir conventions — pulled up here so every add_subdirectory()
# below can use CMAKE_INSTALL_BINDIR / CMAKE_INSTALL_DATADIR / _LIBDIR in
# its own install() calls.
# ---------------------------------------------------------------------------
include(GNUInstallDirs)

# ---------------------------------------------------------------------------
# HTTP/WebSocket dependencies (json, websocketpp, asio — all header-only).
# Loaded *before* the core libraries so prad2dec / prad2det can pull
# nlohmann_json into their own targets (e.g. PulseTemplateStore.cpp).
# ---------------------------------------------------------------------------
include(WebDeps)

# ---------------------------------------------------------------------------
# Core libraries (always built)
# ---------------------------------------------------------------------------
add_subdirectory(prad2dec)                    # libprad2dec.a + WITH_ET option
add_subdirectory(prad2det)                    # libprad2det.a

# ---------------------------------------------------------------------------
# Diagnostic command-line tools
# ---------------------------------------------------------------------------
add_subdirectory(test)    # evio_dump, ped_calc + test/dev tools
add_subdirectory(gem)     # gem_dump + Python GEM tools

# ---------------------------------------------------------------------------
# Optional feature: ROOT-based analysis tools + viewer data source
# ---------------------------------------------------------------------------
prad_option(BUILD_ANALYSIS ON
    "Build analysis tools and ROOT data-source (requires ROOT)")

if(BUILD_ANALYSIS)
    find_package(ROOT 6.0 QUIET CONFIG COMPONENTS Core Tree RIO)
    prad_optional_probe(BUILD_ANALYSIS FOUND ROOT_FOUND
        HINT "Install ROOT 6.x (https://root.cern) and add its bin/ to PATH.")
endif()
if(BUILD_ANALYSIS)
    include(${ROOT_USE_FILE})
    message(STATUS "ROOT: using ${ROOT_VERSION}")
endif()
# prad2srv reads WITH_ROOT to switch root_data_source.cpp in/out.
set(WITH_ROOT ${BUILD_ANALYSIS})

# ---------------------------------------------------------------------------
# Server library (linked by prad2_server)
# ---------------------------------------------------------------------------
find_package(Threads REQUIRED)
# zlib for HTTP gzip on big JSON payloads (gem_apv ~1.3 MB/event).  Comes
# pre-installed on most Linux distros and Windows toolchains.
find_package(ZLIB REQUIRED)

set(SRV_SOURCES
    src/viewer_server.cpp src/viewer_server_http.cpp src/viewer_server_et.cpp
    src/app_state.cpp     src/app_state_init.cpp     src/app_state_api.cpp
    src/app_state_filter.cpp
    src/data_source.cpp   src/evio_data_source.cpp)
if(WITH_ROOT)
    list(APPEND SRV_SOURCES src/root_data_source.cpp)
endif()

add_library(prad2srv STATIC ${SRV_SOURCES})
target_include_directories(prad2srv PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(prad2srv PUBLIC
    prad2dec prad2det
    nlohmann_json::nlohmann_json websocketpp_lib Threads::Threads ZLIB::ZLIB)
if(WITH_ET)
    target_compile_definitions(prad2srv PUBLIC WITH_ET)
endif()
if(WITH_ROOT)
    target_compile_definitions(prad2srv PUBLIC WITH_ROOT)
    target_link_libraries(prad2srv PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)
endif()

# ---------------------------------------------------------------------------
# Unified HTTP server
# ---------------------------------------------------------------------------
add_executable(prad2_server src/prad2_server.cpp)
target_compile_definitions(prad2_server PRIVATE
    DATABASE_DIR="${DATABASE_DIR}" RESOURCE_DIR="${RESOURCE_DIR}")
target_link_libraries(prad2_server PRIVATE prad2srv)
set_target_properties(prad2_server PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------------------------------------------------------------------
# Analysis subdirectory (ROOT-dependent producers)
# ---------------------------------------------------------------------------
if(BUILD_ANALYSIS)
    add_subdirectory(analysis)
endif()

# ---------------------------------------------------------------------------
# Optional feature: Python bindings (prad2py)
# ---------------------------------------------------------------------------
prad_option(BUILD_PYTHON ON "Build Python bindings (prad2py)")

if(BUILD_PYTHON)
    find_package(Python 3.7 QUIET COMPONENTS Interpreter Development)
    set(_PY_OK FALSE)
    if(Python_FOUND AND Python_Development_FOUND)
        set(_PY_OK TRUE)
    endif()
    prad_optional_probe(BUILD_PYTHON FOUND _PY_OK HINT
        "RHEL/Rocky:  sudo dnf install python3-devel\n"
        "Debian/Ubuntu:  sudo apt install python3-dev\n"
        "Or point CMake at a specific interpreter with "
        "-DPython_EXECUTABLE=$(which python3).")
endif()
if(BUILD_PYTHON)
    message(STATUS "Python: ${Python_EXECUTABLE} "
                   "(headers: ${Python_INCLUDE_DIRS})")
    add_subdirectory(python)
endif()

# ---------------------------------------------------------------------------
# Optional feature: Qt GUI applications (prad2_client)
# ---------------------------------------------------------------------------
# The Python GUIs (hycal_event_viewer, gem_event_viewer) cover per-event
# inspection; the remaining Qt-native binary is prad2_client, a thin
# WebEngine wrapper for a remote prad2_server.
prad_option(BUILD_GUI ON "Build Qt GUI application (prad2_client)")

if(BUILD_GUI)
    find_package(Qt6 QUIET COMPONENTS Widgets WebEngineWidgets)
    if(Qt6_FOUND)
        set(QT_VERSION_MAJOR 6)
    else()
        find_package(Qt5 QUIET COMPONENTS Widgets WebEngineWidgets)
        if(Qt5_FOUND)
            set(QT_VERSION_MAJOR 5)
        endif()
    endif()
    set(_QT_OK ${Qt6_FOUND})
    if(NOT _QT_OK)
        set(_QT_OK ${Qt5_FOUND})
    endif()
    prad_optional_probe(BUILD_GUI FOUND _QT_OK HINT
        "Install Qt6 Widgets + WebEngineWidgets (preferred) or Qt5.")
endif()
if(BUILD_GUI)
    message(STATUS "Qt: using Qt${QT_VERSION_MAJOR}")

    # Thin remote client.
    add_executable(prad2_client src/prad2_client.cpp)
    target_link_libraries(prad2_client PRIVATE
        Qt${QT_VERSION_MAJOR}::Widgets
        Qt${QT_VERSION_MAJOR}::WebEngineWidgets)
    set_target_properties(prad2_client PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif()

# ---------------------------------------------------------------------------
# Install
# ---------------------------------------------------------------------------
install(TARGETS prad2_server evio_dump ped_calc
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
if(BUILD_GUI)
    install(TARGETS prad2_client
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

install(DIRECTORY database/  DESTINATION ${CMAKE_INSTALL_DATADIR}/prad2evviewer/database)
install(DIRECTORY resources/ DESTINATION ${CMAKE_INSTALL_DATADIR}/prad2evviewer/resources)

# prad2_setup.csh is generated from scripts/prad2_setup.csh.in so
# @CMAKE_INSTALL_PREFIX@ can be substituted — csh lacks a BASH_SOURCE
# equivalent, so the script can't reliably find itself under nested
# `source`; the baked-in prefix is the fallback when self-detection
# fails.  prad2_setup.sh uses BASH_SOURCE and doesn't need this treatment.
configure_file(${CMAKE_SOURCE_DIR}/scripts/shell/prad2_setup.csh.in
    ${CMAKE_BINARY_DIR}/bin/prad2_setup.csh @ONLY)
install(FILES scripts/shell/prad2_setup.sh ${CMAKE_BINARY_DIR}/bin/prad2_setup.csh
    DESTINATION ${CMAKE_INSTALL_BINDIR})

# ---------------------------------------------------------------------------
# Install the scripts/ tree (non-GEM Python tools + auxiliary data files)
# under share/prad2evviewer/scripts/.  GEM tools live in gem/ and get their
# own install rule in gem/CMakeLists.txt.  Sibling imports still work
# because each tool runs from inside its own directory, which Python puts
# on sys.path.
# ---------------------------------------------------------------------------
install(DIRECTORY scripts/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/prad2evviewer/scripts
    FILES_MATCHING
        PATTERN "*.py"
        PATTERN "*.json"
        PATTERN "*.csv"
        PATTERN "*.txt"
        PATTERN "README.md"
        PATTERN "__pycache__" EXCLUDE
        PATTERN "daq_tool"    EXCLUDE   # dev-only GUIs, not installed
        PATTERN "dev_tool"    EXCLUDE   # dev-only one-shot helpers
        PATTERN "shell"       EXCLUDE   # operator shell scripts + prad2_setup.*
                                        #  (setup files live here and are
                                        #   installed explicitly via install(FILES …))
)

# ---------------------------------------------------------------------------
# Generate + install thin bin/ wrappers for each user-facing Python tool so
# they can be launched without needing a sourced prad2_setup.sh (the wrappers set
# PRAD2_DATABASE_DIR / PRAD2_RESOURCE_DIR / PYTHONPATH on their own).
#
# The script path is relative to <prefix>/share/prad2evviewer/ — GEM tools
# live directly under gem/, the rest under scripts/.
# ---------------------------------------------------------------------------
if(BUILD_PYTHON)
    set(_PRAD2_PY_WRAPPERS
        # "wrapper_name:relative_script_path"
        "gem_event_viewer:gem/gem_event_viewer.py"
        "gem_cluster_view:gem/gem_cluster_view.py"
        "gem_layout:gem/gem_layout.py"
        "tagger_viewer:scripts/tagger_viewer.py"
        "hycal_scaler_map:scripts/hycal_scaler_map.py"
        "hycal_pedestal_monitor:scripts/hycal_pedestal_monitor.py"
        "hycal_gain_monitor:scripts/hycal_gain_monitor.py"
        "online_gain_monitor:scripts/online_gain_monitor.py"
        "hycal_map_builder:scripts/hycal_map_builder.py"
        "hycal_event_viewer:scripts/hycal_event_viewer.py"
        "gem_hycal_match_viewer:scripts/gem_hycal_match_viewer.py"
    )
    foreach(_entry ${_PRAD2_PY_WRAPPERS})
        string(REPLACE ":" ";" _parts "${_entry}")
        list(GET _parts 0 _name)
        list(GET _parts 1 PRAD2_WRAPPER_SCRIPT_REL)
        configure_file(
            ${CMAKE_SOURCE_DIR}/cmake/prad2_py_wrapper.sh.in
            ${CMAKE_BINARY_DIR}/bin/${_name}
            @ONLY
        )
        install(PROGRAMS ${CMAKE_BINARY_DIR}/bin/${_name}
            DESTINATION ${CMAKE_INSTALL_BINDIR}
        )
    endforeach()
    unset(PRAD2_WRAPPER_SCRIPT_REL)
endif()

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message(STATUS "----------------------------------------------------")
message(STATUS "  ${PROJECT_NAME} ${PROJECT_VERSION}")
message(STATUS "  Build type   : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Database     : ${DATABASE_DIR}")
message(STATUS "  Resources    : ${RESOURCE_DIR}")
message(STATUS "  ET support   : ${WITH_ET}")
message(STATUS "  Analysis+ROOT: ${BUILD_ANALYSIS}")
message(STATUS "  Python bind. : ${BUILD_PYTHON}")
message(STATUS "  Qt GUI       : ${BUILD_GUI}")
message(STATUS "----------------------------------------------------")
