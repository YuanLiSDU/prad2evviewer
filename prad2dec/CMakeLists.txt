# prad2dec/CMakeLists.txt — static library: libprad2dec.a
#
# Can be used standalone:  add_subdirectory(prad2dec)
# Or as part of the top-level project.

cmake_minimum_required(VERSION 3.14)

# If built standalone, set up project and standards
if(NOT CMAKE_PROJECT_NAME)
    project(prad2dec VERSION 1.0.0 LANGUAGES CXX)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

# ---------------------------------------------------------------------------
# Dependencies: evio and et
# ---------------------------------------------------------------------------
set(EVIO_SOURCE "hallb" CACHE STRING "evio provider: hallb | fetch | system")
set(ET_SOURCE   "hallb" CACHE STRING "et provider:   hallb | fetch | system")

set(CODA_ROOT "/usr/clas12/release/2.0.0/coda/Linux_x86_64_RHEL9"
    CACHE PATH "Root of the Hall-B CODA installation")
set(CODA_LIB_DIR     "${CODA_ROOT}/lib"     CACHE PATH "Hall-B lib dir")
set(CODA_INCLUDE_DIR "/usr/clas12/release/2.0.0/coda/common/include"
    CACHE PATH "CODA include dir")

include(FetchContent)
set(FETCHCONTENT_QUIET OFF CACHE BOOL "" FORCE)

# --- evio -----------------------------------------------------------------
if(EVIO_SOURCE STREQUAL "hallb")
    find_library(EVIO_LIB NAMES evio PATHS ${CODA_LIB_DIR} NO_DEFAULT_PATH)
    if(EVIO_LIB)
        find_library(EXPAT_LIB NAMES expat REQUIRED)
        add_library(evio STATIC IMPORTED GLOBAL)
        set_target_properties(evio PROPERTIES
            IMPORTED_LOCATION             "${EVIO_LIB}"
            INTERFACE_INCLUDE_DIRECTORIES "${CODA_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES      "${EXPAT_LIB}"
        )
        message(STATUS "evio: using Hall-B installation at ${EVIO_LIB}")
    else()
        message(STATUS "evio: Hall-B library not found, falling back to fetch")
        set(EVIO_SOURCE "fetch")
    endif()
endif()
if(EVIO_SOURCE STREQUAL "fetch")
    # We only need the C library; skip C++ eviocc (which requires Boost,
    # LZ4, Disruptor).  evio-6.0 has a bug: find_package(Boost REQUIRED)
    # runs unconditionally even with C_ONLY=ON, so we patch it out after
    # download but before add_subdirectory.
    set(C_ONLY ON CACHE BOOL "" FORCE)
    set(MAKE_EXAMPLES OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(evio
        GIT_REPOSITORY https://github.com/JeffersonLab/evio.git
        GIT_TAG        evio-6.0
        GIT_PROGRESS   TRUE
    )
    FetchContent_GetProperties(evio)
    if(NOT evio_POPULATED)
        FetchContent_Populate(evio)
        # Patch out unconditional REQUIRED on Boost (only C++ lib needs it)
        file(READ "${evio_SOURCE_DIR}/CMakeLists.txt" _evio_cml)
        string(REPLACE "find_package(Boost REQUIRED" "find_package(Boost"
               _evio_cml "${_evio_cml}")
        file(WRITE "${evio_SOURCE_DIR}/CMakeLists.txt" "${_evio_cml}")
    endif()
    add_subdirectory(${evio_SOURCE_DIR} ${evio_BINARY_DIR})
elseif(EVIO_SOURCE STREQUAL "system")
    find_package(evio REQUIRED)
endif()

# --- et (optional — only needed for live streaming, not file reading) -----
option(WITH_ET "Build ET (Event Transfer) support for live data" ON)

if(WITH_ET)
    if(ET_SOURCE STREQUAL "hallb")
        find_library(ET_LIB NAMES et PATHS ${CODA_LIB_DIR} NO_DEFAULT_PATH)
        if(ET_LIB)
            add_library(et STATIC IMPORTED GLOBAL)
            set_target_properties(et PROPERTIES
                IMPORTED_LOCATION             "${ET_LIB}"
                INTERFACE_INCLUDE_DIRECTORIES "${CODA_INCLUDE_DIR}"
            )
            set(ET_INCLUDE_DIR "${CODA_INCLUDE_DIR}" CACHE INTERNAL "")
            message(STATUS "et: using Hall-B installation at ${ET_LIB}")
        else()
            message(STATUS "et: Hall-B library not found, falling back to fetch")
            set(ET_SOURCE "fetch")
        endif()
    endif()
    if(ET_SOURCE STREQUAL "fetch")
        FetchContent_Declare(et
            GIT_REPOSITORY https://github.com/JeffersonLab/et.git
            GIT_TAG        main
            GIT_PROGRESS   TRUE
        )
        FetchContent_MakeAvailable(et)
        set(ET_INCLUDE_DIR "${et_SOURCE_DIR}/src/libsrc" CACHE INTERNAL "")
    elseif(ET_SOURCE STREQUAL "system")
        find_package(coda_et REQUIRED)
    endif()
endif()

# ---------------------------------------------------------------------------
# Library target
# ---------------------------------------------------------------------------
set(PRAD2DEC_SOURCES
    src/EvChannel.cpp
    src/Fadc250Decoder.cpp
    src/Fadc250RawDecoder.cpp
    src/Adc1881mDecoder.cpp
    src/SspDecoder.cpp
    src/VtpDecoder.cpp
    src/TdcDecoder.cpp
    src/Dsc2Decoder.cpp
    src/EpicsData.cpp
    src/EpicsStore.cpp
    src/HVDecoder.cpp
    src/WaveAnalyzer.cpp
    src/PulseTemplateStore.cpp
    src/Fadc250FwAnalyzer.cpp
    src/InstallPaths.cpp
)
if(WITH_ET)
    list(APPEND PRAD2DEC_SOURCES src/EtChannel.cpp)
endif()

add_library(prad2dec STATIC ${PRAD2DEC_SOURCES})

# Required so this archive can be linked into shared objects (pybind11
# module, ACLiC-compiled ROOT macros).  Static libs don't need PIC for
# plain executable linking, but it's free on x86-64 and avoids a class of
# "relocation R_X86_64_32S can not be used when making a shared object"
# errors downstream.
set_target_properties(prad2dec PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(prad2dec
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/prad2dec>
)

target_link_libraries(prad2dec PUBLIC evio)

# nlohmann::json — header-only, used by PulseTemplateStore.cpp to load
# the per-channel pulse-template JSON.  load_daq_config.h also pulls it
# in but only via downstream callers, so until now prad2dec didn't need
# the dep itself.
target_link_libraries(prad2dec PRIVATE nlohmann_json::nlohmann_json)

# dladdr is used by InstallPaths.cpp for runtime module-path discovery.
# CMAKE_DL_LIBS is the portable "dl" link name — empty on Windows,
# "-ldl" on old glibc (silently no-op on glibc ≥ 2.34).
if(NOT WIN32)
    target_link_libraries(prad2dec PUBLIC ${CMAKE_DL_LIBS})
endif()

if(WITH_ET)
    target_include_directories(prad2dec PUBLIC $<BUILD_INTERFACE:${ET_INCLUDE_DIR}>)
    target_link_libraries(prad2dec PUBLIC et)
    target_compile_definitions(prad2dec PUBLIC USE_ET=1)
else()
    target_compile_definitions(prad2dec PUBLIC USE_ET=0)
endif()

# ---------------------------------------------------------------------------
# Install
# ---------------------------------------------------------------------------
include(GNUInstallDirs)

install(TARGETS prad2dec ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/prad2dec)

# ---------------------------------------------------------------------------
# Tests (opt-in)
# ---------------------------------------------------------------------------
# Building the firmware-emulation regression test does not pull any external
# dep — just the prad2dec lib (which is the system under test).  Wire it
# into CTest when the parent project enables testing.
option(PRAD2DEC_BUILD_TESTS "Build prad2dec unit tests" ON)
if(PRAD2DEC_BUILD_TESTS)
    enable_testing()
    add_executable(prad2dec_test_fadc_fw test/test_fadc_fw.cpp)
    target_link_libraries(prad2dec_test_fadc_fw PRIVATE prad2dec)
    add_test(NAME prad2dec.fadc_fw COMMAND prad2dec_test_fadc_fw)
endif()
