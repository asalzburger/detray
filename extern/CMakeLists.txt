add_subdirectory(algebra-plugins)

if (DETRAY_UNIT_TESTS)
    set(CMAKE_CXX_STANDARD 17)
    add_subdirectory(googletest)
endif()

if (DETRAY_IO_CSV)
    add_subdirectory(dfelibs)
endif()

if(DETRAY_DISPLAY)
    add_subdirectory(matplotplusplus)
endif()


find_package( vecmem QUIET )

if( NOT vecmem_FOUND )
    message( STATUS "Building vecmem..." )

    include( FetchContent )
    FetchContent_Declare( 
        Vecmem
        GIT_REPOSITORY "https://github.com/acts-project/vecmem.git"
        GIT_TAG "v0.2.0"
    )

    # Do not build google test in vecmem again
    set(BUILD_TESTING OFF CACHE INTERNAL "")

    # Get it into the current directory.
    FetchContent_Populate( Vecmem )
    add_subdirectory( "${vecmem_SOURCE_DIR}" "${vecmem_BINARY_DIR}"
    EXCLUDE_FROM_ALL )
else ()
    find_package( vecmem )
endif()
