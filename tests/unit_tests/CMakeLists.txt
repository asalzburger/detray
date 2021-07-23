message(STATUS "Testing: 'detray::unit_tests' enabled")

mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

include(GoogleTest)

macro(add_detray_test TESTNAME FILES PLUGIN_LIBRARY)
    add_executable(${TESTNAME} ${FILES})
    target_link_libraries(${TESTNAME} PRIVATE gtest gmock gtest_main)
    target_link_libraries(${TESTNAME} PRIVATE detray::core detray::io detray::tests_common vecmem::core)
    target_link_libraries(${TESTNAME} PRIVATE  ${PLUGIN_LIBRARY} ${UNIT_TEST_EXTRA_LIBRARIES})
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
    add_test(unit_test_${TESTNAME} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${TESTNAME})
    set_tests_properties(unit_test_${TESTNAME} PROPERTIES 
        ENVIRONMENT DETRAY_TEST_DATA_DIR=${PROJECT_SOURCE_DIR}/data/)
endmacro()

set(all_unit_tests "core")
list(APPEND all_unit_tests "annulus2;cylinder3;rectangle2;ring2;single3;unmasked;trapezoid2")
list(APPEND all_unit_tests "propagator;line_stepper;navigator;cylinder_intersection;planar_intersection;intersection_kernel")
list(APPEND all_unit_tests "detector;transform_store;masks_container")

add_subdirectory(core)
if(DETRAY_ARRAY_PLUGIN)
    add_subdirectory(array)
endif()
if(DETRAY_EIGEN_PLUGIN)
    add_subdirectory(eigen)
endif()
if(DETRAY_SMATRIX_PLUGIN)
    add_subdirectory(smatrix)
endif()
if(DETRAY_VC_ALGEBRA_PLUGIN)
    add_subdirectory(vc)
endif()
