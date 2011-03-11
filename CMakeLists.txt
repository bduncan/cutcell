project( cutcell )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif()

find_package(CGAL QUIET COMPONENTS Core Qt3 )
include(${CGAL_USE_FILE})

if ( CGAL_FOUND )

    set(CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed -Wl,--no-undefined ${CMAKE_EXE_LINKER_FLAGS}")
    include_directories(BEFORE include)
    add_executable( "make_cube_3" "make_cube_3.cpp")
    add_to_cached_list(CGAL_EXECUTABLE_TARGETS "make_cube_3")
    # Link the executable to CGAL and third-party libraries
    if ( CGAL_AUTO_LINK_ENABLED )
      target_link_libraries("make_cube_3" ${CGAL_3RD_PARTY_LIBRARIES})
    else()
      target_link_libraries("make_cube_3" ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
    endif()

    set(CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed -Wl,--no-undefined ${CMAKE_EXE_LINKER_FLAGS}")
    add_executable( "cutcelltest" "cutcelltest.cpp" "cutcell.cpp")
    add_to_cached_list(CGAL_EXECUTABLE_TARGETS "cutcelltest")
    # Link the executable to CGAL and third-party libraries
    if ( CGAL_AUTO_LINK_ENABLED )
      target_link_libraries("cutcelltest" ${CGAL_3RD_PARTY_LIBRARIES} cgns)
    else()
      target_link_libraries("cutcelltest" ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} cgns)
    endif()
    install( TARGETS cutcelltest RUNTIME DESTINATION bin )

    find_package(Qt3-patched)
    if (QT3_FOUND AND CGAL_Qt3_FOUND)
        set(CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed -Wl,--no-undefined ${CMAKE_EXE_LINKER_FLAGS}")
        include_directories(BEFORE ${QT3_INCLUDE_DIR})
        add_executable("visualization_SNC" "visualization_SNC.cpp")
        add_to_cached_list(CGAL_EXECUTABLE_TARGETS "visualization_SNC")
        # Link the executable to CGAL and third-party libraries
        if ( CGAL_AUTO_LINK_ENABLED )
          target_link_libraries("visualization_SNC" ${CGAL_3RD_PARTY_LIBRARIES} ${QT3_LIBRARIES} )
        else()
          target_link_libraries("visualization_SNC" ${CGAL_LIBRARIES} ${CGAL_QT_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT3_LIBRARIES})
        endif()

    else (QT3_FOUND AND CGAL_Qt3_FOUND)

      message(STATUS "NOTICE: This demo requires CGAL and Qt3, and will not be compiled.")

    endif (QT3_FOUND AND CGAL_Qt3_FOUND)

    find_package(Boost COMPONENTS program_options)
    if (Boost_FOUND AND Boost_PROGRAM_OPTIONS_FOUND)
        set(CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed -Wl,--no-undefined ${CMAKE_EXE_LINKER_FLAGS}")
        include_directories(${Boost_INCLUDE_DIRS})
        add_executable("cutcell" "cutcellcli.cpp" "cutcell.cpp")
        target_link_libraries("cutcell" ${Boost_LIBRARIES} ${Boost_program_options_LIBRARY})
        # Link the executable to CGAL and third-party libraries
        if ( CGAL_AUTO_LINK_ENABLED )
          target_link_libraries("cutcell" ${CGAL_3RD_PARTY_LIBRARIES} cgns)
        else()
          target_link_libraries("cutcell" ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} cgns)
        endif()
        install( TARGETS cutcell RUNTIME DESTINATION bin )
    else (Boost_FOUND AND Boost_PROGRAM_OPTIONS_FOUND)
        message(STATUS "NOTICE: This program requires Boost, and will not be compiled.")
    endif(Boost_FOUND AND Boost_PROGRAM_OPTIONS_FOUND)

endif()
