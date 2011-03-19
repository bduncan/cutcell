project( cutcell )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif()

find_package(CGAL QUIET COMPONENTS Core )
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

    find_package(Qt4 COMPONENTS QtCore QtGui REQUIRED)
    if (QT4_FOUND)
        include(${QT_USE_FILE})
        set(CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed -Wl,--no-undefined ${CMAKE_EXE_LINKER_FLAGS}")
        include_directories(BEFORE ${QT_INCLUDE_DIR})
        QT4_WRAP_UI(cutcelluih "cutcell.ui")
        QT4_WRAP_CPP(cutcellguimoc "include/cutcellgui.hpp")
        include_directories(${CMAKE_CURRENT_BINARY_DIR})
        add_executable("cutcellgui" ${cutcellguimoc} ${cutcelluih} "cutcellgui.cpp" "cutcell.cpp")
        add_to_cached_list(CGAL_EXECUTABLE_TARGETS "cutcellgui")

        # Link the executable to CGAL and third-party libraries
        target_link_libraries("cutcellgui" ${CGAL_LIBRARIES} ${QT_LIBRARIES} cgns)
        install( TARGETS cutcellgui RUNTIME DESTINATION bin )
    else (QT4_FOUND)
        message(STATUS "NOTICE: This program requires Qt4, and will not be compiled.")
    endif (QT4_FOUND)

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
