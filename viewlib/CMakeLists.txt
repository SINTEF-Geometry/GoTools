PROJECT(viewlib)


# Find modules

FIND_PACKAGE(OpenGL REQUIRED)
FIND_PACKAGE(Qt4 REQUIRED)
SET(Boost_ADDITIONAL_VERSIONS 
	"1.42" "1.42.0" 
	"1.43" "1.43.0" 
	"1.44" "1.44.0"
	"1.46" "1.46.0"
	"1.47" "1.47.0")
FIND_PACKAGE(Boost REQUIRED)


MACRO (QT4_WRAP_UI_GO outfiles )
  QT4_EXTRACT_OPTIONS(ui_files ui_options ${ARGN})
  FOREACH (it ${ui_files})
    GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
    GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
    SET(outfile ${CMAKE_CURRENT_BINARY_DIR}/include/GoTools/viewlib/ui_${outfile}.h) # Here we set output
    ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
      COMMAND ${QT_UIC_EXECUTABLE}
      ARGS ${ui_options} -o ${outfile} ${infile}
      MAIN_DEPENDENCY ${infile})
    SET(${outfiles} ${${outfiles}} ${outfile})
  ENDFOREACH (it)
ENDMACRO (QT4_WRAP_UI_GO)

SET(QT_USE_QTOPENGL 1)
SET(QT_USE_QT3SUPPORT 1)
INCLUDE(${QT_USE_FILE})


# Include directories

INCLUDE_DIRECTORIES(
  ${viewlib_BINARY_DIR}/include
  ${viewlib_SOURCE_DIR}/include
  ${GoIgeslib_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${GoTools_COMMON_INCLUDE_DIRS}
  ${QT_INCLUDES}
  ${OPENGL_INCLUDE_DIR}
  )


# Linked in libraries

SET(DEPLIBS
  GoIgeslib
  GoToolsCore
  sisl
  ${QT_LIBRARIES}
  ${OPENGL_LIBRARIES}
  )


# Qt uic
FILE(GLOB_RECURSE viewlib_UIS include/*.ui)
QT4_WRAP_UI_GO(viewlib_UIS_H ${viewlib_UIS})


# Qt moc. The option '-nw' is "not recommended" according to the Qt
# documentation, but it gets rid of warnings...
FILE(GLOB_RECURSE viewlib_MOC_HDRS_TMP include/*.h)
SET(viewlib_MOC_HDRS ${viewlib_MOC_HDRS_TMP} ${viewlib_UIS_H})
QT4_WRAP_CPP(viewlib_MOC_SRCS ${viewlib_MOC_HDRS}
  OPTIONS "-nw")


# Make the viewlib library

FILE(GLOB_RECURSE viewlib_SRCS src/*.C include/*.h)
ADD_LIBRARY(viewlib ${viewlib_SRCS} ${viewlib_MOC_SRCS})
SET_PROPERTY(TARGET viewlib
  PROPERTY FOLDER "viewlib/Libs")


# Set various compiler flags

IF(NOT WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-conversion")
ENDIF(NOT WIN32)


# Apps, examples, tests, ...?

FILE(GLOB viewlib_APPS app/*.C)
FOREACH(app ${viewlib_APPS})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} viewlib ${DEPLIBS})
  SET_TARGET_PROPERTIES(${appname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "viewlib/Apps")
ENDFOREACH(app)

# Copy data
ADD_CUSTOM_COMMAND(
  TARGET viewlib
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy_directory ${viewlib_SOURCE_DIR}/data
  ${viewlib_BINARY_DIR}/data
  )


# 'install' target
IF(WIN32)
  # Windows
  # lib
  INSTALL(TARGETS viewlib DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  INSTALL(DIRECTORY include/GoTools/viewlib
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
  # include ui files
  INSTALL(DIRECTORY ${viewlib_BINARY_DIR}/include/GoTools/viewlib/
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools/viewlib
    FILES_MATCHING PATTERN "ui_*_form.h"
    PATTERN ".svn" EXCLUDE
    )
  # data
  INSTALL(DIRECTORY data/
    DESTINATION ${GoTools_INSTALL_PREFIX}/share/viewlib
    FILES_MATCHING PATTERN "*.g2" PATTERN "*.igs"
    PATTERN ".svn" EXCLUDE
    )
  # bin
  SET(goview_INSTALL_PREFIX "$ENV{PROGRAMFILES}/SINTEF/goview"
    CACHE PATH "Path to install goview")
  INSTALL(TARGETS goview DESTINATION ${goview_INSTALL_PREFIX})
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS viewlib DESTINATION lib COMPONENT viewlib)
  # include
  INSTALL(DIRECTORY include/GoTools/viewlib
    DESTINATION include/GoTools
    COMPONENT viewlib-dev
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
  # include ui files
  INSTALL(DIRECTORY ${viewlib_BINARY_DIR}/include/GoTools/viewlib/
    DESTINATION include/GoTools/viewlib
    COMPONENT viewlib-dev
    FILES_MATCHING PATTERN "ui_*_form.h"
    PATTERN ".svn" EXCLUDE
    )
  # data
  INSTALL(DIRECTORY data/
    DESTINATION share/GoTools/viewlib
    COMPONENT viewlib
    FILES_MATCHING PATTERN "*.g2" PATTERN "*.igs"
    PATTERN ".svn" EXCLUDE
    )
  # bin
  INSTALL(TARGETS goview DESTINATION bin COMPONENT viewlib)
ENDIF(WIN32)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libviewlib.so)
