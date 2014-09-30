PROJECT(viewlib)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Find modules

FIND_PACKAGE(OpenGL REQUIRED)
FIND_PACKAGE(Qt4 REQUIRED)

# SET(Boost_ADDITIONAL_VERSIONS 
# 	"1.42" "1.42.0" 
# 	"1.43" "1.43.0" 
# 	"1.44" "1.44.0"
# 	"1.47" "1.47.0")

IF(WIN32)
  SET(GLUT_INCLUDE_DIR
    "$ENV{PROGRAMW6432}/Microsoft SDKs/Windows/v6.0A/Include"
    )
ELSE(WIN32)
  FIND_PACKAGE(GLUT REQUIRED)
ENDIF(WIN32)


MACRO (QT4_WRAP_UI_GO outfiles outpath)
  IF (${CMAKE_MAJOR_VERSION} GREATER 1 AND ${CMAKE_MINOR_VERSION} GREATER 7 AND ${CMAKE_PATCH_VERSION} GREATER 11)
    QT4_EXTRACT_OPTIONS(ui_files ui_options ui_target ${ARGN})
  ELSE() 
    QT4_EXTRACT_OPTIONS(ui_files ui_options ${ARGN})
  ENDIF()
  FOREACH (it ${ui_files})
    GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
    GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
    SET(outfile ${viewlib_BINARY_DIR}/${outpath}/ui_${outfile}.h) # Here we set output
    ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
      COMMAND ${QT_UIC_EXECUTABLE}
      ARGS ${ui_options} -o ${outfile} ${infile}
      MAIN_DEPENDENCY ${infile})
    SET(${outfiles} ${${outfiles}} ${outfile})
  ENDFOREACH (it)
ENDMACRO (QT4_WRAP_UI_GO)


OPTION(ViewlibVolAndLRSupport "Include spline volumes & LRSplineSurface?" OFF)
#OPTION(LRSplinesSupport "Include LRSplines?" OFF)


SET(QT_USE_QTOPENGL 1)
SET(QT_USE_QT3SUPPORT 1)
INCLUDE(${QT_USE_FILE})


# Include directories

INCLUDE_DIRECTORIES(
  ${viewlib_BINARY_DIR}/include
  ${viewlib_SOURCE_DIR}/include
  ${GoIgeslib_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${sisl_SOURCE_DIR}/include
  ${QT_INCLUDES}
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${GoTools_COMMON_INCLUDE_DIRS}
#  ${Boost_INCLUDE_DIR}
  )


# Linked in libraries

SET(DEPLIBS
  GoIgeslib
  GoToolsCore
  sisl
  ${QT_LIBRARIES}
  ${OPENGL_LIBRARIES}
  )


# # We also support LRSplines.
# IF (LRSplinesSupport)
# 	INCLUDE_DIRECTORIES(
# 		${CUDA_TOOLKIT_INCLUDE}
# 		${CUSPARSE_INCLUDE})
# 	CUDA_INCLUDE_DIRECTORIES(
# 		${LinSolver_SOURCE_DIR}/src
# 		${LinSolver_SOURCE_DIR}/include)
# ENDIF (LRSplinesSupport)


# Qt uic
# We not use GLOB_RECURSE to avoid including volume directory.
# For the file to be generated we need to pass it through the QT4_WRAP_CPP macro.
FILE(GLOB viewlib_UIS include/GoTools/viewlib/*.ui)
QT4_WRAP_UI_GO(viewlib_UIS_H "include/GoTools/viewlib" ${viewlib_UIS})
#MESSAGE(STATUS "viewlib_UIS: " ${viewlib_UIS}) # The source ui-files.
#MESSAGE(STATUS "viewlib_UIS_H: " ${viewlib_UIS_H}) # The generated h-files.

# Qt moc. The option '-nw' is "not recommended" according to the Qt
# documentation, but it gets rid of warnings...
FILE(GLOB viewlib_MOC_HDRS_TMP ${viewlib_SOURCE_DIR}/include/GoTools/viewlib/*.h)
SET(viewlib_MOC_HDRS ${viewlib_MOC_HDRS_TMP} ${viewlib_UIS_H})
QT4_WRAP_CPP(viewlib_MOC_SRCS ${viewlib_MOC_HDRS}
  OPTIONS "-nw")
# MESSAGE(STATUS "viewlib_MOC_HDRS:"${viewlib_MOC_HDRS})
# MESSAGE(STATUS "viewlib_MOC_SRCS:"${viewlib_MOC_SRCS})

# Make the viewlib library

FILE(GLOB viewlib_SRCS src/*.C include/GoTools/viewlib/*.h)
ADD_LIBRARY(viewlib ${viewlib_SRCS} ${viewlib_MOC_SRCS})
SET_PROPERTY(TARGET viewlib
  PROPERTY FOLDER "viewlib/Libs")


# For volume models, a few more includes and libs are needed.
IF (ViewlibVolAndLRSupport)
  INCLUDE_DIRECTORIES(
    ${GoTrivariate_SOURCE_DIR}/include
    ${GoLRspline2D_SOURCE_DIR}/include
    include/GoTools/viewlib/vol_and_lr/DataHandlerVolAndLR.h
    include/GoTools/viewlib/vol_and_lr/gvRectangularVolumePaintable.h
    include/GoTools/viewlib/vol_and_lr/RectangularVolumePropertySheet.h
    )

  SET(DEPLIBS
    GoIgeslib
    GoToolsCore
    sisl
    GoTrivariate
    GoLRspline2D
    ${QT_LIBRARIES}
    ${OPENGL_LIBRARIES}
    )

  # Qt uic
  FILE(GLOB viewlib_volandlr_UIS ${viewlib_SOURCE_DIR}/include/GoTools/viewlib/vol_and_lr/*.ui)
  QT4_WRAP_UI_GO(viewlib_volandlr_UIS_H "include/GoTools/viewlib/vol_and_lr" ${viewlib_volandlr_UIS})
  #  QT4_WRAP_UI_GO_VOLUME(viewlib_volandlr_UIS_H ${viewlib_volandlr_UIS})
  # MESSAGE(STATUS "Wrapping volume!")
  # MESSAGE(STATUS "Done wrapping volume!")
  # MESSAGE(STATUS "viewlib_volandlr_UIS: "${viewlib_volandlr_UIS}) # The source ui-files.
  # MESSAGE(STATUS "viewlib_volandlr_UIS_H: "${viewlib_volandlr_UIS_H}) # The generated h-files.

  # Qt moc. The option '-nw' is "not recommended" according to the Qt
  # documentation, but it gets rid of warnings...
  # MESSAGE(STATUS "Starting with moc stuff!")
  FILE(GLOB viewlib_volandlr_MOC_HDRS_TMP ${viewlib_SOURCE_DIR}/include/GoTools/viewlib/vol_and_lr/*.h)
  # MESSAGE(STATUS "viewlib_volandlr_MOC_HDRS_TMP:"${viewlib_volandlr_MOC_HDRS_TMP})
  SET(viewlib_volandlr_MOC_HDRS ${viewlib_volandlr_MOC_HDRS_TMP} ${viewlib_volandlr_UIS_H})
  QT4_WRAP_CPP(viewlib_volandlr_MOC_SRCS ${viewlib_volandlr_MOC_HDRS}
    OPTIONS "-nw")
  # MESSAGE(STATUS ${viewlib_volandlr_MOC_HDRS})
  # MESSAGE(STATUS ${viewlib_volandlr_MOC_SRCS})
  # MESSAGE(STATUS "Done with moc stuff!")

  # Make the viewlibVolume library
  FILE(GLOB_RECURSE viewlib_SRCS src/*.C include/*.h)
  ADD_LIBRARY(viewlibVolAndLR ${viewlib_SRCS} ${viewlib_MOC_SRCS} ${viewlib_volandlr_SRCS} ${viewlib_volandlr_MOC_SRCS})
  SET_PROPERTY(TARGET viewlibVolAndLR
    PROPERTY FOLDER "viewlib/Libs")

  # # Make the viewlibVolume library.
  # FILE(GLOB viewlibVolume_src ${viewlib_SOURCE_DIR}/src/volume/*.C)
  # ADD_LIBRARY(viewlibVolume STATIC ${viewlibVolume_src})
  # MESSAGE(STATUS viewlibVolume_src: ${viewlibVolume_src})
  # # TARGET_LINK_LIBRARIES(app_BenchSolvers 
  # #   LinSolverCBlasDep  ${CBLAS_LIBRARY} ${ATLAS_LIBRARY})

  FILE(GLOB viewlibVolAndLR_APPS app/vol_and_lr/*.C)
  FOREACH(app ${viewlibVolAndLR_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    # MESSAGE(STATUS "app: "${app})
    # MESSAGE(STATUS "appname: "${appname})
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} viewlibVolAndLR ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "viewlib/Apps")
  ENDFOREACH(app)

ENDIF (ViewlibVolAndLRSupport)


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

#ADD_SUBDIRECTORY(src)

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
