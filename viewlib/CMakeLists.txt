PROJECT(viewlib)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

set(CMAKE_MODULE_PATH
	${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

FIND_PACKAGE(OpenMP REQUIRED)

FIND_PACKAGE(JsonCpp REQUIRED)

OPTION(GoTools_USE_Qt4 "Use Qt4 (not Qt5)?" ON)

MESSAGE("HOME: " $ENV{HOME})

IF(GoTools_USE_Qt4)
  # MESSAGE("Using Qt4!")
  FIND_PACKAGE(Qt4 REQUIRED)
#  set(CMAKE_AUTOMOC ON)
#  set(CMAKE_AUTOUIC ON)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
ELSE()
  SET(CMAKE_PREFIX_PATH "C:\\Qt\\Qt5.5.1\\5.5\\msvc2010\\")
  #  SET(QT_VERSION_REQ "5")
  # MESSAGE("Using Qt5!")
  FIND_PACKAGE(Qt5 REQUIRED COMPONENTS Core Widgets Gui OpenGL)
  # Imperative that the QTDIR variable is set (from .bashrc for instance) and pointing to qt5.
  FIND_PROGRAM(QT5_UIC_EXECUTABLE uic
    $ENV{QTDIR}/bin
#    "C:\\Qt\\Qt5.5.1\\5.5\\msvc2010\\"
    )
  IF(POLICY CMP0020)
    cmake_policy(SET CMP0020 OLD)
  ENDIF(POLICY CMP0020)
  IF(POLICY CMP0043)
    cmake_policy(SET CMP0043 OLD)
  ENDIF(POLICY CMP0043)
#  FIND_PACKAGE(Qt5OpenGL REQUIRED)
#  FIND_PACKAGE(Qt5Widgets REQUIRED)
#  set(CMAKE_AUTOMOC ON)
#  set(CMAKE_AUTOUIC ON)
#  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  macro (QT5_EXTRACT_OPTIONS _qt5_files _qt5_options)
    set(${_qt5_files})
    set(${_qt5_options})
    set(_QT5_DOING_OPTIONS FALSE)
    foreach(_currentArg ${ARGN})
      if ("${_currentArg}" STREQUAL "OPTIONS")
        set(_QT5_DOING_OPTIONS TRUE)
      else ()
        if(_QT5_DOING_OPTIONS)
          list(APPEND ${_qt5_options} "${_currentArg}")
        else()
          list(APPEND ${_qt5_files} "${_currentArg}")
        endif()
      endif ()
    endforeach()
  endmacro ()
ENDIF(GoTools_USE_Qt4)
MESSAGE("QT5_UIC_EXECUTABLE: " ${QT5_UIC_EXECUTABLE})
MESSAGE("QT_UIC_EXECUTABLE: " ${QT_UIC_EXECUTABLE})
MESSAGE("QT_MOC_EXECUTABLE: " ${QT_MOC_EXECUTABLE})
MESSAGE("QTDIR: " $ENV{QTDIR})
# Find modules

FIND_PACKAGE(OpenGL REQUIRED)

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

IF(GoTools_USE_Qt4)
  MACRO (QT4_WRAP_UI_GO outfiles outpath)
    IF (${CMAKE_MAJOR_VERSION} GREATER 2)
      QT4_EXTRACT_OPTIONS(ui_files ui_options ui_target ${ARGN})
    ELSEIF (${CMAKE_MAJOR_VERSION} GREATER 1 AND ${CMAKE_MINOR_VERSION} GREATER 7 AND ${CMAKE_PATCH_VERSION} GREATER 11)
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
ELSE()
  # MESSAGE("Inside Qt5 macro!")
  # MESSAGE("QT5_UIC_EXECUTABLE: " ${QT5_UIC_EXECUTABLE})
  # MESSAGE("QTDIR: " $ENV{QTDIR})
  MACRO (QT5_WRAP_UI_GO outfiles outpath)
    IF (${CMAKE_MAJOR_VERSION} GREATER 2)
      QT5_EXTRACT_OPTIONS(ui_files ui_options ${ARGN})
      #QT5_EXTRACT_OPTIONS(ui_files ui_options ui_target ${ARGN})
    ELSEIF (${CMAKE_MAJOR_VERSION} GREATER 1 AND ${CMAKE_MINOR_VERSION} GREATER 7 AND ${CMAKE_PATCH_VERSION} GREATER 11)
      QT5_EXTRACT_OPTIONS(ui_files ui_options ui_target ${ARGN})
    ELSE() 
      QT5_EXTRACT_OPTIONS(ui_files ui_options ${ARGN})
    ENDIF()
    FOREACH (it ${ui_files})
      GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
      GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
      SET(outfile ${viewlib_BINARY_DIR}/${outpath}/ui_${outfile}.h) # Here we set output
      ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
        COMMAND ${QT5_UIC_EXECUTABLE}
        ARGS ${ui_options} -o ${outfile} ${infile}
        MAIN_DEPENDENCY ${infile})
      SET(${outfiles} ${${outfiles}} ${outfile})
    ENDFOREACH (it)
  ENDMACRO (QT5_WRAP_UI_GO)
ENDIF()

OPTION(ViewlibVolAndLRSupport "Include spline volumes & LRSplineSurface?" OFF)
#OPTION(LRSplinesSupport "Include LRSplines?" OFF)


SET(QT_USE_QTOPENGL 1)
#SET(QT_USE_QT3SUPPORT 1)
IF(GoTools_USE_Qt4)
  INCLUDE(${QT_USE_FILE})
ELSE()
  #INCLUDE(${QT_USE_FILE})
ENDIF()

# Include directories


IF(GoTools_USE_Qt4)
  INCLUDE_DIRECTORIES(
    ${viewlib_BINARY_DIR}/include
    ${viewlib_SOURCE_DIR}/include
    ${GoIgeslib_SOURCE_DIR}/include
    ${GoToolsCore_SOURCE_DIR}/include
    ${sisl_SOURCE_DIR}/include
    ${OPENGL_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}
    ${GoTools_COMMON_INCLUDE_DIRS}
    ${QT_INCLUDES}
    #  ${Boost_INCLUDE_DIR}
    )
  # MESSAGE("QT_INCLUDES: " ${QT_INCLUDES})
ELSE()
  INCLUDE_DIRECTORIES(
    ${viewlib_BINARY_DIR}/include
    ${viewlib_SOURCE_DIR}/include
    ${GoIgeslib_SOURCE_DIR}/include
    ${GoToolsCore_SOURCE_DIR}/include
    ${sisl_SOURCE_DIR}/include
    ${OPENGL_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}
    ${GoTools_COMMON_INCLUDE_DIRS}
    #  ${Boost_INCLUDE_DIR}
    ${QT5_INCLUDES}
    ${Qt5Widgets_INCLUDE_DIRS}
    ${Qt5Gui_OPENGL_INCLUDE_DIRS}
    )
  # MESSAGE("QT5_INCLUDES: " ${QT5_INCLUDES})
  # MESSAGE("Qt5Widgets_INCLUDE_DIRS: " ${Qt5Widgets_INCLUDE_DIRS})
  # MESSAGE("Qt5Gui_OPENGL_INCLUDE_DIRS: " ${Qt5Gui_OPENGL_INCLUDE_DIRS})
  # MESSAGE("Qt5Gui_EGL_INCLUDE_DIRS: " ${Qt5Gui_EGL_INCLUDE_DIRS})
ENDIF()

# Linked in libraries

IF(GoTools_USE_Qt4)
  SET(DEPLIBS
    GoIgeslib
    GoToolsCore
    sisl
    ${OPENGL_LIBRARIES}
    ${QT_LIBRARIES}
    )
ELSE()
  SET(DEPLIBS
    GoIgeslib
    GoToolsCore
    sisl
    ${OPENGL_LIBRARIES}
    ${Qt5Widgets_LIBRARIES}
    ${Qt5Gui_OPENGL_LIBRARIES}
    ${Qt5Gui_EGL_LIBRARIES}
#    Qt5OpenGL
    )
  MESSAGE("Qt5Widgets_LIBRARIES: " ${Qt5Widgets_LIBRARIES})
  MESSAGE("Qt5Gui_OPENGL_LIBRARIES: " ${Qt5Gui_OPENGL_LIBRARIES})
  MESSAGE("Qt5Gui_EGL_LIBRARIES: " ${Qt5Gui_EGL_LIBRARIES})
  # MESSAGE("Qt5_OPENGL_LIBRARIES: " ${Qt5_OPENGL_LIBRARIES})
ENDIF()
#MESSAGE("QT_LIBRARIES: " ${QT_LIBRARIES})

IF (NOT WIN32)
  SET(DEPLIBS ${DEPLIBS}
    Qt5OpenGL
    )
ENDIF()

MESSAGE("DEPLIBS: " ${DEPLIBS})

# Qt uic
# We not use GLOB_RECURSE to avoid including volume directory.
# For the file to be generated we need to pass it through the QT4_WRAP_CPP macro.
FILE(GLOB viewlib_UIS include/GoTools/viewlib/*.ui)
IF(GoTools_USE_Qt4)
  QT4_WRAP_UI_GO(viewlib_UIS_H "include/GoTools/viewlib" ${viewlib_UIS})
ELSE()
  # MESSAGE("Calling QT5 UI_GO macro!")
  QT5_WRAP_UI_GO(viewlib_UIS_H "include/GoTools/viewlib" ${viewlib_UIS})
ENDIF()
#MESSAGE(STATUS "viewlib_UIS: " ${viewlib_UIS}) # The source ui-files.
#MESSAGE(STATUS "viewlib_UIS_H: " ${viewlib_UIS_H}) # The generated h-files.

# Qt moc. The option '-nw' is "not recommended" according to the Qt
# documentation, but it gets rid of warnings...
FILE(GLOB viewlib_MOC_HDRS_TMP ${viewlib_SOURCE_DIR}/include/GoTools/viewlib/*.h)
SET(viewlib_MOC_HDRS ${viewlib_MOC_HDRS_TMP} ${viewlib_UIS_H})
#MESSAGE("\nviewlib_MOC_HDRS: " ${viewlib_MOC_HDRS})
IF(GoTools_USE_Qt4)
  QT4_WRAP_CPP(viewlib_MOC_SRCS ${viewlib_MOC_HDRS}
    OPTIONS "-nw")
ELSE()
  # MESSAGE("Calling QT5 CPP macro!")
  QT5_WRAP_CPP(viewlib_MOC_SRCS ${viewlib_MOC_HDRS}
    OPTIONS "-nw")
ENDIF()
#MESSAGE(STATUS "viewlib_MOC_HDRS:" ${viewlib_MOC_HDRS})
#MESSAGE(STATUS "viewlib_MOC_SRCS:" ${viewlib_MOC_SRCS})

# Make the viewlib library

FILE(GLOB viewlib_SRCS src/*.C include/GoTools/viewlib/*.h)
if (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(viewlib SHARED ${viewlib_SRCS} ${viewlib_MOC_SRCS})
else (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(viewlib ${viewlib_SRCS} ${viewlib_MOC_SRCS})
endif (BUILD_AS_SHARED_LIBRARY)

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
  MESSAGE("appname: " ${appname})
  MESSAGE("app: " ${app})
  IF(NOT GoTools_USE_Qt4)
     qt5_use_modules(${appname} Core Widgets Gui OpenGL)
  ENDIF(NOT GoTools_USE_Qt4)
  TARGET_LINK_LIBRARIES(${appname} viewlib ${DEPLIBS})
  SET_TARGET_PROPERTIES(${appname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "viewlib/Apps")
ENDFOREACH(app)

# For volume models, a few more includes and libs are needed.
IF (ViewlibVolAndLRSupport)

  IF(GoTools_USE_Qt4)
    INCLUDE_DIRECTORIES(
      ${GoTrivariate_SOURCE_DIR}/include
      ${GoLRspline2D_SOURCE_DIR}/include
      include/GoTools/viewlib/vol_and_lr/DataHandlerVolAndLR.h
      include/GoTools/viewlib/vol_and_lr/gvRectangularVolumePaintable.h
      include/GoTools/viewlib/vol_and_lr/RectangularVolumePropertySheet.h
      ${JSONCPP_INCLUDE_DIR}
      )
  ELSE()
    INCLUDE_DIRECTORIES(
      ${GoTrivariate_SOURCE_DIR}/include
      ${GoLRspline2D_SOURCE_DIR}/include
      ${JSONCPP_INCLUDE_DIR}
      ${QT5_INCLUDES}
      ${Qt5Widgets_INCLUDE_DIRS}
      ${Qt5Gui_OPENGL_INCLUDE_DIRS}
      )
  ENDIF()
  
  IF(GoTools_USE_Qt4)
    SET(DEPLIBS
      GoIgeslib
      GoToolsCore
      sisl
      GoTrivariate
      GoLRspline2D
      ${OPENGL_LIBRARIES}
      ${JSONCPP_LIBRARY}
      ${QT_LIBRARIES}
      )
  ELSE()
    SET(DEPLIBS
      GoIgeslib
      GoToolsCore
      sisl
      GoTrivariate
      GoLRspline2D
      Qt5OpenGL
      ${OPENGL_LIBRARIES}
      ${JSONCPP_LIBRARY}
      ${Qt5Widgets_LIBRARIES}
      ${Qt5Gui_OPENGL_LIBRARIES}
      )
  ENDIF()
  
  # Qt uic
  FILE(GLOB viewlib_volandlr_UIS ${viewlib_SOURCE_DIR}/include/GoTools/viewlib/vol_and_lr/*.ui)
  IF(GoTools_USE_Qt4)
    QT4_WRAP_UI_GO(viewlib_volandlr_UIS_H "include/GoTools/viewlib/vol_and_lr" ${viewlib_volandlr_UIS})
  ELSE()
    MESSAGE("Vol & LR: Calling QT5 UI_GO macro!")
    QT5_WRAP_UI_GO(viewlib_volandlr_UIS_H "include/GoTools/viewlib/vol_and_lr" ${viewlib_volandlr_UIS})
  ENDIF()
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
  IF(GoTools_USE_Qt4)
    QT4_WRAP_CPP(viewlib_volandlr_MOC_SRCS ${viewlib_volandlr_MOC_HDRS}
      OPTIONS "-nw")
  ELSE()
    MESSAGE("Vol & LR: Calling QT5 CPP macro!")
    QT5_WRAP_CPP(viewlib_volandlr_MOC_SRCS ${viewlib_volandlr_MOC_HDRS}
      OPTIONS "-nw")
    # QT5_WRAP_CPP(viewlib_volandlr_MOC_SRCS ${viewlib_volandlr_MOC_HDRS}
    #   OPTIONS "-nw"))
  ENDIF()
  # MESSAGE(STATUS ${viewlib_volandlr_MOC_HDRS})
  # MESSAGE(STATUS ${viewlib_volandlr_MOC_SRCS})
  # MESSAGE(STATUS "Done with moc stuff!")

  # Make the viewlibVolume library
  FILE(GLOB_RECURSE viewlib_SRCS src/*.C include/*.h)
  if (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(viewlibVolAndLR SHARED ${viewlib_SRCS} ${viewlib_MOC_SRCS} ${viewlib_volandlr_SRCS} ${viewlib_volandlr_MOC_SRCS})
  else (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(viewlibVolAndLR ${viewlib_SRCS} ${viewlib_MOC_SRCS} ${viewlib_volandlr_SRCS} ${viewlib_volandlr_MOC_SRCS})
  endif (BUILD_AS_SHARED_LIBRARY)
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
    MESSAGE("app: " ${app})
    MESSAGE("appname: " ${appname})
    ADD_EXECUTABLE(${appname} ${app})
    IF(NOT GoTools_USE_Qt4)
     qt5_use_modules(${appname} Core Widgets Gui OpenGL)
    ENDIF(NOT GoTools_USE_Qt4)
    TARGET_LINK_LIBRARIES(${appname} viewlibVolAndLR ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "viewlib/Apps")
  ENDFOREACH(app)

ENDIF (ViewlibVolAndLRSupport)

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
