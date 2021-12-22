project(viewlib)

set(CMAKE_MODULE_PATH
	${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

if(GoTools_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif(GoTools_ENABLE_OPENMP)

find_package(OpenGL REQUIRED)
find_path(OPENGL_INCLUDE_DIR
  gl.h
  PATH_PREFIXES GL
  PATHS
  /usr/local/include
  /usr/include
  "$ENV{HOME}\\Install\\include" # For some reason HOME is not found in VS 2019.
  "C:\\Users\\$ENV{USERNAME}\\Install\\include")
#  "C:\\Users\\sbr\\Install\\include")
option(GoTools_USE_Qt4 "Use Qt4 (not Qt5)?" OFF)

if(GoTools_USE_Qt4)
  find_package(Qt4 REQUIRED)
#  set(CMAKE_AUTOMOC ON)
#  set(CMAKE_AUTOUIC ON)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
else()
  if(WIN32 OR APPLE)
    # Imperative that the QTDIR variable is set (from .bashrc for instance) and pointing to qt5.
    if (QTDIR OR DEFINED ENV{QTDIR})
      # The variable was already set.
    else()
      set(QTDIR "" CACHE PATH "Path to Qt (e.g. C:/Qt/Qt5.6.2/5.6/msvc2015_64/)")
      message(SEND_ERROR "QTDIR variable is missing but required. Set to QT install dir with the required lib dir.")
    endif()
    set(CMAKE_PREFIX_PATH "$ENV{QTDIR}")
  endif()
  find_package(Qt5 REQUIRED COMPONENTS Core Widgets Gui OpenGL)
  # set(QT_VERSION_REQ "5")
  if(WIN32 OR APPLE)
    find_program(QT5_UIC_EXECUTABLE bin/uic
      	         PATHS "${QTDIR}")
 #   message("QT5_UIC_EXECUTABLE: " ${QT5_UIC_EXECUTABLE})
else()
     find_program(QT5_UIC_EXECUTABLE uic
      PATHS /usr/lib/x86_64-linux-gnu/qt5/bin)
    # $ENV{QTDIR}/bin
#    PATH_SUFFIXES qt5/bin/)
  # NO_SYSTEM_ENVIRONMENT_PATH)
  # NO_DEFAULT_PATH)
    #    "C:\\Qt\\Qt5.5.1\\5.5\\msvc2010\\"
    #)
  endif()
  #message("QT5_UIC_EXECUTABLE: " ${QT5_UIC_EXECUTABLE})
  #message("QT_UIC_EXECUTABLE: " ${QT_UIC_EXECUTABLE})

  if(POLICY CMP0020)
    cmake_policy(SET CMP0020 OLD)
  endif(POLICY CMP0020)
  if(POLICY CMP0043)
    cmake_policy(SET CMP0043 OLD)
  endif(POLICY CMP0043)
#  set(CMAKE_AUTOUIC ON)
#  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  macro(QT5_EXTRACT_OPTIONS _qt5_files _qt5_options)
    set(${_qt5_files})
    set(${_qt5_options})
    set(_QT5_DOING_OPTIONS FALSE)
    foreach(_currentArg ${ARGN})
      if("${_currentArg}" STREQUAL "OPTIONS")
        set(_QT5_DOING_OPTIONS TRUE)
      else()
        if(_QT5_DOING_OPTIONS)
          list(APPEND ${_qt5_options} "${_currentArg}")
        else()
          list(APPEND ${_qt5_files} "${_currentArg}")
        endif()
      endif()
    endforeach()
  endmacro()
endif()
# message("QT5_UIC_EXECUTABLE: " ${QT5_UIC_EXECUTABLE})
# message("QT_UIC_EXECUTABLE: " ${QT_UIC_EXECUTABLE})
# message("QT_MOC_EXECUTABLE: " ${QT_MOC_EXECUTABLE})
# Find modules

if(WIN32)
  set(GLUT_INCLUDE_DIR
    "$ENV{PROGRAMW6432}/Microsoft SDKs/Windows/v6.0A/Include"
    )
else()
  find_package(GLUT REQUIRED)
endif()

if(GoTools_USE_Qt4)
  macro(QT4_WRAP_UI_GO outfiles outpath)
    if(${CMAKE_MAJOR_VERSION} GREATER 2)
      qt4_extract_options(ui_files ui_options ui_target ${ARGN})
    elseif(${CMAKE_MAJOR_VERSION} GREATER 1 AND ${CMAKE_MINOR_VERSION} GREATER 7 AND ${CMAKE_PATCH_VERSION} GREATER 11)
      qt4_extract_options(ui_files ui_options ui_target ${ARGN})
    else() 
      qt4_extract_options(ui_files ui_options ${ARGN})
    endif()
    foreach(it ${ui_files})
      get_filename_component(outfile ${it} NAME_WE)
      get_filename_component(infile ${it} ABSOLUTE)
      set(outfile ${viewlib_BINARY_DIR}/${outpath}/ui_${outfile}.h) # Here we set output
      add_custom_command(OUTPUT ${outfile}
        COMMAND ${QT_UIC_EXECUTABLE}
        ARGS ${ui_options} -o ${outfile} ${infile}
        MAIN_DEPENDENCY ${infile})
      set(${outfiles} ${${outfiles}} ${outfile})
    endforeach(it)
  endmacro(QT4_WRAP_UI_GO)
else()
  # message("Inside Qt5 macro!")
  # message("QT5_UIC_EXECUTABLE: " ${QT5_UIC_EXECUTABLE})
  # message("QTDIR: " $ENV{QTDIR})
  macro(QT5_WRAP_UI_GO outfiles outpath)
    if(${CMAKE_MAJOR_VERSION} GREATER 2)
      qt5_extract_options(ui_files ui_options ${ARGN})
      #qt5_extract_options(ui_files ui_options ui_target ${ARGN})
    elseif(${CMAKE_MAJOR_VERSION} GREATER 1 AND ${CMAKE_MINOR_VERSION} GREATER 7 AND ${CMAKE_PATCH_VERSION} GREATER 11)
      qt5_extract_options(ui_files ui_options ${ARGN})
    else() 
      qt5_extract_options(ui_files ui_options ${ARGN})
    endif()
    foreach(it ${ui_files})
      get_filename_component(outfile ${it} NAME_WE)
      get_filename_component(infile ${it} ABSOLUTE)
      set(outfile ${viewlib_BINARY_DIR}/${outpath}/ui_${outfile}.h) # Here we set output
      add_custom_command(OUTPUT ${outfile}
        COMMAND ${QT5_UIC_EXECUTABLE}
        ARGS ${ui_options} -o ${outfile} ${infile}
        MAIN_DEPENDENCY ${infile})
      set(${outfiles} ${${outfiles}} ${outfile})
    endforeach(it)
  endmacro(QT5_WRAP_UI_GO)
endif()

option(ViewlibVolAndLRSupport "Include spline volumes & LRSplineSurface?" ON)
#option(LRSplinesSupport "Include LRSplines?" OFF)


set(QT_USE_QTOPENGL 1)
#set(QT_USE_QT3SUPPORT 1)
if(GoTools_USE_Qt4)
  include(${QT_USE_FILE})
else()
  #include(${QT_USE_FILE})
endif()

# Include directories

include_directories(
  ${viewlib_BINARY_DIR}/include
  ${viewlib_SOURCE_DIR}/include
  ${GoIgeslib_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${sisl_SOURCE_DIR}/include
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${GoTools_COMMON_INCLUDE_DIRS}
  #  ${Boost_INCLUDE_DIR}
  )
if(GoTools_USE_Qt4)
  include_directories(
    ${QT_INCLUDES}
    )
  # message("QT_INCLUDES: " ${QT_INCLUDES})
else()
  include_directories(
    ${QT5_INCLUDES}
    ${Qt5Widgets_INCLUDE_DIRS}
    ${Qt5Gui_OPENGL_INCLUDE_DIRS}
    )
  # message("QT5_INCLUDES: " ${QT5_INCLUDES})
  # message("Qt5Widgets_INCLUDE_DIRS: " ${Qt5Widgets_INCLUDE_DIRS})
  # message("Qt5Gui_OPENGL_INCLUDE_DIRS: " ${Qt5Gui_OPENGL_INCLUDE_DIRS})
  # message("Qt5Gui_EGL_INCLUDE_DIRS: " ${Qt5Gui_EGL_INCLUDE_DIRS})
endif()

# Linked in libraries

set(DEPLIBS
  GoIgeslib
  GoToolsCore
  sisl
  ${OPENGL_LIBRARIES}
)

#message("OPENGL_LIBRARIES" ${OPENGL_LIBRARIES})
  
#message("QT_LIBRARIES: " ${QT_LIBRARIES})
if(GoTools_USE_Qt4)
  set(DEPLIBS ${DEPLIBS}
    ${QT_LIBRARIES}
    )
else()
    #message("OPENGL_LIBRARIES: " ${OPENGL_LIBRARIES})
  if(NOT (WIN32 OR APPLE))
    set(DEPLIBS ${DEPLIBS}
      Qt5OpenGL
      )
  endif()
  if(NOT APPLE)
    set(DEPLIBS ${DEPLIBS}
      ${Qt5Widgets_LIBRARIES}
      ${Qt5Core_LIBRARIES}
      ${Qt5Gui_OPENGL_LIBRARIES}
      ${Qt5Gui_EGL_LIBRARIES}
      )
  endif()
 # message("Qt5Widgets_LIBRARIES: " ${Qt5Widgets_LIBRARIES})
 # message("Qt5Gui_OPENGL_LIBRARIES: " ${Qt5Gui_OPENGL_LIBRARIES})
 # message("Qt5Gui_EGL_LIBRARIES: " ${Qt5Gui_EGL_LIBRARIES})
 # message("Qt5_OPENGL_LIBRARIES: " ${Qt5_OPENGL_LIBRARIES})
endif()

#message("DEPLIBS: " ${DEPLIBS})

# Qt uic
# We not use GLOB_RECURSE to avoid including volume directory.
# For the file to be generated we need to pass it through the QT4_WRAP_CPP macro.
file(GLOB viewlib_UIS include/GoTools/viewlib/*.ui)
if(GoTools_USE_Qt4)
  qt4_wrap_ui_go(viewlib_UIS_H "include/GoTools/viewlib" ${viewlib_UIS})
else()
  # message("Calling QT5 UI_GO macro!")
  qt5_wrap_ui_go(viewlib_UIS_H "include/GoTools/viewlib" ${viewlib_UIS})
endif()
#message(STATUS "viewlib_UIS: " ${viewlib_UIS}) # The source ui-files.
#message(STATUS "viewlib_UIS_H: " ${viewlib_UIS_H}) # The generated h-files.

# Qt moc. The option '-nw' is "not recommended" according to the Qt
# documentation, but it gets rid of warnings...
file(GLOB viewlib_MOC_HDRS_TMP ${viewlib_SOURCE_DIR}/include/GoTools/viewlib/*.h)
set(viewlib_MOC_HDRS ${viewlib_MOC_HDRS_TMP} ${viewlib_UIS_H})
#message("\nviewlib_MOC_HDRS: " ${viewlib_MOC_HDRS})
if(GoTools_USE_Qt4)
  qt4_wrap_cpp(viewlib_MOC_SRCS ${viewlib_MOC_HDRS}
    OPTIONS "-nw")
else()
  # message("Calling QT5 CPP macro!")
  qt5_wrap_cpp(viewlib_MOC_SRCS ${viewlib_MOC_HDRS}
    OPTIONS "-nw")
endif()
#message(STATUS "viewlib_MOC_HDRS:" ${viewlib_MOC_HDRS})
#message(STATUS "viewlib_MOC_SRCS:" ${viewlib_MOC_SRCS})

# Make the viewlib library

file(GLOB viewlib_SRCS src/*.C include/GoTools/viewlib/*.h)
if(BUILD_AS_SHARED_LIBRARY)
    add_library(viewlib SHARED ${viewlib_SRCS} ${viewlib_MOC_SRCS})
else(BUILD_AS_SHARED_LIBRARY)
    add_library(viewlib ${viewlib_SRCS} ${viewlib_MOC_SRCS})
endif(BUILD_AS_SHARED_LIBRARY)

set_property(TARGET viewlib
  PROPERTY FOLDER "viewlib/Libs")

# Set various compiler flags

if(NOT WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-conversion")
endif()

# Apps, examples, tests, ...?

file(GLOB viewlib_APPS app/*.C)
foreach(app ${viewlib_APPS})
  get_filename_component(appname ${app} NAME_WE)
  add_executable(${appname} ${app})
  # message("appname: " ${appname})
  # message("app: " ${app})
  if(NOT GoTools_USE_Qt4)
    qt5_use_modules(${appname} Core Widgets Gui OpenGL)
  endif()
  # message("DEPLIBS app linker: " ${DEPLIBS})
  # message("QT_LIBRARIES: " ${QT_LIBRARIES})

  target_link_libraries(${appname} viewlib ${DEPLIBS})
  set_target_properties(${appname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    set_property(TARGET ${appname}
      PROPERTY FOLDER "viewlib/Apps")
endforeach(app)

# Copy Qt dll's.
IF(WIN32)
  if(GoTools_USE_Qt4)
    set(QT_DLLS QtWidgets4.dll QtGui4.dll QtCore4.dll QtOpenGL4.dll opengl32sw.dll) # QtCore4.dll for Qt4 I think.
    set(QT_DEBUG_DLLS QtWidgets4d.dll QtGui4d.dll QtCore4d.dll QtOpenGL4d.dll opengl32sw.dll) # QtCore4d.dll for Qt4 I think.
  else()
    set(QT_DLLS Qt5Widgets.dll Qt5Gui.dll Qt5Core.dll Qt5OpenGL.dll opengl32sw.dll)
    set(QT_DEBUG_DLLS Qt5Widgetsd.dll Qt5Guid.dll Qt5Cored.dll Qt5OpenGLd.dll opengl32sw.dll)
  endif()
  # First the release dll's.
  file(MAKE_DIRECTORY ${viewlib_BINARY_DIR}\\app\\Release) # We make sure that the directory exists.
  foreach( file_i ${QT_DLLS})
    ADD_CUSTOM_COMMAND(
      TARGET viewlib
      POST_BUILD
      COMMAND ${CMAKE_COMMAND}
      ARGS -E copy $ENV{QTDIR}\\bin\\${file_i} ${viewlib_BINARY_DIR}\\app\\Release\\
    )
  endforeach(file_i)
  # Then the debug dll's.
    file(MAKE_DIRECTORY ${viewlib_BINARY_DIR}\\app\\Debug) # We make sure that the directory exists.
  foreach( file_i ${QT_DEBUG_DLLS})
    ADD_CUSTOM_COMMAND(
      TARGET viewlib
      POST_BUILD
      COMMAND ${CMAKE_COMMAND}
      ARGS -E copy $ENV{QTDIR}\\bin\\${file_i} ${viewlib_BINARY_DIR}\\app\\Debug\\
    )
  endforeach(file_i)
ENDIF(WIN32)

# For volume models, a few more includes and libs are needed.
if(ViewlibVolAndLRSupport)

  find_package(JsonCpp REQUIRED)

  include_directories(
    ${GoTrivariate_SOURCE_DIR}/include
    ${GoLRspline2D_SOURCE_DIR}/include
    ${JSONCPP_INCLUDE_DIR}
    )
  
  set(DEPLIBS ${DEPLIBS}
    GoTrivariate
    GoLRspline2D
    ${JSONCPP_LIBRARIES}
    )

  # Qt uic
  file(GLOB viewlib_volandlr_UIS ${viewlib_SOURCE_DIR}/include/GoTools/viewlib/vol_and_lr/*.ui)
  if(GoTools_USE_Qt4)
    qt4_wrap_ui_go(viewlib_volandlr_UIS_H "include/GoTools/viewlib/vol_and_lr" ${viewlib_volandlr_UIS})
  else()
    # message("Vol & LR: Calling QT5 UI_GO macro!")
    qt5_wrap_ui_go(viewlib_volandlr_UIS_H "include/GoTools/viewlib/vol_and_lr" ${viewlib_volandlr_UIS})
  endif()
  #  qt4_wrap_ui_go_volume(viewlib_volandlr_UIS_H ${viewlib_volandlr_UIS})
  # message(STATUS "Wrapping volume!")
  # message(STATUS "Done wrapping volume!")
  # message(STATUS "viewlib_volandlr_UIS: "${viewlib_volandlr_UIS}) # The source ui-files.
  # message(STATUS "viewlib_volandlr_UIS_H: "${viewlib_volandlr_UIS_H}) # The generated h-files.

  # Qt moc. The option '-nw' is "not recommended" according to the Qt
  # documentation, but it gets rid of warnings...
  # message(STATUS "Starting with moc stuff!")
  file(GLOB viewlib_volandlr_MOC_HDRS_TMP ${viewlib_SOURCE_DIR}/include/GoTools/viewlib/vol_and_lr/*.h)
  # message(STATUS "viewlib_volandlr_MOC_HDRS_TMP:"${viewlib_volandlr_MOC_HDRS_TMP})
  set(viewlib_volandlr_MOC_HDRS ${viewlib_volandlr_MOC_HDRS_TMP} ${viewlib_volandlr_UIS_H})
  if(GoTools_USE_Qt4)
    qt4_wrap_cpp(viewlib_volandlr_MOC_SRCS ${viewlib_volandlr_MOC_HDRS}
      OPTIONS "-nw")
  else()
    # message("Vol & LR: Calling QT5 CPP macro!")
    qt5_wrap_cpp(viewlib_volandlr_MOC_SRCS ${viewlib_volandlr_MOC_HDRS}
      OPTIONS "-nw")
    # qt5_wrap_cpp(viewlib_volandlr_MOC_SRCS ${viewlib_volandlr_MOC_HDRS}
    #   OPTIONS "-nw"))
  endif()
  # message(STATUS ${viewlib_volandlr_MOC_HDRS})
  # message(STATUS ${viewlib_volandlr_MOC_SRCS})
  # message(STATUS "Done with moc stuff!")

  # Make the viewlibVolume library
  file(GLOB_RECURSE viewlib_SRCS src/*.C include/*.h)
  if(BUILD_AS_SHARED_LIBRARY)
    add_library(viewlibVolAndLR SHARED ${viewlib_SRCS} ${viewlib_MOC_SRCS} ${viewlib_volandlr_SRCS} ${viewlib_volandlr_MOC_SRCS})
  else()
    add_library(viewlibVolAndLR ${viewlib_SRCS} ${viewlib_MOC_SRCS} ${viewlib_volandlr_SRCS} ${viewlib_volandlr_MOC_SRCS})
  endif()
  set_property(TARGET viewlibVolAndLR
    PROPERTY FOLDER "viewlib/Libs")

  # # Make the viewlibVolume library.
  # file(GLOB viewlibVolume_src ${viewlib_SOURCE_DIR}/src/volume/*.C)
  # add_library(viewlibVolume STATIC ${viewlibVolume_src})
  # message(STATUS viewlibVolume_src: ${viewlibVolume_src})
  # # target_link_libraries(app_BenchSolvers 
  # #   LinSolverCBlasDep  ${CBLAS_LIBRARY} ${ATLAS_LIBRARY})

  file(GLOB viewlibVolAndLR_APPS app/vol_and_lr/*.C)
  foreach(app ${viewlibVolAndLR_APPS})
    get_filename_component(appname ${app} NAME_WE)
    # message("app: " ${app})
    # message("appname: " ${appname})
    add_executable(${appname} ${app})
    if(NOT GoTools_USE_Qt4)
     qt5_use_modules(${appname} Core Widgets Gui OpenGL)
    endif()
    target_link_libraries(${appname} viewlibVolAndLR ${DEPLIBS})
    set_target_properties(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    set_property(TARGET ${appname}
      PROPERTY FOLDER "viewlib/Apps")
  endforeach(app)

endif(ViewlibVolAndLRSupport)

# Create a tmp dir for generated data.
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/tmp)

#ADD_SUBDIRECTORY(src)

# 'install' target
if(WIN32)
  # Windows
  # lib
  install(TARGETS viewlib DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  install(DIRECTORY include/GoTools/viewlib
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
  # include ui files
  install(DIRECTORY ${viewlib_BINARY_DIR}/include/GoTools/viewlib/
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools/viewlib
    FILES_MATCHING PATTERN "ui_*_form.h"
    PATTERN ".svn" EXCLUDE
    )
  # data
  install(DIRECTORY data/
    DESTINATION ${GoTools_INSTALL_PREFIX}/share/viewlib
    FILES_MATCHING PATTERN "*.g2" PATTERN "*.igs"
    PATTERN ".svn" EXCLUDE
    )
  # bin
  set(goview_INSTALL_PREFIX "$ENV{PROGRAMFILES}/SINTEF/goview"
    CACHE PATH "Path to install goview")
  install(TARGETS goview DESTINATION ${goview_INSTALL_PREFIX})
else()
  # Linux
  # lib
  install(TARGETS viewlib DESTINATION lib COMPONENT viewlib)
  # include
  install(DIRECTORY include/GoTools/viewlib
    DESTINATION include/GoTools
    COMPONENT viewlib-dev
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
  # include ui files
  install(DIRECTORY ${viewlib_BINARY_DIR}/include/GoTools/viewlib/
    DESTINATION include/GoTools/viewlib
    COMPONENT viewlib-dev
    FILES_MATCHING PATTERN "ui_*_form.h"
    PATTERN ".svn" EXCLUDE
    )
  # data
  install(DIRECTORY data/
    DESTINATION share/GoTools/viewlib
    COMPONENT viewlib
    FILES_MATCHING PATTERN "*.g2" PATTERN "*.igs"
    PATTERN ".svn" EXCLUDE
    )
  # bin
  install(TARGETS goview DESTINATION bin COMPONENT viewlib)
endif()

set(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libviewlib.so)
