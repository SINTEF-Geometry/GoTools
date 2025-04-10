project(viewlib)

set(CMAKE_MODULE_PATH
	${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

if(GoTools_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif(GoTools_ENABLE_OPENMP)

find_package(OpenGL REQUIRED)
find_path(OPENGL_INCLUDE_DIR
  GL/gl.h
  # PATH_PREFIXES GL
  PATHS
  "$ENV{HOME}\\Install\\include" # For some reason HOME is not found in VS 2019.
  "C:\\Users\\$ENV{USERNAME}\\Install\\include"
  "C:\\local\\include"
  /usr/local/include
  /usr/include
  )
option(GoTools_USE_Qt4 "Use Qt4 (not Qt5)?" OFF)

# Handle Qt4 and Qt5 specific configurations
if(GoTools_USE_Qt4)
  find_package(Qt4 REQUIRED)
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
  
  # Automatically handle moc, uic, and rcc files for Qt5
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif()

if(WIN32)
  set(GLUT_INCLUDE_DIR
    "$ENV{PROGRAMW6432}/Microsoft SDKs/Windows/v6.0A/Include"
    )
else()
  find_package(GLUT REQUIRED)
endif()

# Qt4 macros for UI and MOC processing
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
  include(${QT_USE_FILE})
else()
  # Qt5 does not require a macro, as CMake_AUTOUIC handles it automatically
endif()

set(QT_USE_QTOPENGL 1)

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
)

if(GoTools_USE_Qt4)
  include_directories(
    ${QT_INCLUDES}
  )
else()
  include_directories(
    ${Qt5Widgets_INCLUDE_DIRS}
    ${Qt5Gui_INCLUDE_DIRS}
  )
endif()

# Linked in libraries

set(DEPLIBS
  GoIgeslib
  GoToolsCore
  sisl
  ${OPENGL_LIBRARIES}
)

if(GoTools_USE_Qt4)
  set(DEPLIBS ${DEPLIBS}
    ${QT_LIBRARIES}
  )
else()
  if(NOT (WIN32 OR APPLE))
    set(DEPLIBS ${DEPLIBS}
      Qt5OpenGL
    )
  endif()
  if(NOT APPLE)
    set(DEPLIBS ${DEPLIBS}
      ${Qt5Widgets_LIBRARIES}
      ${Qt5Core_LIBRARIES}
      ${Qt5Gui_LIBRARIES}
    )
  endif()
endif()

# Qt uic

# We do not use GLOB_RECURSE to avoid including volume directory.
# For the file to be generated we need to pass it through the QT4_WRAP_CPP macro.

# Manually specify the paths for AUTOUIC if needed
set(CMAKE_AUTOUIC_SEARCH_PATHS
    "${viewlib_SOURCE_DIR}/include/GoTools/viewlib"
    "${viewlib_SOURCE_DIR}/include/GoTools/viewlib/vol_and_lr"
)
file(GLOB viewlib_UIS ${viewlib_SOURCE_DIR}/include/GoTools/viewlib/*.ui)

if(GoTools_USE_Qt4)
  qt4_wrap_ui_go(viewlib_UIS_H "include/GoTools/viewlib" ${viewlib_UIS})
endif()

# Qt moc
file(GLOB viewlib_MOC_HDRS_TMP ${viewlib_SOURCE_DIR}/include/GoTools/viewlib/*.h)
set(viewlib_MOC_HDRS ${viewlib_MOC_HDRS_TMP} ${viewlib_UIS_H})

if(GoTools_USE_Qt4)
  qt4_wrap_cpp(viewlib_MOC_SRCS ${viewlib_MOC_HDRS}
    OPTIONS "-nw")
endif()

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
  
  if(NOT GoTools_USE_Qt4)
    qt5_use_modules(${appname} Core Widgets Gui OpenGL)
  endif()

  target_link_libraries(${appname} viewlib ${DEPLIBS})
  set_target_properties(${appname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    set_property(TARGET ${appname}
      PROPERTY FOLDER "viewlib/Apps")
endforeach(app)

# Copy Qt dll's.
IF(WIN32)
  if(GoTools_USE_Qt4)
    set(QT_DLLS QtWidgets4.dll QtGui4.dll QtCore4.dll QtOpenGL4.dll opengl32sw.dll)
    set(QT_DEBUG_DLLS QtWidgets4d.dll QtGui4d.dll QtCore4d.dll QtOpenGL4d.dll opengl32sw.dll)
  else()
    set(QT_DLLS Qt5Widgets.dll Qt5Gui.dll Qt5Core.dll Qt5OpenGL.dll opengl32sw.dll)
    set(QT_DEBUG_DLLS Qt5Widgetsd.dll Qt5Guid.dll Qt5Cored.dll Qt5OpenGLd.dll opengl32sw.dll)
  endif()
  # First the release dll's.
  file(MAKE_DIRECTORY ${viewlib_BINARY_DIR}\\app\\Release)
  foreach( file_i ${QT_DLLS})
    ADD_CUSTOM_COMMAND(
      TARGET viewlib
      POST_BUILD
      COMMAND ${CMAKE_COMMAND}
      ARGS -E copy $ENV{QTDIR}\\bin\\${file_i} ${viewlib_BINARY_DIR}\\app\\Release\\
    )
  endforeach(file_i)
  # Then the debug dll's.
  file(MAKE_DIRECTORY ${viewlib_BINARY_DIR}\\app\\Debug)
  foreach( file_i ${QT_DEBUG_DLLS})
    ADD_CUSTOM_COMMAND(
      TARGET viewlib
      POST_BUILD
      COMMAND ${CMAKE_COMMAND}
      ARGS -E copy $ENV{QTDIR}\\bin\\${file_i} ${viewlib_BINARY_DIR}\\app\\Debug\\
    )
  endforeach(file_i)
ENDIF(WIN32)

OPTION(ViewlibVolAndLRSupport "Enable support for volumes and lr models?" ON)

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
  endif()

  file(GLOB viewlib_volandlr_MOC_HDRS_TMP ${viewlib_SOURCE_DIR}/include/GoTools/viewlib/vol_and_lr/*.h)
  set(viewlib_volandlr_MOC_HDRS ${viewlib_volandlr_MOC_HDRS_TMP} ${viewlib_volandlr_UIS_H})

  # Qt moc. The option '-nw' is "not recommended" according to the Qt
  # documentation, but it gets rid of warnings...
  if(GoTools_USE_Qt4)
    qt4_wrap_cpp(viewlib_volandlr_MOC_SRCS ${viewlib_volandlr_MOC_HDRS}
      OPTIONS "-nw")
  endif()

  # Make the viewlibVolume library.
  file(GLOB_RECURSE viewlib_SRCS src/*.C include/*.h)
  if(BUILD_AS_SHARED_LIBRARY)
    add_library(viewlibVolAndLR SHARED ${viewlib_SRCS} ${viewlib_MOC_SRCS} ${viewlib_volandlr_SRCS} ${viewlib_volandlr_MOC_SRCS})
  else()
    add_library(viewlibVolAndLR ${viewlib_SRCS} ${viewlib_MOC_SRCS} ${viewlib_volandlr_SRCS} ${viewlib_volandlr_MOC_SRCS})
  endif()
  set_property(TARGET viewlibVolAndLR
    PROPERTY FOLDER "viewlib/Libs")

  file(GLOB viewlibVolAndLR_APPS app/vol_and_lr/*.C)
  foreach(app ${viewlibVolAndLR_APPS})
    get_filename_component(appname ${app} NAME_WE)
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

# === Install Targets and Export Package ===

include(CMakePackageConfigHelpers)

# Setup standard install directories
if(WIN32)
  set(CMAKE_INSTALL_LIBDIR ${GoTools_INSTALL_PREFIX}/lib)
  set(CMAKE_INSTALL_INCLUDEDIR ${GoTools_INSTALL_PREFIX}/include)
  set(CMAKE_INSTALL_BINDIR ${GoTools_INSTALL_PREFIX}/bin)
else()
  include(GNUInstallDirs)
endif()

# Install viewlib target
install(TARGETS viewlib
    EXPORT viewlibTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install viewlib binaries
# Path to your app build folder
set(APP_BUILD_DIR "${viewlib_BINARY_DIR}/app")

# Determine configuration subdir (e.g., Debug/Release) — used by multi-config generators
if(CMAKE_CONFIGURATION_TYPES)
  # Multi-config generator (e.g., Visual Studio, Xcode)
  set(CONFIG_SUBDIR "$<CONFIG>")
else()
  # Single-config generator (e.g., Ninja, Unix Makefiles)
  set(CONFIG_SUBDIR "")
endif()

# Path to use at install time — app/[config]/ if multi-config, app/ otherwise
set(APP_RUNTIME_DIR "${APP_BUILD_DIR}/${CONFIG_SUBDIR}")
#message("APP_RUNTIME_DIR: ${APP_RUNTIME_DIR}")
# Install everything from app/[config]/ (or app/) into install bin directory
install(DIRECTORY "${APP_RUNTIME_DIR}/"
  DESTINATION ${CMAKE_INSTALL_BINDIR}
  FILES_MATCHING PATTERN "*"
  PATTERN ".svn" EXCLUDE
)

#message(STATUS "Will install app runtime files from: ${APP_RUNTIME_DIR}")

# Install headers (including generated ui_*.h headers)
install(DIRECTORY include/GoTools/viewlib
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
)

install(DIRECTORY ${viewlib_BINARY_DIR}/include/GoTools/viewlib/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/GoTools/viewlib
    FILES_MATCHING PATTERN "ui_*_form.h"
    PATTERN ".svn" EXCLUDE
)

# Export the targets
install(EXPORT viewlibTargets
    FILE viewlibTargets.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/viewlib
)

# Generate version file
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/viewlibConfigVersion.cmake"
    VERSION ${GoTools_VERSION_MAJOR}.${GoTools_VERSION_MINOR}.${GoTools_VERSION_PATCH}
    COMPATIBILITY SameMajorVersion
)

# Generate config file from template
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/viewlibConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/viewlibConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/viewlib
)

# Install config files
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/viewlibConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/viewlibConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/viewlib
)

# Optional stripping
set(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libviewlib.so)
