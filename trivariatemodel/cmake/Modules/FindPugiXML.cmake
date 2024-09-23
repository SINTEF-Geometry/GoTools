# Find the pugixml XML parsing library.
#
# Sets the usual variables expected for find_package scripts:
#
# PUGIXML_INCLUDE_DIR - header location
# PUGIXML_LIBRARIES - library to link against
# PUGIXML_FOUND - true if pugixml was found.

#Find pugixml header
find_path(PUGIXML_INCLUDE_DIR
  NAMES pugixml.hpp
  PATHS
  ${PUGIXML_HOME}/include/pugi
  /usr/local/include/pugixml-1.8/
  "/usr/include"
  "$ENV{ProgramFiles}/Microsoft Visual Studio 8/VC/PlatformSDK/Include"
  "$ENV{ProgramFiles}/Microsoft Visual Studio 9.0/VC/include/"
  "$ENV{PROGRAMW6432}/Microsoft SDKs/Windows/v6.0A/Include"
  "~/mylibs/pugixml/include"
  "~/Install/include/pugi"
  "C:\\local\\include\\pugi"
  "$ENV{ProgramFiles}/Microsoft SDKs/Windows/v7.0A/Include"
  "${PUGIXML_ROOT}/include"
)

if(WIN32)
  if(${MSVC_VERSION} EQUAL 1900)
    set(MSVC_NAME "msvc2015_")
    # MESSAGE("Visual Studio 2015!")
  elseif((${MSVC_VERSION} GREATER_EQUAL 1920) AND (${MSVC_VERSION} LESS 1930))
    # MESSAGE("Visual Studio 2019!")
    set(MSVC_NAME "msvc2019_")
  elseif((${MSVC_VERSION} GREATER_EQUAL 1930) AND (${MSVC_VERSION} LESS 1950))
    set(MSVC_NAME "msvc2022_")
  else()
    message("MSVC version not supported or not installed!")
  endif()
  if(CMAKE_CL_64)
    set(WIN_LIB_TYPE "64")
  else()
    set(WIN_LIB_TYPE "32")
  endif()
endif()

# Find pugixml lib
find_library(PUGIXML_LIBRARY 
  NAMES pugixml
  PATHS
  "~/Install/${MSVC_NAME}lib${WIN_LIB_TYPE}"
  "C:/local/${MSVC_NAME}lib${WIN_LIB_TYPE}"
   ${PUGIXML_HOME}/lib
  /usr/local/lib/pugixml-1.8
  "/usr/lib"
  "/usr/lib64"
  "$ENV{ProgramFiles}/Microsoft Visual Studio 8/VC/PlatformSDK/Lib"
  "$ENV{ProgramFiles}/Microsoft Visual Studio 9.0/VC/lib/"
  "$ENV{PROGRAMW6432}/Microsoft SDKs/Windows/v6.0A/Lib"
  "~/mylibs/pugixml/lib"
  "$ENV{ProgramFiles}/Microsoft SDKs/Windows/v7.0A/Lib"
  "${PUGIXML_ROOT}/lib"
  )

# Find pugixml debug
find_library(PUGIXML_LIBRARY_DEBUG
  NAMES pugixml
  PATHS
  "~/Install/${MSVC_NAME}lib${WIN_LIB_TYPE}/Debug"
  "C:/local/${MSVC_NAME}lib${WIN_LIB_TYPE}/Debug"
  )

# Find pugixml release lib
find_library(PUGIXML_LIBRARY_RELEASE
  NAMES pugixml
  PATHS
  "~/Install/${MSVC_NAME}lib${WIN_LIB_TYPE}/Release"
  "C:/local/${MSVC_NAME}lib${WIN_LIB_TYPE}/Release"
  )
# message("PUGIXML_LIBRARY_RELEASE:" ${PUGIXML_LIBRARY_RELEASE})
if(PUGIXML_LIBRARY_DEBUG)
  set(PUGIXML_LIBRARIES ${PUGIXML_LIBRARIES} debug ${PUGIXML_LIBRARY_DEBUG})
endif()
if(PUGIXML_LIBRARY_RELEASE)
  set(PUGIXML_LIBRARIES ${PUGIXML_LIBRARIES} optimized ${PUGIXML_LIBRARY_RELEASE})
endif()
if(PUGIXML_LIBRARY)
  set(PUGIXML_LIBRARIES ${PUGIXML_LIBRARIES} ${PUGIXML_LIBRARY})
endif()

# Check that we have found everything
set(PUGIXML_FOUND "NO")
if(PUGIXML_INCLUDE_DIR AND PUGIXML_LIBRARIES)
    set(PUGIXML_FOUND "YES")
endif()

# Mark options as advanced (i.e. they will not be shown in the cmake GUI)
mark_as_advanced(
  PUGIXML_INCLUDE_DIR
  PUGIXML_LIBRARY
  PUGIXML_LIBRARIES
)