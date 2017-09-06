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
  "$ENV{ProgramFiles}/Microsoft SDKs/Windows/v7.0A/Include"
  "${PUGIXML_ROOT}/include"
)

if(WIN32)
  if(CMAKE_CL_64)
    set(WIN_LIB_DIR "win64")
#    message("The project is set to 64 bits!")
  else()
    set(WIN_LIB_DIR "win32")
#    message("The project is set to 32 bits!")
  endif()
endif()

# Find pugixml lib
find_library(PUGI_LIBRARY 
  NAMES pugixml
  PATHS
  "~/Install/lib/${WIN_LIB_DIR}"
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
find_library(PUGI_LIBRARY_DEBUG
  NAMES pugixml
  PATHS
  "~/Install/lib/${WIN_LIB_DIR}/Debug"
  )

# Find pugixml release lib
find_library(PUGI_LIBRARY_RELEASE
  NAMES pugixml
  PATHS
  "~/Install/lib/${WIN_LIB_DIR}/Release"
  )

#check that we have found everything
set(PUGIXML_FOUND "NO")
if(PUGIXML_LIBRARY)
  if(PUGIXML_INCLUDE_DIR)
    set(PUGIXML_FOUND "YES")
	set(PUGIXML_LIBRARIES ${PUGIXML_LIBRARY})
  endif()
endif()

#Mark options as advanced
mark_as_advanced(
  PUGIXML_INCLUDE_DIR
  PUGIXML_LIBRARY
  PUGIXML_LIBRARIES
)