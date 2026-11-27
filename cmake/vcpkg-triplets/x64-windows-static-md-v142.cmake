# Like x64-windows-static-md, but built with the Visual Studio 2019 toolset
# (v142). Selected by gotools/cmake/VcpkgSetup.cmake when configuring with
# "-T v142".
set(VCPKG_TARGET_ARCHITECTURE x64)
set(VCPKG_CRT_LINKAGE dynamic)
set(VCPKG_LIBRARY_LINKAGE static)
set(VCPKG_PLATFORM_TOOLSET v142)
