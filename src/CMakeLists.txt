
set(fluidx3d_HEADERS
  ${fluidx3d_HEADERS}
  defines.hpp
  ghraphics.hpp
  info.hpp
  kernel.hpp
  lbm.hpp
  lodepng.hpp
  opencl.hpp
  setup.hpp
  shapes.hpp
  units.hpp
  utilities.hpp
  )


set(fluidx3d_SOURCES
  ${fluidx3d_SOURCES}
  graphics.cpp
  info.cpp
  kernel.cpp
  lbm.cpp
  lodepng.cpp
  main.cpp
  setup.cpp
  shapes.cpp
)


set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable(${PROJECT_NAME} ${fluidx3d_SOURCES})


target_include_directories( ${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/OpenCL/include>
  )

target_link_libraries(${PROJECT_NAME} PUBLIC
  Threads::Threads
  ${CMAKE_DL_LIBS}
  -L${PROJECT_SOURCE_DIR}/src/OpenCL/lib -lOpenCL
  )
