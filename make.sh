# command line argument $1: device ID; if empty, FluidX3D will automatically choose the fastest available device

mkdir -p bin # create directory for executable
rm -f ./bin/FluidX3D.exe # prevent execution of old version if compiling fails

#g++ ./src/*.cpp -o ./bin/FluidX3D.exe -std=c++17 -pthread -I./src/OpenCL/include -L./src/OpenCL/lib -lOpenCL -I./src/X11/include -L./src/X11/lib -lX11 # compile on Linux with X11

g++ ./src/*.cpp -o ./bin/FluidX3D.exe -std=c++17 -pthread -I./src/OpenCL/include -L./src/OpenCL/lib -lOpenCL # compile on Linux (without X11)
#g++ ./src/*.cpp -o ./bin/FluidX3D -std=c++17 -pthread -I./src/OpenCL/include -framework OpenCL # compile on macOS (without X11)
#g++ ./src/*.cpp -o ./bin/FluidX3D.exe -std=c++17 -pthread -I./src/OpenCL/include -L/system/vendor/lib64 -lOpenCL # compile on Android (without X11)

./bin/FluidX3D.exe $1 # run FluidX3D
