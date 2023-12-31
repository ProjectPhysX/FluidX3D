# command line argument(s): device ID(s); if empty, FluidX3D will automatically choose the fastest available device(s)

#TARGET=Linux-X11 # compile on Linux with X11 graphics
TARGET=Linux # compile on Linux (without X11)
#TARGET=macOS # compile on macOS (without X11)
#TARGET=Android # compile on Android (without X11)

if command -v make &>/dev/null; then # if make is available, compile FluidX3D with multiple CPU cores
	make $TARGET
else # else (make is not installed), compile FluidX3D with a single CPU core
	mkdir -p bin # create directory for executable
	rm -rf temp bin/FluidX3D # prevent execution of old executable if compiling fails
	if [ $TARGET == Linux-X11 ]; then g++ src/*.cpp -o bin/FluidX3D -std=c++17 -pthread -I./src/OpenCL/include -L./src/OpenCL/lib -lOpenCL -I./src/X11/include -L./src/X11/lib -lX11; fi
	if [ $TARGET == Linux     ]; then g++ src/*.cpp -o bin/FluidX3D -std=c++17 -pthread -I./src/OpenCL/include -L./src/OpenCL/lib -lOpenCL; fi
	if [ $TARGET == macOS     ]; then g++ src/*.cpp -o bin/FluidX3D -std=c++17 -pthread -I./src/OpenCL/include -framework OpenCL; fi
	if [ $TARGET == Android   ]; then g++ src/*.cpp -o bin/FluidX3D -std=c++17 -pthread -I./src/OpenCL/include -L/system/vendor/lib64 -lOpenCL; fi
fi

bin/FluidX3D "$@" # run FluidX3D