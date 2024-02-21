#!/usr/bin/env bash
# command line argument(s) for make.sh: device ID(s); if empty, FluidX3D will automatically choose the fastest available device(s)

case "$(uname -a)" in # automatically detect operating system and X11 support on Linux
	 Darwin*) target=macOS                                                       ;;
	*Android) target=Android                                                     ;;
	 Linux* ) if xhost >&/dev/null; then target=Linux-X11; else target=Linux; fi ;;
	*       ) target=Linux                                                       ;;
esac

#target=Linux-X11 # manually set to compile on Linux with X11 graphics
#target=Linux     # manually set to compile on Linux (without X11)
#target=macOS     # manually set to compile on macOS (without X11)
#target=Android   # manually set to compile on Android (without X11)

echo -e "\033[92mInfo\033[0m: Detected Operating System: "${target}
echo_and_execute() { echo "$@"; "$@"; }
if command -v make &>/dev/null; then # if make is available, compile FluidX3D with multiple CPU cores
	echo -e "\033[92mInfo\033[0m: Compiling with "$(nproc)" CPU cores."
	make ${target} -j$(nproc) # compile FluidX3D with makefile
else # else (make is not installed), compile FluidX3D with a single CPU core
	echo -e "\033[92mInfo\033[0m: Compiling with 1 CPU core. For faster multi-core compiling, install make with \"sudo apt install make\"."
	mkdir -p bin # create directory for executable
	rm -rf temp bin/FluidX3D # prevent execution of old executable if compiling fails
	case "${target}" in
		Linux-X11) echo_and_execute g++ src/*.cpp -o bin/FluidX3D -std=c++17 -pthread -Wno-comment -I./src/OpenCL/include -L./src/OpenCL/lib -lOpenCL -I./src/X11/include -L./src/X11/lib -lX11 -lXrandr ;;
		Linux    ) echo_and_execute g++ src/*.cpp -o bin/FluidX3D -std=c++17 -pthread -Wno-comment -I./src/OpenCL/include -L./src/OpenCL/lib -lOpenCL                                                    ;;
		macOS    ) echo_and_execute g++ src/*.cpp -o bin/FluidX3D -std=c++17 -pthread -Wno-comment -I./src/OpenCL/include -framework OpenCL                                                              ;;
		Android  ) echo_and_execute g++ src/*.cpp -o bin/FluidX3D -std=c++17 -pthread -Wno-comment -I./src/OpenCL/include -L/system/vendor/lib64 -lOpenCL                                                ;;
	esac
fi

if [[ $? == 0 ]]; then bin/FluidX3D "$@"; fi # run FluidX3D only if last compilation was successful
