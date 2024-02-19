MAKEFLAGS = -j$(nproc)
CC = g++
CFLAGS = -std=c++17 -pthread -Wno-comment

.PHONY: no-target
no-target:
	@echo "\033[91mError\033[0m: Please select one of these targets: make Linux-X11, make Linux, make macOS, make Android"

Linux-X11 Linux macOS Android: LDFLAGS_OPENCL = -I./src/OpenCL/include

Linux-X11 Linux: LDLIBS_OPENCL = -L./src/OpenCL/lib -lOpenCL
macOS: LDLIBS_OPENCL = -framework OpenCL
Android: LDLIBS_OPENCL = -L/system/vendor/lib64 -lOpenCL

Linux-X11: LDFLAGS_X11 = -I./src/X11/include
Linux macOS Android: LDFLAGS_X11 =

Linux-X11: LDLIBS_X11 = -L./src/X11/lib -lX11 -lXrandr
Linux macOS Android: LDLIBS_X11 =

Linux-X11 Linux macOS Android: bin/FluidX3D

bin/FluidX3D: temp/graphics.o temp/info.o temp/kernel.o temp/lbm.o temp/lodepng.o temp/main.o temp/setup.o temp/shapes.o make.sh
	@mkdir -p bin
	$(CC) temp/*.o -o bin/FluidX3D $(CFLAGS) $(LDFLAGS_OPENCL) $(LDLIBS_OPENCL) $(LDFLAGS_X11) $(LDLIBS_X11)

temp/graphics.o: src/graphics.cpp src/defines.hpp src/graphics.hpp src/lodepng.hpp src/utilities.hpp make.sh
	@mkdir -p temp
	$(CC) -c src/graphics.cpp -o temp/graphics.o $(CFLAGS) $(LDFLAGS_X11)

temp/info.o: src/info.cpp src/defines.hpp src/graphics.hpp src/info.hpp src/lbm.hpp src/lodepng.hpp src/opencl.hpp src/units.hpp src/utilities.hpp make.sh
	@mkdir -p temp
	$(CC) -c src/info.cpp -o temp/info.o $(CFLAGS) $(LDFLAGS_OPENCL)

temp/kernel.o: src/kernel.cpp src/kernel.hpp src/lodepng.hpp src/utilities.hpp
	@mkdir -p temp
	$(CC) -c src/kernel.cpp -o temp/kernel.o $(CFLAGS)

temp/lbm.o: src/lbm.cpp src/defines.hpp src/graphics.hpp src/info.hpp src/lbm.hpp src/lodepng.hpp src/opencl.hpp src/units.hpp src/utilities.hpp make.sh
	@mkdir -p temp
	$(CC) -c src/lbm.cpp -o temp/lbm.o $(CFLAGS) $(LDFLAGS_OPENCL)

temp/lodepng.o: src/lodepng.cpp src/lodepng.hpp
	@mkdir -p temp
	$(CC) -c src/lodepng.cpp -o temp/lodepng.o $(CFLAGS)

temp/main.o: src/main.cpp src/defines.hpp src/graphics.hpp src/info.hpp src/lbm.hpp src/lodepng.hpp src/opencl.hpp src/setup.hpp src/shapes.hpp src/units.hpp src/utilities.hpp make.sh
	@mkdir -p temp
	$(CC) -c src/main.cpp -o temp/main.o $(CFLAGS) $(LDFLAGS_OPENCL)

temp/setup.o: src/setup.cpp src/defines.hpp src/graphics.hpp src/info.hpp src/lbm.hpp src/lodepng.hpp src/opencl.hpp src/setup.hpp src/shapes.hpp src/units.hpp src/utilities.hpp make.sh
	@mkdir -p temp
	$(CC) -c src/setup.cpp -o temp/setup.o $(CFLAGS) $(LDFLAGS_OPENCL)

temp/shapes.o: src/shapes.cpp src/defines.hpp src/graphics.hpp src/info.hpp src/lbm.hpp src/lodepng.hpp src/opencl.hpp src/shapes.hpp src/units.hpp src/utilities.hpp make.sh
	@mkdir -p temp
	$(CC) -c src/shapes.cpp -o temp/shapes.o $(CFLAGS) $(LDFLAGS_OPENCL)

.PHONY: clean
clean:
	@rm -rf temp bin/FluidX3D
