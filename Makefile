help:
	@echo "Welcome to FluidX3D! \n\nTo get started run 'make' with any of the options: '[linux,mac,android]' depending on your platform. \nIf you'd like to benchmark your OpenCL devices, you may either run 'make [benchmark]' after having run make for your desired plaform or run 'make [platform] benchmark' to do both at once. \n\nIf a device ID is provided by passing 'device=[device ID number]' the benchmark will run on the selected device. Otherwise the most powerful device will be auto selected."


linux:
	@echo "Compile for Linux!"
	@g++ ./src/*.cpp -o ./bin/FluidX3D -std=c++17 -pthread -I./src/OpenCL/include -L./src/OpenCL/lib -lOpenCL -O3 -march=native 


mac:
	@echo "Compile for Mac!"
	@g++ ./src/*.cpp -o ./bin/FluidX3D -std=c++17 -pthread -I./src/OpenCL/include -framework OpenCL # compile on macOS


android:
	@echo "Compile for Android!"
	@g++ ./src/*.cpp -o ./bin/FluidX3D -std=c++17 -pthread -I./src/OpenCL/include -L/system/vendor/lib64 -lOpenCL # compile on Android


clean:
	@echo 'Clean up and remove the old binary; prevent execution of old version if compiling fails'
	@rm ./bin/FluidX3D

benchmark:
	@echo "Running benchmark:"
	@./bin/FluidX3D $(device)
