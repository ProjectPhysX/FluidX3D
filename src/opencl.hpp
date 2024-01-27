#pragma once

#define WORKGROUP_SIZE 64 // needs to be 64 to fully use AMD GPUs
//#define PTX
//#define LOG

#ifndef _WIN32
#pragma GCC diagnostic ignored "-Wignored-attributes" // ignore compiler warnings for CL/cl.hpp with g++
#endif // _WIN32
#include <CL/cl.hpp> // OpenCL 1.0, 1.1, 1.2
#include "utilities.hpp"
using cl::Event;

struct Device_Info {
	cl::Device cl_device; // OpenCL device
	cl::Context cl_context; // multiple devices in the same context can communicate buffers
	uint id = 0u; // unique device ID assigned by get_devices()
	string name, vendor; // device name, vendor
	string driver_version, opencl_c_version; // device driver version, OpenCL C version
	uint memory=0u; // global memory in MB
	uint memory_used=0u; // track global memory usage in MB
	uint global_cache=0u, local_cache=0u; // global cache in KB, local cache in KB
	uint max_global_buffer=0u, max_constant_buffer=0u; // maximum global buffer size in MB, maximum constant buffer size in KB
	uint compute_units=0u; // compute units (CUs) can contain multiple cores depending on the microarchitecture
	uint clock_frequency=0u; // in MHz
	bool is_cpu=false, is_gpu=false;
	bool intel_gpu_above_4gb_patch = false; // memory allocations greater than 4GB need to be specifically enabled on Intel GPUs
	uint is_fp64_capable=0u, is_fp32_capable=0u, is_fp16_capable=0u, is_int64_capable=0u, is_int32_capable=0u, is_int16_capable=0u, is_int8_capable=0u;
	uint cores=0u; // for CPUs, compute_units is the number of threads (twice the number of cores with hyperthreading)
	float tflops=0.0f; // estimated device FP32 floating point performance in TeraFLOPs/s
	inline Device_Info(const cl::Device& cl_device, const cl::Context& cl_context, const uint id) {
		this->cl_device = cl_device; // see https://www.khronos.org/registry/OpenCL/sdk/1.2/docs/man/xhtml/clGetDeviceInfo.html
		this->cl_context = cl_context;
		this->id = id;
		name = trim(cl_device.getInfo<CL_DEVICE_NAME>()); // device name
		vendor = trim(cl_device.getInfo<CL_DEVICE_VENDOR>()); // device vendor
		driver_version = trim(cl_device.getInfo<CL_DRIVER_VERSION>()); // device driver version
		opencl_c_version = trim(cl_device.getInfo<CL_DEVICE_OPENCL_C_VERSION>()); // device OpenCL C version
		memory = (uint)(cl_device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>()/1048576ull); // global memory in MB
		global_cache = (uint)(cl_device.getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_SIZE>()/1024ull); // global cache in KB
		local_cache = (uint)(cl_device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>()/1024ull); // local cache in KB
		max_global_buffer = (uint)(cl_device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>()/1048576ull); // maximum global buffer size in MB
		max_constant_buffer = (uint)(cl_device.getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>()/1024ull); // maximum constant buffer size in KB
		compute_units = (uint)cl_device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>(); // compute units (CUs) can contain multiple cores depending on the microarchitecture
		clock_frequency = (uint)cl_device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>(); // in MHz
		is_fp64_capable = (uint)cl_device.getInfo<CL_DEVICE_NATIVE_VECTOR_WIDTH_DOUBLE>()*(uint)contains(cl_device.getInfo<CL_DEVICE_EXTENSIONS>(), "cl_khr_fp64");
		is_fp32_capable = (uint)cl_device.getInfo<CL_DEVICE_NATIVE_VECTOR_WIDTH_FLOAT>();
		is_fp16_capable = (uint)cl_device.getInfo<CL_DEVICE_NATIVE_VECTOR_WIDTH_HALF>()*(uint)contains(cl_device.getInfo<CL_DEVICE_EXTENSIONS>(), "cl_khr_fp16");
		is_int64_capable = (uint)cl_device.getInfo<CL_DEVICE_NATIVE_VECTOR_WIDTH_LONG>();
		is_int32_capable = (uint)cl_device.getInfo<CL_DEVICE_NATIVE_VECTOR_WIDTH_INT>();
		is_int16_capable = (uint)cl_device.getInfo<CL_DEVICE_NATIVE_VECTOR_WIDTH_SHORT>();
		is_int8_capable = (uint)cl_device.getInfo<CL_DEVICE_NATIVE_VECTOR_WIDTH_CHAR>();
		is_cpu = cl_device.getInfo<CL_DEVICE_TYPE>()==CL_DEVICE_TYPE_CPU;
		is_gpu = cl_device.getInfo<CL_DEVICE_TYPE>()==CL_DEVICE_TYPE_GPU;
		const uint ipc = is_gpu?2u:32u; // IPC (instructions per cycle) is 2 for GPUs and 32 for most modern CPUs
		const bool nvidia_192_cores_per_cu = contains_any(to_lower(name), {"gt 6", "gt 7", "gtx 6", "gtx 7", "quadro k", "tesla k"}) || (clock_frequency<1000u&&contains(to_lower(name), "titan")); // identify Kepler GPUs
		const bool nvidia_64_cores_per_cu = contains_any(to_lower(name), {"p100", "v100", "a100", "a30", " 16", " 20", "titan v", "titan rtx", "quadro t", "tesla t", "quadro rtx"}) && !contains(to_lower(name), "rtx a"); // identify P100, Volta, Turing, A100, A30
		const bool amd_128_cores_per_dualcu = contains(to_lower(name), "gfx10"); // identify RDNA/RDNA2 GPUs where dual CUs are reported
		const bool amd_256_cores_per_dualcu = contains(to_lower(name), "gfx11"); // identify RDNA3 GPUs where dual CUs are reported
		const bool intel_16_cores_per_cu = contains(to_lower(name), "gpu max"); // identify PVC GPUs
		const float nvidia = (float)(contains(to_lower(vendor), "nvidia"))*(nvidia_64_cores_per_cu?64.0f:nvidia_192_cores_per_cu?192.0f:128.0f); // Nvidia GPUs have 192 cores/CU (Kepler), 128 cores/CU (Maxwell, Pascal, Ampere, Hopper, Ada) or 64 cores/CU (P100, Volta, Turing, A100, A30)
		const float amd = (float)(contains_any(to_lower(vendor), {"amd", "advanced"}))*(is_gpu?(amd_256_cores_per_dualcu?256.0f:amd_128_cores_per_dualcu?128.0f:64.0f):0.5f); // AMD GPUs have 64 cores/CU (GCN, CDNA), 128 cores/dualCU (RDNA, RDNA2) or 256 cores/dualCU (RDNA3), AMD CPUs (with SMT) have 1/2 core/CU
		const float intel = (float)(contains(to_lower(vendor), "intel"))*(is_gpu?(intel_16_cores_per_cu?16.0f:8.0f):0.5f); // Intel GPUs have 16 cores/CU (PVC) or 8 cores/CU (integrated/Arc), Intel CPUs (with HT) have 1/2 core/CU
		const float apple = (float)(contains(to_lower(vendor), "apple"))*(128.0f); // Apple ARM GPUs usually have 128 cores/CU
		const float arm = (float)(contains(to_lower(vendor), "arm"))*(is_gpu?8.0f:1.0f); // ARM GPUs usually have 8 cores/CU, ARM CPUs have 1 core/CU
		cores = to_uint((float)compute_units*(nvidia+amd+intel+apple+arm)); // for CPUs, compute_units is the number of threads (twice the number of cores with hyperthreading)
		tflops = 1E-6f*(float)cores*(float)ipc*(float)clock_frequency; // estimated device floating point performance in TeraFLOPs/s
		if(intel==8.0f) { // fix wrong global memory reporting for Intel Arc GPUs
			if((contains_any(name, {"A770", "0x56a0"})&&memory>=11739u&&memory<14168u)||(contains_any(name, {"A770", "A750", "A580", "0x56a0", "0x56a1", "0x56a2"})&&memory>=5869u&&memory<7084u)||(contains_any(name, {"A380", "0x56a5"})&&memory>=4402u&&memory<5313u)) { // 72.5%-87.5% reporting -> /0.8
				memory = (uint)((cl_device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>()*5ull/4ull)/1048576ull); // fix wrong (80% on Windows) memory reporting on Intel Arc
			}
			if((contains_any(name, {"A770", "0x56a0"})&&memory>=14168u&&memory<15625u)||(contains_any(name, {"A770", "A750", "A580", "0x56a0", "0x56a1", "0x56a2"})&&memory>=7084u&&memory<7812u)||(contains_any(name, {"A380", "0x56a5"})&&memory>=5313u&&memory<5859u)) { // 87.5%-96.5% reporting -> /0.95
				memory = (uint)((cl_device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>()*20ull/19ull)/1048576ull); // fix wrong (95% on Linux) memory reporting on Intel Arc
			}
			if((contains_any(name, {"A770", "0x56a0"})&&memory>=15625u&&memory<16030u)||(contains_any(name, {"A770", "A750", "A580", "0x56a0", "0x56a1", "0x56a2"})&&memory>=7812u&&memory<8015u)||(contains_any(name, {"A380", "0x56a5"})&&memory>=5859u&&memory<6011u)) { // 96.5%-99.0% reporting -> /0.98
				memory = (uint)((cl_device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>()*50ull/49ull)/1048576ull); // fix wrong (98% on Windows) memory reporting on Intel Arc
			}
		}
		intel_gpu_above_4gb_patch = (intel==8.0f)&&(memory>4096); // enable memory allocations greater than 4GB for Intel GPUs with >4GB VRAM
	}
	inline Device_Info() {}; // default constructor
};

string get_opencl_c_code(); // implemented in kernel.hpp
inline void print_device_info(const Device_Info& d) { // print OpenCL device info
#if defined(_WIN32)
	const string os = "Windows";
#elif defined(__linux__)
	const string os = "Linux";
#elif defined(__APPLE__)
	const string os = "macOS";
#else // unknown operating system
	const string os = "unknown operating system";
#endif // operating system
	println("\r|----------------.------------------------------------------------------------|");
	println("| Device ID      | "+alignl(58, to_string(d.id)             )+" |");
	println("| Device Name    | "+alignl(58, d.name                      )+" |");
	println("| Device Vendor  | "+alignl(58, d.vendor                    )+" |");
	println("| Device Driver  | "+alignl(58, d.driver_version+" ("+os+")")+" |");
	println("| OpenCL Version | "+alignl(58, d.opencl_c_version          )+" |");
	println("| Compute Units  | "+alignl(58, to_string(d.compute_units)+" at "+to_string(d.clock_frequency)+" MHz ("+to_string(d.cores)+" cores, "+to_string(d.tflops, 3)+" TFLOPs/s)")+" |");
	println("| Memory, Cache  | "+alignl(58, to_string(d.memory)+" MB, "+to_string(d.global_cache)+" KB global / "+to_string(d.local_cache)+" KB local")+" |");
	println("| Buffer Limits  | "+alignl(58, to_string(d.max_global_buffer)+" MB global, "+to_string(d.max_constant_buffer)+" KB constant")+" |");
	println("|----------------'------------------------------------------------------------|");
}
inline vector<Device_Info> get_devices(const bool print_info=true) { // returns a vector of all available OpenCL devices
	vector<Device_Info> devices; // get all devices of all platforms
	vector<cl::Platform> cl_platforms; // get all platforms (drivers)
	cl::Platform::get(&cl_platforms);
	uint id = 0u;
	for(uint i=0u; i<(uint)cl_platforms.size(); i++) {
		vector<cl::Device> cl_devices;
		cl_platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &cl_devices);
		//cl::Context cl_context(cl_devices); // same cl::Context for all devices (allocates extra VRAM on all other unused Nvidia GPUs)
		for(uint j=0u; j<(uint)cl_devices.size(); j++) {
			cl::Context cl_context(cl_devices[j]); // separate cl::Context for each device
			devices.push_back(Device_Info(cl_devices[j], cl_context, id++));
		}
	}
	if((uint)cl_platforms.size()==0u||(uint)devices.size()==0u) {
		print_error("There are no OpenCL devices available. Make sure that the OpenCL 1.2 Runtime for your device is installed. For GPUs it comes by default with the graphics driver, for CPUs it has to be installed separately.");
	}
	if(print_info) {
		println("\r|----------------.------------------------------------------------------------|");
		for(uint i=0u; i<(uint)devices.size(); i++) println("| Device ID "+alignr(4u, i)+" | "+alignl(58u, devices[i].name)+" |");
		println("|----------------'------------------------------------------------------------|");
	}
	return devices;
}
inline Device_Info select_device_with_most_flops(const vector<Device_Info>& devices=get_devices()) { // returns device with best floating-point performance
	float best_value = 0.0f;
	uint best_i = 0u;
	for(uint i=0u; i<(uint)devices.size(); i++) { // find device with highest (estimated) floating point performance
		if(devices[i].tflops>best_value) {
			best_value = devices[i].tflops;
			best_i = i;
		}
	}
	return devices[best_i];
}
inline Device_Info select_device_with_most_memory(const vector<Device_Info>& devices=get_devices()) { // returns device with largest memory capacity
	uint best_value = 0u;
	uint best_i = 0u;
	for(uint i=0u; i<(uint)devices.size(); i++) { // find device with most memory
		if(devices[i].memory>best_value) {
			best_value = devices[i].memory;
			best_i = i;
		}
	}
	return devices[best_i];
}
inline Device_Info select_device_with_id(const uint id, const vector<Device_Info>& devices=get_devices()) { // returns device with specified ID
	if(id<(uint)devices.size()) {
		return devices[id];
	} else {
		print_error("Your selected Device ID ("+to_string(id)+") is wrong.");
		return devices[0]; // is never executed, just to avoid compiler warnings
	}
}

class Device {
private:
	cl::Program cl_program;
	cl::CommandQueue cl_queue;
	bool exists = false;
	inline string enable_device_capabilities() const { return // enable FP64/FP16 capabilities if available
		"\n	#define def_workgroup_size "+to_string(WORKGROUP_SIZE)+"u"
		"\n	#ifdef cl_khr_fp64"
		"\n	#pragma OPENCL EXTENSION cl_khr_fp64 : enable" // make sure cl_khr_fp64 extension is enabled
		"\n	#endif"
		"\n	#ifdef cl_khr_fp16"
		"\n	#pragma OPENCL EXTENSION cl_khr_fp16 : enable" // make sure cl_khr_fp16 extension is enabled
		"\n	#endif"
		"\n	#ifdef cl_khr_int64_base_atomics"
		"\n	#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable" // make sure cl_khr_int64_base_atomics extension is enabled
		"\n	#endif"
	;}
public:
	Device_Info info;
	inline Device(const Device_Info& info, const string& opencl_c_code=get_opencl_c_code()) {
		print_device_info(info);
		this->info = info;
		this->cl_queue = cl::CommandQueue(info.cl_context, info.cl_device); // queue to push commands for the device
		cl::Program::Sources cl_source;
		const string kernel_code = enable_device_capabilities()+"\n"+opencl_c_code;
		cl_source.push_back({ kernel_code.c_str(), kernel_code.length() });
		this->cl_program = cl::Program(info.cl_context, cl_source);
		const string build_options = string("-cl-finite-math-only -cl-no-signed-zeros -cl-mad-enable")+(info.intel_gpu_above_4gb_patch ? " -cl-intel-greater-than-4GB-buffer-required" : "");
#ifndef LOG
		int error = cl_program.build({ info.cl_device }, (build_options+" -w").c_str()); // compile OpenCL C code, disable warnings
		if(error) print_warning(cl_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(info.cl_device)); // print build log
#else // LOG, generate logfile for OpenCL code compilation
		int error = cl_program.build({ info.cl_device }, build_options.c_str()); // compile OpenCL C code
		const string log = cl_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(info.cl_device);
		write_file("bin/kernel.log", log); // save build log
		if((uint)log.length()>2u) print_warning(log); // print build log
#endif // LOG
		if(error) print_error("OpenCL C code compilation failed with error code "+to_string(error)+". Make sure there are no errors in kernel.cpp.");
		else print_info("OpenCL C code successfully compiled.");
#ifdef PTX // generate assembly (ptx) file for OpenCL code
		write_file("bin/kernel.ptx", cl_program.getInfo<CL_PROGRAM_BINARIES>()[0]); // save binary (ptx file)
#endif // PTX
		this->exists = true;
	}
	inline Device() {} // default constructor
	inline void barrier(const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) { cl_queue.enqueueBarrierWithWaitList(event_waitlist, event_returned); }
	inline void finish_queue() { cl_queue.finish(); }
	inline cl::Context get_cl_context() const { return info.cl_context; }
	inline cl::Program get_cl_program() const { return cl_program; }
	inline cl::CommandQueue get_cl_queue() const { return cl_queue; }
	inline bool is_initialized() const { return exists; }
};

template<typename T> class Memory {
private:
	ulong N = 0ull; // buffer length
	uint d = 1u; // buffer dimensions
	bool host_buffer_exists = false;
	bool device_buffer_exists = false;
	bool external_host_buffer = false;
	T* host_buffer = nullptr; // host buffer
	cl::Buffer device_buffer; // device buffer
	Device* device = nullptr; // pointer to linked Device
	cl::CommandQueue cl_queue; // command queue
	inline void initialize_auxiliary_pointers() {
		/********/ x = s0 = host_buffer; /******/ if(d>0x4u) s4 = host_buffer+N*0x4ull; if(d>0x8u) s8 = host_buffer+N*0x8ull; if(d>0xCu) sC = host_buffer+N*0xCull;
		if(d>0x1u) y = s1 = host_buffer+N; /****/ if(d>0x5u) s5 = host_buffer+N*0x5ull; if(d>0x9u) s9 = host_buffer+N*0x9ull; if(d>0xDu) sD = host_buffer+N*0xDull;
		if(d>0x2u) z = s2 = host_buffer+N*0x2ull; if(d>0x6u) s6 = host_buffer+N*0x6ull; if(d>0xAu) sA = host_buffer+N*0xAull; if(d>0xEu) sE = host_buffer+N*0xEull;
		if(d>0x3u) w = s3 = host_buffer+N*0x3ull; if(d>0x7u) s7 = host_buffer+N*0x7ull; if(d>0xBu) sB = host_buffer+N*0xBull; if(d>0xFu) sF = host_buffer+N*0xFull;
	}
	inline void allocate_device_buffer(Device& device, const bool allocate_device) {
		this->device = &device;
		this->cl_queue = device.get_cl_queue();
		if(allocate_device) {
			device.info.memory_used += (uint)(capacity()/1048576ull); // track device memory usage
			if(device.info.memory_used>device.info.memory) print_error("Device \""+device.info.name+"\" does not have enough memory. Allocating another "+to_string((uint)(capacity()/1048576ull))+" MB would use a total of "+to_string(device.info.memory_used)+" MB / "+to_string(device.info.memory)+" MB.");
			int error = 0;
			device_buffer = cl::Buffer(device.get_cl_context(), CL_MEM_READ_WRITE|((int)device.info.intel_gpu_above_4gb_patch<<23), capacity(), nullptr, &error); // for Intel GPUs, set flag CL_MEM_ALLOW_UNRESTRICTED_SIZE_INTEL = (1<<23)
			if(error==-61) print_error("Memory size is too large at "+to_string((uint)(capacity()/1048576ull))+" MB. Device \""+device.info.name+"\" accepts a maximum buffer size of "+to_string(device.info.max_global_buffer)+" MB.");
			else if(error) print_error("Device buffer allocation failed with error code "+to_string(error)+".");
			device_buffer_exists = true;
		}
	}
public:
	T *x=nullptr, *y=nullptr, *z=nullptr, *w=nullptr; // host buffer auxiliary pointers for multi-dimensional array access (array of structures)
	T *s0=nullptr, *s1=nullptr, *s2=nullptr, *s3=nullptr, *s4=nullptr, *s5=nullptr, *s6=nullptr, *s7=nullptr, *s8=nullptr, *s9=nullptr, *sA=nullptr, *sB=nullptr, *sC=nullptr, *sD=nullptr, *sE=nullptr, *sF=nullptr;
	inline Memory(Device& device, const ulong N, const uint dimensions=1u, const bool allocate_host=true, const bool allocate_device=true, const T value=(T)0) {
		if(!device.is_initialized()) print_error("No Device selected. Call Device constructor.");
		if(N*(ulong)dimensions==0ull) print_error("Memory size must be larger than 0.");
		this->N = N;
		this->d = dimensions;
		allocate_device_buffer(device, allocate_device);
		if(allocate_host) {
			host_buffer = new T[N*(ulong)d];
			initialize_auxiliary_pointers();
			host_buffer_exists = true;
		}
		reset(value);
	}
	inline Memory(Device& device, const ulong N, const uint dimensions, T* const host_buffer, const bool allocate_device=true) {
		if(!device.is_initialized()) print_error("No Device selected. Call Device constructor.");
		if(N*(ulong)dimensions==0ull) print_error("Memory size must be larger than 0.");
		this->N = N;
		this->d = dimensions;
		allocate_device_buffer(device, allocate_device);
		this->host_buffer = host_buffer;
		initialize_auxiliary_pointers();
		host_buffer_exists = true;
		external_host_buffer = true;
		write_to_device();
	}
	inline Memory() {} // default constructor
	inline ~Memory() {
		delete_buffers();
	}
	inline Memory& operator=(Memory&& memory) noexcept { // move assignment
		delete_buffers(); // delete existing buffers and restore default state
		N = memory.length(); // copy values/pointers from memory
		d = memory.dimensions();
		device = memory.device;
		cl_queue = memory.device->get_cl_queue();
		if(memory.device_buffer_exists) {
			device_buffer = memory.get_cl_buffer(); // transfer device_buffer pointer
			device->info.memory_used += (uint)(capacity()/1048576ull); // track device memory usage
			device_buffer_exists = true;
		}
		if(memory.host_buffer_exists) {
			host_buffer = memory.exchange_host_buffer(nullptr); // transfer host_buffer pointer
			initialize_auxiliary_pointers();
			host_buffer_exists = true;
		}
		return *this; // destructor of memory will be called automatically
	}
	inline T* const exchange_host_buffer(T* const host_buffer) { // sets host_buffer to new pointer and returns old pointer
		T* const swap = this->host_buffer;
		this->host_buffer = host_buffer;
		return swap;
	}
	inline void add_host_buffer() { // makes only sense if there is no host buffer yet but an existing device buffer
		if(!host_buffer_exists&&device_buffer_exists) {
			host_buffer = new T[N*(ulong)d];
			initialize_auxiliary_pointers();
			read_from_device();
			host_buffer_exists = true;
		} else if(!device_buffer_exists) {
			print_error("There is no existing device buffer, so can't add host buffer.");
		}
	}
	inline void add_device_buffer() { // makes only sense if there is no device buffer yet but an existing host buffer
		if(!device_buffer_exists&&host_buffer_exists) {
			allocate_device_buffer(*device, true);
			write_to_device();
		} else if(!host_buffer_exists) {
			print_error("There is no existing host buffer, so can't add device buffer.");
		}
	}
	inline void delete_host_buffer() {
		host_buffer_exists = false;
		if(!external_host_buffer) delete[] host_buffer;
		if(!device_buffer_exists) {
			N = 0ull;
			d = 1u;
		}
	}
	inline void delete_device_buffer() {
		if(device_buffer_exists) device->info.memory_used -= (uint)(capacity()/1048576ull); // track device memory usage
		device_buffer_exists = false;
		device_buffer = nullptr;
		if(!host_buffer_exists) {
			N = 0ull;
			d = 1u;
		}
	}
	inline void delete_buffers() {
		delete_device_buffer();
		delete_host_buffer();
	}
	inline void reset(const T value=(T)0) {
		if(host_buffer_exists) {
			std::fill(host_buffer, host_buffer+range(), value); // faster than "for(ulong i=0ull; i<range(); i++) host_buffer[i] = value;"
		}
		write_to_device(); // enqueueFillBuffer is broken for large buffers on Nvidia GPUs!
		//if(device_buffer_exists) {
		//	cl_queue.enqueueFillBuffer(device_buffer, value, 0ull, capacity()); // faster than "write_to_device();"
		//	cl_queue.finish();
		//}
	}
	inline const ulong length() const { return N; }
	inline const uint dimensions() const { return d; }
	inline const ulong range() const { return N*(ulong)d; }
	inline const ulong capacity() const { return N*(ulong)d*sizeof(T); } // returns capacity of the buffer in Byte
	inline T* const data() { return host_buffer; }
	inline const T* const data() const { return host_buffer; }
	inline T* const operator()() { return host_buffer; }
	inline const T* const operator()() const { return host_buffer; }
	inline T& operator[](const ulong i) { return host_buffer[i]; }
	inline const T& operator[](const ulong i) const { return host_buffer[i]; }
	inline const T operator()(const ulong i) const { return host_buffer[i]; }
	inline const T operator()(const ulong i, const uint dimension) const { return host_buffer[i+(ulong)dimension*N]; } // array of structures
	inline void read_from_device(const bool blocking=true, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) {
		if(host_buffer_exists&&device_buffer_exists) cl_queue.enqueueReadBuffer(device_buffer, blocking, 0ull, capacity(), (void*)host_buffer, event_waitlist, event_returned);
	}
	inline void write_to_device(const bool blocking=true, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) {
		if(host_buffer_exists&&device_buffer_exists) cl_queue.enqueueWriteBuffer(device_buffer, blocking, 0ull, capacity(), (void*)host_buffer, event_waitlist, event_returned);
	}
	inline void read_from_device(const ulong offset, const ulong length, const bool blocking=true, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) {
		if(host_buffer_exists&&device_buffer_exists) {
			const ulong safe_offset=min(offset, range()), safe_length=min(length, range()-safe_offset);
			if(safe_length>0ull) cl_queue.enqueueReadBuffer(device_buffer, blocking, safe_offset*sizeof(T), safe_length*sizeof(T), (void*)(host_buffer+safe_offset), event_waitlist, event_returned);
		}
	}
	inline void write_to_device(const ulong offset, const ulong length, const bool blocking=true, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) {
		if(host_buffer_exists&&device_buffer_exists) {
			const ulong safe_offset=min(offset, range()), safe_length=min(length, range()-safe_offset);
			if(safe_length>0ull) cl_queue.enqueueWriteBuffer(device_buffer, blocking, safe_offset*sizeof(T), safe_length*sizeof(T), (void*)(host_buffer+safe_offset), event_waitlist, event_returned);
		}
	}
	inline void read_from_device_1d(const ulong x0, const ulong x1, const int dimension=-1, const bool blocking=true, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) { // read 1D domain from device, either for all vector dimensions (-1) or for a specified dimension
		if(host_buffer_exists&&device_buffer_exists) {
			const uint i0=(uint)max(0, dimension), i1=dimension<0 ? d : i0+1u;
			for(uint i=i0; i<i1; i++) {
				const ulong safe_offset=min((ulong)i*N+x0, range()), safe_length=min(x1-x0, range()-safe_offset);
				if(safe_length>0ull) cl_queue.enqueueReadBuffer(device_buffer, false, safe_offset*sizeof(T), safe_length*sizeof(T), (void*)(host_buffer+safe_offset), event_waitlist, event_returned);
			}
			if(blocking) cl_queue.finish();
		}
	}
	inline void write_to_device_1d(const ulong x0, const ulong x1, const int dimension=-1, const bool blocking=true, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) { // write 1D domain to device, either for all vector dimensions (-1) or for a specified dimension
		if(host_buffer_exists&&device_buffer_exists) {
			const uint i0=(uint)max(0, dimension), i1=dimension<0 ? d : i0+1u;
			for(uint i=i0; i<i1; i++) {
				const ulong safe_offset=min((ulong)i*N+x0, range()), safe_length=min(x1-x0, range()-safe_offset);
				if(safe_length>0ull) cl_queue.enqueueWriteBuffer(device_buffer, false, safe_offset*sizeof(T), safe_length*sizeof(T), (void*)(host_buffer+safe_offset), event_waitlist, event_returned);
			}
			if(blocking) cl_queue.finish();
		}
	}
	inline void read_from_device_2d(const ulong x0, const ulong x1, const ulong y0, const ulong y1, const ulong Nx, const ulong Ny, const int dimension=-1, const bool blocking=true, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) { // read 2D domain from device, either for all vector dimensions (-1) or for a specified dimension
		if(host_buffer_exists&&device_buffer_exists) {
			for(uint y=y0; y<y1; y++) {
				const ulong n = x0+y*Nx;
				const uint i0=(uint)max(0, dimension), i1=dimension<0 ? d : i0+1u;
				for(uint i=i0; i<i1; i++) {
					const ulong safe_offset=min((ulong)i*N+n, range()), safe_length=min(x1-x0, range()-safe_offset);
					if(safe_length>0ull) cl_queue.enqueueReadBuffer(device_buffer, false, safe_offset*sizeof(T), safe_length*sizeof(T), (void*)(host_buffer+safe_offset), event_waitlist, event_returned);
				}
			}
			if(blocking) cl_queue.finish();
		}
	}
	inline void write_to_device_2d(const ulong x0, const ulong x1, const ulong y0, const ulong y1, const ulong Nx, const ulong Ny, const int dimension=-1, const bool blocking=true, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) { // write 2D domain to device, either for all vector dimensions (-1) or for a specified dimension
		if(host_buffer_exists&&device_buffer_exists) {
			for(uint y=y0; y<y1; y++) {
				const ulong n = x0+y*Nx;
				const uint i0=(uint)max(0, dimension), i1=dimension<0 ? d : i0+1u;
				for(uint i=i0; i<i1; i++) {
					const ulong safe_offset=min((ulong)i*N+n, range()), safe_length=min(x1-x0, range()-safe_offset);
					if(safe_length>0ull) cl_queue.enqueueWriteBuffer(device_buffer, false, safe_offset*sizeof(T), safe_length*sizeof(T), (void*)(host_buffer+safe_offset), event_waitlist, event_returned);
				}
			}
			if(blocking) cl_queue.finish();
		}
	}
	inline void read_from_device_3d(const ulong x0, const ulong x1, const ulong y0, const ulong y1, const ulong z0, const ulong z1, const ulong Nx, const ulong Ny, const ulong Nz, const int dimension=-1, const bool blocking=true, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) { // read 3D domain from device, either for all vector dimensions (-1) or for a specified dimension
		if(host_buffer_exists&&device_buffer_exists) {
			for(uint z=z0; z<z1; z++) {
				for(uint y=y0; y<y1; y++) {
					const ulong n = x0+(y+z*Ny)*Nx;
					const uint i0=(uint)max(0, dimension), i1=dimension<0 ? d : i0+1u;
					for(uint i=i0; i<i1; i++) {
						const ulong safe_offset=min((ulong)i*N+n, range()), safe_length=min(x1-x0, range()-safe_offset);
						if(safe_length>0ull) cl_queue.enqueueReadBuffer(device_buffer, false, safe_offset*sizeof(T), safe_length*sizeof(T), (void*)(host_buffer+safe_offset), event_waitlist, event_returned);
					}
				}
			}
			if(blocking) cl_queue.finish();
		}
	}
	inline void write_to_device_3d(const ulong x0, const ulong x1, const ulong y0, const ulong y1, const ulong z0, const ulong z1, const ulong Nx, const ulong Ny, const ulong Nz, const int dimension=-1, const bool blocking=true, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) { // write 3D domain to device, either for all vector dimensions (-1) or for a specified dimension
		if(host_buffer_exists&&device_buffer_exists) {
			for(uint z=z0; z<z1; z++) {
				for(uint y=y0; y<y1; y++) {
					const ulong n = x0+(y+z*Ny)*Nx;
					const uint i0=(uint)max(0, dimension), i1=dimension<0 ? d : i0+1u;
					for(uint i=i0; i<i1; i++) {
						const ulong safe_offset=min((ulong)i*N+n, range()), safe_length=min(x1-x0, range()-safe_offset);
						if(safe_length>0ull) cl_queue.enqueueWriteBuffer(device_buffer, false, safe_offset*sizeof(T), safe_length*sizeof(T), (void*)(host_buffer+safe_offset), event_waitlist, event_returned);
					}
				}
			}
			if(blocking) cl_queue.finish();
		}
	}
	inline void enqueue_read_from_device(const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) { read_from_device(false, event_waitlist, event_returned); }
	inline void enqueue_write_to_device(const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) { write_to_device(false, event_waitlist, event_returned); }
	inline void enqueue_read_from_device(const ulong offset, const ulong length, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) { read_from_device(offset, length, false, event_waitlist, event_returned); }
	inline void enqueue_write_to_device(const ulong offset, const ulong length, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) { write_to_device(offset, length, false, event_waitlist, event_returned); }
	inline void finish_queue() { cl_queue.finish(); }
	inline const cl::Buffer& get_cl_buffer() const { return device_buffer; }
};

class Kernel {
private:
	ulong N = 0ull; // kernel range
	uint number_of_parameters = 0u;
	cl::Kernel cl_kernel;
	cl::NDRange cl_range_global, cl_range_local;
	cl::CommandQueue cl_queue;
	template<typename T> inline void link_parameter(const uint position, const Memory<T>& memory) {
		cl_kernel.setArg(position, memory.get_cl_buffer());
	}
	template<typename T> inline void link_parameter(const uint position, const T& constant) {
		cl_kernel.setArg(position, sizeof(T), (void*)&constant);
	}
	inline void link_parameters(const uint starting_position) {
		number_of_parameters = max(number_of_parameters, starting_position);
	}
	template<class T, class... U> inline void link_parameters(const uint starting_position, const T& parameter, const U&... parameters) {
		link_parameter(starting_position, parameter);
		link_parameters(starting_position+1u, parameters...);
	}
public:
	template<class... T> inline Kernel(const Device& device, const ulong N, const string& name, const T&... parameters) { // accepts Memory<T> objects and fundamental data type constants
		if(!device.is_initialized()) print_error("No Device selected. Call Device constructor.");
		cl_kernel = cl::Kernel(device.get_cl_program(), name.c_str());
		link_parameters(number_of_parameters, parameters...); // expand variadic template to link kernel parameters
		set_ranges(N);
		cl_queue = device.get_cl_queue();
	}
	template<class... T> inline Kernel(const Device& device, const ulong N, const uint workgroup_size, const string& name, const T&... parameters) { // accepts Memory<T> objects and fundamental data type constants
		if(!device.is_initialized()) print_error("No Device selected. Call Device constructor.");
		cl_kernel = cl::Kernel(device.get_cl_program(), name.c_str());
		link_parameters(number_of_parameters, parameters...); // expand variadic template to link kernel parameters
		set_ranges(N, (ulong)workgroup_size);
		cl_queue = device.get_cl_queue();
	}
	inline Kernel() {} // default constructor
	inline Kernel& set_ranges(const ulong N, const ulong workgroup_size=(ulong)WORKGROUP_SIZE) {
		this->N = N;
		cl_range_global = cl::NDRange(((N+workgroup_size-1ull)/workgroup_size)*workgroup_size); // make global range a multiple of local range
		cl_range_local = cl::NDRange(workgroup_size);
		return *this;
	}
	inline const ulong range() const { return N; }
	inline uint get_number_of_parameters() const { return number_of_parameters; }
	template<class... T> inline Kernel& add_parameters(const T&... parameters) { // add parameters to the list of existing parameters
		link_parameters(number_of_parameters, parameters...); // expand variadic template to link kernel parameters
		return *this;
	}
	template<class... T> inline Kernel& set_parameters(const uint starting_position, const T&... parameters) { // set parameters starting at specified position
		link_parameters(starting_position, parameters...); // expand variadic template to link kernel parameters
		return *this;
	}
	inline Kernel& enqueue_run(const uint t=1u, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) {
		for(uint i=0u; i<t; i++) {
			cl_queue.enqueueNDRangeKernel(cl_kernel, cl::NullRange, cl_range_global, cl_range_local, event_waitlist, event_returned);
		}
		return *this;
	}
	inline Kernel& run(const uint t=1u, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) {
		enqueue_run(t, event_waitlist, event_returned);
		finish_queue();
		return *this;
	}
	inline Kernel& operator()(const uint t=1u, const vector<Event>* event_waitlist=nullptr, Event* event_returned=nullptr) {
		return run(t, event_waitlist, event_returned);
	}
	inline Kernel& finish_queue() {
		cl_queue.finish();
		return *this;
	}
};