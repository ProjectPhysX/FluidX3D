#pragma once

#include "defines.hpp"
#include "opencl.hpp"
#include "graphics.hpp"

class LBM {
private:
	Device device; // OpenCL device associated with this LBM object
	Kernel kernel_initialize; // initialization kernel
	Kernel kernel_stream_collide; // main LBM kernel
	Kernel kernel_update_fields; // reads DDFs and updates (rho, u, T) in device memory
	Memory<fpxx> fi; // LBM density distribution functions (DDFs); only exist in device memory
	ulong t_last_update_fields = 0ull; // optimization to not call kernel_update_fields multiple times if (rho, u, T) are already up-to-date

#if defined(D2Q9)
	const uint velocity_set = 9u;
	const uint dimensions = 2u;
#elif defined(D3Q15)
	const uint velocity_set = 15u;
	const uint dimensions = 3u;
#elif defined(D3Q19)
	const uint velocity_set = 19u;
	const uint dimensions = 3u;
#elif defined(D3Q27)
	const uint velocity_set = 27u;
	const uint dimensions = 3u;
#endif // D3Q27

#ifdef FORCE_FIELD
	Kernel kernel_calculate_force_on_boundaries; // calculate forces from fluid on TYPE_S nodes
#endif // FORCE_FIELD

#ifdef MOVING_BOUNDARIES
	Kernel kernel_update_moving_boundaries; // mark/unmark nodes next to TYPE_S nodes with velocity!=0 with TYPE_MS
#endif // MOVING_BOUNDARIES

#ifdef SURFACE
	Kernel kernel_surface_0; // additional kernel for computing mass conservation and mass flux computation
	Kernel kernel_surface_1; // additional kernel for flag handling
	Kernel kernel_surface_2; // additional kernel for flag handling
	Kernel kernel_surface_3; // additional kernel for flag handling and mass conservation
	Memory<float> mass; // fluid mass; phi=mass/rho
	Memory<float> massex; // excess mass; used for mass conservation
#endif // SURFACE

#ifdef TEMPERATURE
	Memory<fpxx> gi; // thermal DDFs
#endif // TEMPERATURE

	bool initialized = false;
	void sanity_checks_constructor(); // sanity checks during constructor call on grid resolution and parameters
	void sanity_checks_initialization(); // sanity checks during initialization on used extensions based on used flags
	void allocate(Device& device); // allocate all memory for data fields on host and device and set up kernels
	void initialize(); // write all data fields to device and call kernel_initialize
	void do_time_step(); // call kernel_stream_collide to perform one LBM time step
	string device_defines() const; // returns preprocessor constants for embedding in OpenCL C code

public:
	uint Nx=1u, Ny=1u, Nz=1u; // lattice dimensions
	float nu = 1.0f/6.0f; // kinematic shear viscosity
	float fx=0.0f, fy=0.0f, fz=0.0f; // global force per volume
	float sigma=0.0f; // surface tension coefficient
	float T_avg=1.0f, alpha=1.0f, beta=1.0f; // T_avg = 1 = average temperature, alpha = thermal diffusion coefficient, beta = thermal expansion coefficient
	ulong t = 0ull; // discrete time step in LBM units

	Memory<float> rho; // density of every node
	Memory<float> u; // velocity of every node
	Memory<uchar> flags; // flags of every node

#ifdef FORCE_FIELD
	Memory<float> F; // individual force for every node
#endif // FORCE_FIELD

#ifdef SURFACE
	Memory<float> phi; // fill level of every node
#endif // SURFACE

#ifdef TEMPERATURE
	Memory<float> T; // temperature of every node
#endif // TEMPERATURE

	LBM(const uint Nx, const uint Ny, const uint Nz, const float nu, const float fx=0.0f, const float fy=0.0f, const float fz=0.0f, const float sigma=0.0f, const float alpha=0.0f, const float beta=0.0f);
	~LBM();

	void set_fx(const float fx) { this->fx = fx; } // set global froce per volume
	void set_fy(const float fy) { this->fy = fy; } // set global froce per volume
	void set_fz(const float fz) { this->fz = fz; } // set global froce per volume
	void set_f(const float fx, const float fy, const float fz) { this->fx = fx; this->fy = fy; this->fz = fz; } // set global froce per volume

	uint get_Nx() const { return Nx; }
	uint get_Ny() const { return Ny; }
	uint get_Nz() const { return Nz; }
	uint get_N() const { return Nx*Ny*Nz; }
	float get_nu() const { return nu; }
	float get_tau() const { return 3.0f*nu+0.5f; }
	float get_fx() const { return fx; }
	float get_fy() const { return fy; }
	float get_fz() const { return fz; }
	float get_sigma() const { return sigma; }
	float get_T_avg() const { return T_avg; }
	float get_alpha() const { return alpha; }
	float get_beta() const { return beta; }
	ulong get_t() const { return t; }
	uint get_velocity_set() const { return velocity_set; }
	float get_Re_max() const { return 0.57735027f*(float)min(min(Nx, Ny), Nz)/nu; } // Re < c*L/nu
	void coordinates(const uint n, uint& x, uint& y, uint& z) const { // disassemble 1D linear index to 3D coordinates (n -> x,y,z)
		const uint t = n%(Nx*Ny); // n = x+(y+z*Ny)*Nx
		x = t%Nx;
		y = t/Nx;
		z = n/(Nx*Ny);
	}
	uint index(const uint x, const uint y, const uint z) const { // turn 3D coordinates into 1D linear index
		return x+(y+z*Ny)*Nx;
	}
	float3 position(const uint x, const uint y, const uint z) const { // returns position in box [-Nx/2, Nx/2] x [-Ny/2, Ny/2] x [-Nz/2, Nz/2]
		return float3((float)x-0.5f*(float)Nx+0.5f, (float)y-0.5f*(float)Ny+0.5f, (float)z-0.5f*(float)Nz+0.5f);
	}
	float3 position(const uint n) const { // returns position in box [-Nx/2, Nx/2] x [-Ny/2, Ny/2] x [-Nz/2, Nz/2]
		uint x, y, z;
		coordinates(n, x, y, z);
		return position(x, y, z);
	}
	float3 size() const { // returns size of box
		return float3((float)Nx, (float)Ny, (float)Nz);
	}
	float3 center() const { // returns center of box
		return float3(0.5f*(float)Nx-0.5f, 0.5f*(float)Ny-0.5f, 0.5f*(float)Nz-0.5f);
	}
	uint smallest_side_length() const {
		return min(min(Nx, Ny), Nz);
	}
	uint largest_side_length() const {
		return max(max(Nx, Ny), Nz);
	}
	float3 relative_position(const uint x, const uint y, const uint z) const { // returns relative position in box [-0.5, 0.5] x [-0.5, 0.5] x [-0.5, 0.5]
		return float3(((float)x+0.5f)/(float)Nx-0.5f, ((float)y+0.5f)/(float)Ny-0.5f, ((float)z+0.5f)/(float)Nz-0.5f);
	}
	float3 relative_position(const uint n) const { // returns relative position in box [-0.5, 0.5] x [-0.5, 0.5] x [-0.5, 0.5]
		uint x, y, z;
		coordinates(n, x, y, z);
		return relative_position(x, y, z);
	}

	void run(const ulong steps=max_ulong); // initializes the LBM simulation (copies data to device and runs initialize kernel), then runs LBM
	void update_fields(); // update fields (rho, u, T) manually
	void reset(); // reset simulation (takes effect in following run() call)

	string default_filename(const string& path, const string& name, const string& extension); // generate a default filename with timestamp
	string default_filename(const string& name, const string& extension); // generate a default filename with timestamp at exe_path/export/

#ifdef FORCE_FIELD
	void calculate_force_on_boundaries(); // calculate forces from fluid on TYPE_S nodes
#endif // FORCE_FIELD

#ifdef MOVING_BOUNDARIES
	void update_moving_boundaries(); // mark/unmark nodes next to TYPE_S nodes with velocity!=0 with TYPE_MS
#endif // MOVING_BOUNDARIES

	void write_status(const string& path=""); // write LBM status report to a .txt file

	void rho_write_host_to_vtk(const string& path=""); // write binary .vtk file
	void rho_write_device_to_vtk(const string& path=""); // write binary .vtk file
	void u_write_host_to_vtk(const string& path=""); // write binary .vtk file
	void u_write_device_to_vtk(const string& path=""); // write binary .vtk file
	void flags_write_host_to_vtk(const string& path=""); // write binary .vtk file
	void flags_write_device_to_vtk(const string& path=""); // write binary .vtk file
#ifdef FORCE_FIELD
	void F_write_host_to_vtk(const string& path=""); // write binary .vtk file
	void F_write_device_to_vtk(const string& path=""); // write binary .vtk file
#endif // FORCE_FIELD
#ifdef SURFACE
	void phi_write_host_to_vtk(const string& path=""); // write binary .vtk file
	void phi_write_device_to_vtk(const string& path=""); // write binary .vtk file
#endif // SURFACE
#ifdef TEMPERATURE
	void T_write_host_to_vtk(const string& path=""); // write binary .vtk file
	void T_write_device_to_vtk(const string& path=""); // write binary .vtk file
#endif // TEMPERATURE

	void voxelize_mesh(const Mesh* mesh, const uchar flag=TYPE_S); // voxelize mesh
	void voxelize_stl(const string& path, const float3& center, const float3x3& rotation, const float size=-1.0f, const uchar flag=TYPE_S); // read and voxelize binary .stl file
	void voxelize_stl(const string& path, const float3x3& rotation, const float size=-1.0f, const uchar flag=TYPE_S); // read and voxelize binary .stl file (place in box center)
	void voxelize_stl(const string& path, const float3& center, const float size=-1.0f, const uchar flag=TYPE_S); // read and voxelize binary .stl file (no rotation)
	void voxelize_stl(const string& path, const float size=-1.0f, const uchar flag=TYPE_S); // read and voxelize binary .stl file (place in box center, no rotation)

#ifdef GRAPHICS
	class Graphics {
	private:
		Kernel kernel_clear; // reset bitmap and zbuffer
		Memory<uint> bitmap; // bitmap for rendering
		Memory<int> zbuffer; // z-buffer for rendering
		Memory<float> camera_parameters; // contains camera position, rotation, field of view etc.

		LBM* lbm = nullptr;
		ulong t_last_frame = 0ull; // optimization to not call draw_frame() multiple times if camera_parameters and LBM time step are unchanged
		std::atomic_int running_encoders = 0;

		Kernel kernel_graphics_flags; // render flag lattice
		Kernel kernel_graphics_field; // render a colored velocity vector for each node
		Kernel kernel_graphics_streamline; // render streamlines
		Kernel kernel_graphics_q; // render vorticity (Q-criterion)

#ifdef SURFACE
		const string path_skybox = get_exe_path()+"../skybox/skybox8k.png";
		Image* skybox_image = nullptr;
		Memory<uint> skybox; // skybox for free surface raytracing
		Kernel kernel_graphics_rasterize_phi; // rasterize free surface
		Kernel kernel_graphics_raytrace_phi; // raytrace free surface
#endif // SURFACE

		void default_settings(); // set what is viewed at startup
		bool update_camera(); // update camera_parameters and return if they are changed from their previous state

	public:
		Graphics() {} // default constructor
		~Graphics() { // destructor must wait for all encoder threads to finish
			int last_value = running_encoders.load();
			while(last_value>0) {
				const int current_value = running_encoders.load();
				if(last_value!=current_value) {
					print_info("Finishing encoder threads: "+to_string(current_value));
					last_value = current_value;
				}
				sleep(0.016);
			}
		}
		Graphics(LBM* lbm) {
			this->lbm = lbm;
#ifdef SURFACE
			skybox_image = read_png(path_skybox);
#endif // SURFACE
		}
		Graphics& operator=(const Graphics& graphics) { // copy assignment
			lbm = graphics.lbm;
#ifdef SURFACE
			skybox_image = graphics.get_skybox_image();
#endif // SURFACE
			return *this;
		}
#ifdef SURFACE
		Image* get_skybox_image() const { return skybox_image; }
#endif // SURFACE

		void allocate(Device& device); // allocate memory for bitmap and zbuffer
		void* draw_frame(); // main rendering function, calls rendering kernels
		string device_defines() const; // returns preprocessor constants for embedding in OpenCL C code

		void set_camera_centered(const float rx=0.0f, const float ry=0.0f, const float fov=100.0f, const float zoom=1.0f); // set camera centered
		void set_camera_free(const float3& p=float3(0.0f), const float rx=0.0f, const float ry=0.0f, const float fov=100.0f); // set camera free
		void print_frame(); // preview preview of current frame in console
		void write_frame(const string& path="", const string& name="image", const string& extension=".png", bool print_preview=false); // save current frame
		void write_frame(const uint x1, const uint y1, const uint x2, const uint y2, const string& path="", const string& name="image", const string& extension=".png", bool print_preview=false); // save current frame cropped with two corner points (x1,y1) and (x2,y2)
		void write_frame_png(const string& path="", bool print_preview=false); // save current frame as .png file (smallest file size, but slow)
		void write_frame_qoi(const string& path="", bool print_preview=false); // save current frame as .qoi file (small file size, fast)
		void write_frame_bmp(const string& path="", bool print_preview=false); // save current frame as .bmp file (large file size, fast)
		void write_frame_png(const uint x1, const uint y1, const uint x2, const uint y2, const string& path="", bool print_preview=false); // save current frame as .png file (smallest file size, but slow)
		void write_frame_qoi(const uint x1, const uint y1, const uint x2, const uint y2, const string& path="", bool print_preview=false); // save current frame as .qoi file (small file size, fast)
		void write_frame_bmp(const uint x1, const uint y1, const uint x2, const uint y2, const string& path="", bool print_preview=false); // save current frame as .bmp file (large file size, fast)

	}; // Graphics
	Graphics graphics;
#endif // GRAPHICS
}; // LBM