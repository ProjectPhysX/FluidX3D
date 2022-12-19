#include "lbm.hpp"
#include "info.hpp"
#include "graphics.hpp"
#include "units.hpp"

Units units; // for unit conversion

LBM::LBM(const uint Nx, const uint Ny, const uint Nz, const float nu, const float fx, const float fy, const float fz, const float sigma, const float alpha, const float beta) { // compiles OpenCL C code and allocates memory
	this->Nx = Nx; this->Ny = Ny; this->Nz = Nz;
	this->nu = nu;
	this->fx = fx; this->fy = fy; this->fz = fz;
	this->sigma = sigma;
	this->alpha = alpha;
	this->beta = beta;
	sanity_checks_constructor();
	string opencl_c_code;
#ifdef GRAPHICS
	graphics = Graphics(this);
	opencl_c_code = device_defines()+graphics.device_defines()+get_opencl_c_code();
#else // GRAPHICS
	opencl_c_code = device_defines()+get_opencl_c_code();
#endif // GRAPHICS
	int select_device = -1;
	if((int)main_arguments.size()>0) select_device = to_int(main_arguments[0]);
	const vector<Device_Info>& devices = get_devices();
	this->device = Device(select_device<0 ? select_device_with_most_flops(devices) : select_device_with_id((uint)select_device, devices), opencl_c_code);
	allocate(device); // lbm first
#ifdef GRAPHICS
	graphics.allocate(device); // graphics after lbm
#endif // GRAPHICS
	info.initialize(this);
}
LBM::~LBM() {
	info.print_finalize();
}

void LBM::sanity_checks_constructor() { // sanity checks on grid resolution and parameters
	if(Nx*Ny*Nz==0u || (Nx*Ny*Nz)%WORKGROUP_SIZE!=0u) { // sanity checks for simulation box size
		int ws=WORKGROUP_SIZE;
		string error = "";
		if(Nx*Ny*Nz!=0u) {
			error += "Lattice point number is not a multiple of \"WORKGROUP_SIZE\": ("+to_string(Nx)+"x"+to_string(Ny)+"x"+to_string(Nz)+") % "+to_string(ws)+" = "+to_string((Nx*Ny*Nz)%ws)+". ";
		} else {
			error += "Lattice point number is 0: "+to_string(Nx)+"x"+to_string(Ny)+"x"+to_string(Nz)+" = 0. ";
		}
		bool suggestion = false; // find the closest suitable box size based on minimum percentage correction for each component
		Nx = max(Nx, 1u); // box dimensions need to be at least 1 for the algorithm to work
		Ny = max(Ny, 1u);
		Nz = max(Nz, 1u);
		const int ix=(int)Nx, iy=(int)Ny, iz=(int)Nz, imax=max(max(ix,iy),iz); // integer variables for later arithmetic
		while(ws>=32&&!suggestion) { // if search is unsuccessful, try again with threadblock size halved, until it is 32
			for(int k=1; k<imax&&!suggestion; k++) { // slowly extend search range k
				for(int kx=0; abs(kx)<=k*ix/imax&&!suggestion; kx=kx<=0?1-kx:-kx) { // 0  1 -1  2 -2  3 -3 ...
					for(int ky=0; abs(ky)<=k*iy/imax&&!suggestion; ky=ky<=0?1-ky:-ky) { // 0  1 -1  2 -2  3 -3 ...
						for(int kz=0; abs(kz)<=k*iz/imax&&!suggestion; kz=kz<=0?1-kz:-kz) { // 0  1 -1  2 -2  3 -3 ...
							if(((ix+kx)*(iy+ky)*(iz+kz))%ws==0 && ix+kx>0 && iy+ky>0 && iz+kz>0) { // criteria for suitable box size
								Nx = (uint)(ix+kx); // set new box size
								Ny = (uint)(iy+ky);
								Nz = (uint)(iz+kz);
								suggestion = true; // possible box size found -> exit all loops
							}
						}
					}
				}
			}
			if(!suggestion) ws /= 2;
		}
		if(suggestion && ws==WORKGROUP_SIZE) {
			error += "Change grid resolution in setup.cpp to ("+to_string(Nx)+"x"+to_string(Ny)+"x"+to_string(Nz)+").";
		} else {
			if(suggestion&&ws!=WORKGROUP_SIZE) {
				error += "Change \"WORKGROUP_SIZE\" to "+to_string(ws)+", then ("+to_string(Nx)+"x"+to_string(Ny)+"x"+to_string(Nz)+") would be a suitable grid resolution.";
			} else {
				error += "Grid resolution in setup.cpp is too small.";
			}
		}
		print_error(error);
	}
	if(nu==0.0f) print_error("Viscosity cannot be 0. Change it in setup.cpp."); // sanity checks for viscosity
	else if(nu<0.0f) print_error("Viscosity cannot be negative. Remove the \"-\" in setup.cpp.");
#ifndef VOLUME_FORCE
	if(fx!=0.0f||fy!=0.0f||fz!=0.0f) print_error("Volume force is set in LBM constructor in main_setup(), but VOLUME_FORCE is not enabled. Uncomment \"#define VOLUME_FORCE\" in defines.hpp.");
#else // VOLUME_FORCE
#ifndef FORCE_FIELD
	if(fx==0.0f&&fy==0.0f&&fz==0.0f) print_warning("The VOLUME_FORCE extension is enabled but the volume force in LBM constructor is set to zero. You may disable the extension by commenting out \"#define VOLUME_FORCE\" in defines.hpp.");
#endif // FORCE_FIELD
#endif // VOLUME_FORCE
#ifndef SURFACE
	if(sigma!=0.0f) print_error("Surface tension is set in LBM constructor in main_setup(), but SURFACE is not enabled. Uncomment \"#define SURFACE\" in defines.hpp.");
#endif // SURFACE
#ifndef TEMPERATURE
	if(alpha!=0.0f||beta!=0.0f) print_error("Thermal diffusion/expansion coefficients are set in LBM constructor in main_setup(), but TEMPERATURE is not enabled. Uncomment \"#define TEMPERATURE\" in defines.hpp.");
#else // TEMPERATURE
	if(alpha==0.0f&&beta==0.0f) print_warning("The TEMPERATURE extension is enabled but the thermal diffusion/expansion coefficients alpha/beta in the LBM constructor are both set to zero. You may disable the extension by commenting out \"#define TEMPERATURE\" in defines.hpp.");
#endif // TEMPERATURE
}
void LBM::sanity_checks_initialization() { // sanity checks during initialization on used extensions based on used flags
	bool moving_boundaries_used=false, equilibrium_boundaries_used=false, surface_used=false, temperature_used=false; // identify used extensions based used flags
	for(uint n=0u; n<get_N(); n++) {
		const uchar flagsn = flags[n];
		const uchar flagsn_bo = flagsn&(TYPE_S|TYPE_E);
		const uchar flagsn_su = flagsn&(TYPE_F|TYPE_I|TYPE_G);
		moving_boundaries_used = moving_boundaries_used || (((flagsn_bo==TYPE_S)&&(u.x[n]!=0.0f||u.y[n]!=0.0f||u.z[n]!=0.0f))||(flagsn_bo==(TYPE_S|TYPE_E)));
		equilibrium_boundaries_used = equilibrium_boundaries_used || (flagsn_bo==TYPE_E);
		surface_used = surface_used || flagsn_su;
		temperature_used = temperature_used || (flagsn&TYPE_T);
	}
#ifndef MOVING_BOUNDARIES
	if(moving_boundaries_used) print_warning("Some boundary nodes have non-zero velocity, but MOVING_BOUNDARIES is not enabled. If you intend to use moving boundaries, uncomment \"#define MOVING_BOUNDARIES\" in defines.hpp.");
#else // MOVING_BOUNDARIES
	if(!moving_boundaries_used) print_warning("The MOVING_BOUNDARIES extension is enabled but no moving boundary nodes (TYPE_S flag and velocity unequal to zero) are placed in the simulation box. You may disable the extension by commenting out \"#define MOVING_BOUNDARIES\" in defines.hpp.");
#endif // MOVING_BOUNDARIES
#ifndef EQUILIBRIUM_BOUNDARIES
	if(equilibrium_boundaries_used) print_error("Some nodes are set as equilibrium boundaries with the TYPE_E flag, but EQUILIBRIUM_BOUNDARIES is not enabled. Uncomment \"#define EQUILIBRIUM_BOUNDARIES\" in defines.hpp.");
#else // EQUILIBRIUM_BOUNDARIES
	if(!equilibrium_boundaries_used) print_warning("The EQUILIBRIUM_BOUNDARIES extension is enabled but no equilibrium boundary nodes (TYPE_E flag) are placed in the simulation box. You may disable the extension by commenting out \"#define EQUILIBRIUM_BOUNDARIES\" in defines.hpp.");
#endif // EQUILIBRIUM_BOUNDARIES
#ifndef SURFACE
	if(surface_used) print_error("Some nodes are set as fluid/interface/gas with the TYPE_F/TYPE_I/TYPE_G flags, but SURFACE is not enabled. Uncomment \"#define SURFACE\" in defines.hpp.");
#else // SURFACE
	if(!surface_used) print_error("The SURFACE extension is enabled but no fluid/interface/gas nodes (TYPE_F/TYPE_I/TYPE_G flags) are placed in the simulation box. Disable the extension by commenting out \"#define SURFACE\" in defines.hpp.");
#endif // SURFACE
#ifndef TEMPERATURE
	if(temperature_used) print_error("Some nodes are set as temperature boundary with the TYPE_T flag, but TEMPERATURE is not enabled. Uncomment \"#define TEMPERATURE\" in defines.hpp.");
#endif // TEMPERATURE
}

string LBM::default_filename(const string& path, const string& name, const string& extension) {
	string time = "00000000"+to_string(t);
	time = substring(time, length(time)-9u, 9u);
	return create_file_extension((path=="" ? get_exe_path()+"export/" : path)+(name=="" ? "file" : name)+"-"+time, extension);
}
string LBM::default_filename(const string& name, const string& extension) {
	return default_filename("", name, extension);
}

void LBM::allocate(Device& device) {
	const uint N = Nx*Ny*Nz;
	rho = Memory<float>(device, N, 1u, true, true, 1.0f);
	u = Memory<float>(device, N, 3u);
	flags = Memory<uchar>(device, N);
	fi = Memory<fpxx>(device, N, velocity_set, false);
	kernel_initialize = Kernel(device, N, "initialize", fi, rho, u, flags);
	kernel_stream_collide = Kernel(device, N, "stream_collide", fi, rho, u, flags, t, fx, fy, fz);
	kernel_update_fields = Kernel(device, N, "update_fields", fi, rho, u, flags, t, fx, fy, fz);

#ifdef FORCE_FIELD
	F = Memory<float>(device, N, 3u);
	kernel_stream_collide.add_parameters(F);
	kernel_update_fields.add_parameters(F);
	kernel_calculate_force_on_boundaries = Kernel(device, N, "calculate_force_on_boundaries", fi, flags, t, F);
#endif // FORCE_FIELD

#ifdef MOVING_BOUNDARIES
	kernel_update_moving_boundaries = Kernel(device, N, "update_moving_boundaries", u, flags);
#endif // MOVING_BOUNDARIES

#ifdef SURFACE
	phi = Memory<float>(device, N);
	mass = Memory<float>(device, N);
	massex = Memory<float>(device, N);
	kernel_initialize.add_parameters(mass, massex, phi);
	kernel_stream_collide.add_parameters(mass);
	kernel_surface_0 = Kernel(device, N, "surface_0", fi, rho, u, flags, mass, massex, phi, t, fx, fy, fz);
	kernel_surface_1 = Kernel(device, N, "surface_1", flags);
	kernel_surface_2 = Kernel(device, N, "surface_2", fi, rho, u, flags, t);
	kernel_surface_3 = Kernel(device, N, "surface_3", rho, flags, mass, massex, phi);
#endif // SURFACE

#ifdef TEMPERATURE
	T = Memory<float>(device, N, 1u, true, true, 1.0f);
	gi = Memory<fpxx>(device, N, 7u, false);
	kernel_initialize.add_parameters(gi, T);
	kernel_stream_collide.add_parameters(gi, T);
	kernel_update_fields.add_parameters(gi, T);
#endif // TEMPERATURE
}

void LBM::initialize() {
	sanity_checks_initialization();
	rho.write_to_device();
	u.write_to_device();
	flags.write_to_device();
#ifdef FORCE_FIELD
	F.write_to_device();
#endif // FORCE_FIELD
#ifdef SURFACE
	phi.write_to_device();
#endif // SURFACE
#ifdef TEMPERATURE
	T.write_to_device();
#endif // TEMPERATURE
	kernel_initialize.run();
	initialized = true;
}

void LBM::do_time_step() {
#ifdef SURFACE
	kernel_surface_0.set_parameters(7u, t, fx, fy, fz).run();
#endif // SURFACE
	kernel_stream_collide.set_parameters(4u, t, fx, fy, fz).run();
#ifdef SURFACE
	kernel_surface_1.run();
	kernel_surface_2.set_parameters(4u, t).run();
	kernel_surface_3.run();
#endif // SURFACE
	t++; // increment time step
#ifdef UPDATE_FIELDS
	t_last_update_fields = t;
#endif // UPDATE_FIELDS
}

void LBM::run(const ulong steps) { // initializes the LBM simulation (copies data to device and runs initialize kernel), then runs LBM
	info.append(steps, t);
	if(!initialized) {
		initialize(); // only initialize if run() was not called before
		info.print_initialize(); // only print setup info if the setup is new (run() was not called before)
	}
	Clock clock;
	for(ulong i=1ull; i<=steps; i++) { // run LBM in loop, runs infinitely long if steps = max_ulong
#if defined(INTERACTIVE_GRAPHICS)||defined(INTERACTIVE_GRAPHICS_ASCII)
		while(!key_P&&running) sleep(0.016);
		if(!running) break;
#endif // INTERACTIVE_GRAPHICS || INTERACTIVE_GRAPHICS_ASCII
		clock.start();
		do_time_step(); // execute one LBM time step
		info.update(clock.stop());
	}
}

void LBM::update_fields() { // update fields (rho, u, T) manually
#ifndef UPDATE_FIELDS
	if(t>t_last_update_fields) { // only run kernel_update_fields if the time step has changed since last update
		kernel_update_fields.set_parameters(4u, t, fx, fy, fz).run();
		t_last_update_fields = t;
	}
#endif // UPDATE_FIELDS
}

void LBM::reset() { // reset simulation (takes effect in following run() call)
	initialized = false;
}

#ifdef FORCE_FIELD
void LBM::calculate_force_on_boundaries() { // calculate forces from fluid on TYPE_S nodes
	kernel_calculate_force_on_boundaries.set_parameters(2u, t).run();
}
float3 LBM::calculate_force_on_object(const uchar flag_marker) { // add up force for all nodes flagged with flag_marker
	double3 force(0.0, 0.0, 0.0);
	for(uint n=0u; n<get_N(); n++) {
		if(flags[n]==flag_marker) {
			force.x += (double)F.x[n];
			force.y += (double)F.y[n];
			force.z += (double)F.z[n];
		}
	}
	return float3(force.x, force.y, force.z);
}
float3 LBM::calculate_torque_on_object(const uchar flag_marker) { // add up torque around center of mass for all nodes flagged with flag_marker
	double3 center_of_mass(0.0, 0.0, 0.0);
	uint counter = 0u;
	for(uint n=0u; n<get_N(); n++) {
		if(flags[n]==flag_marker) {
			const float3 p = position(n);
			center_of_mass.x += (double)p.x;
			center_of_mass.y += (double)p.y;
			center_of_mass.z += (double)p.z;
			counter++;
		}
	}
	return calculate_torque_on_object(float3(center_of_mass.x/(double)counter, center_of_mass.y/(double)counter, center_of_mass.z/(double)counter)+center(), flag_marker);
}
float3 LBM::calculate_torque_on_object(const float3& rotation_center, const uchar flag_marker) { // add up torque around specified rotation center for all nodes flagged with flag_marker
	double3 torque(0.0, 0.0, 0.0);
	const float3 rotation_center_in_box = rotation_center-center();
	for(uint n=0u; n<get_N(); n++) {
		if(flags[n]==flag_marker) {
			const float3 t = cross(position(n)-rotation_center_in_box, float3(F.x[n], F.y[n], F.z[n]));
			torque.x += (double)t.x;
			torque.y += (double)t.y;
			torque.z += (double)t.z;
		}
	}
	return float3(torque.x, torque.y, torque.z);
}
#endif // FORCE_FIELD

#ifdef MOVING_BOUNDARIES
void LBM::update_moving_boundaries() { // mark/unmark nodes next to TYPE_S nodes with velocity!=0 with TYPE_MS
	kernel_update_moving_boundaries.run();
}
#endif // MOVING_BOUNDARIES

void LBM::write_status(const string& path) { // write LBM status report to a .txt file
	string status = "";
	status += "Grid Resolution = ("+to_string(Nx)+", "+to_string(Ny)+", "+to_string(Nz)+")\n";
	status += "LBM type = D"+string(velocity_set==9 ? "2" : "3")+"Q"+to_string(velocity_set)+" "+info.collision+"\n";
	status += "Memory Usage = "+to_string(info.cpu_mem_required)+" MB (CPU), "+to_string(info.gpu_mem_required)+" MB (GPU)\n";
	status += "Maximum Allocation Size = "+to_string((uint)((ulong)get_N()*(ulong)(velocity_set*sizeof(fpxx))/1048576ull))+" MB\n";
	status += "Time Step = "+to_string(t)+" / "+(info.steps==max_ulong ? "infinite" : to_string(info.steps))+"\n";
	status += "Kinematic Viscosity = "+to_string(nu)+"\n";
	status += "Relaxation Time = "+to_string(get_tau())+"\n";
	status += "Maximum Reynolds Number = "+to_string(get_Re_max())+"\n";
#ifdef VOLUME_FORCE
	status += "Volume Force = ("+to_string(fx)+", "+to_string(fy)+", "+to_string(fz)+")\n";
#endif // VOLUME_FORCE
#ifdef SURFACE
	status += "Surface Tension Coefficient = "+to_string(sigma)+"\n";
#endif // SURFACE
#ifdef TEMPERATURE
	status += "Thermal Diffusion Coefficient = "+to_string(alpha)+"\n";
	status += "Thermal Expansion Coefficient = "+to_string(beta)+"\n";
#endif // TEMPERATURE
	const string filename = default_filename(path, "status", ".txt");
	write_file(filename, status);
}

template<typename T> static string vtk_type() {
	     if constexpr(std::is_same<T,  char >::value) return "char";
	else if constexpr(std::is_same<T, uchar >::value) return "unsigned_char";
	else if constexpr(std::is_same<T,  short>::value) return "short";
	else if constexpr(std::is_same<T, ushort>::value) return "unsigned_short";
	else if constexpr(std::is_same<T,  int  >::value) return "int";
	else if constexpr(std::is_same<T, uint  >::value) return "unsigned_int";
	else if constexpr(std::is_same<T, slong >::value) return "long";
	else if constexpr(std::is_same<T, ulong >::value) return "unsigned_long";
	else if constexpr(std::is_same<T, float >::value) return "float";
	else if constexpr(std::is_same<T, double>::value) return "double";
	else print_error("Error in vtk_type(): Type not supported.");
	return "";
}
template<typename T> void write_vtk(const string& path, const Memory<T>& memory, const uint Nx, const uint Ny, const uint Nz) { // write binary .vtk file
	const float spacing = units.si_x(1.0f);
	const float3 origin = spacing*float3(0.5f-0.5f*(float)Nx, 0.5f-0.5f*(float)Ny, 0.5f-0.5f*(float)Nz);
	const string header =
		"# vtk DataFile Version 3.0\nData\nBINARY\nDATASET STRUCTURED_POINTS\n"
		"DIMENSIONS "+to_string(Nx)+" "+to_string(Ny)+" "+to_string(Nz)+"\n"
		"ORIGIN "+to_string(origin.x)+" "+to_string(origin.y)+" "+to_string(origin.z)+"\n"
		"SPACING "+to_string(spacing)+" "+to_string(spacing)+" "+to_string(spacing)+"\n"
		"POINT_DATA "+to_string(Nx*Ny*Nz)+"\nSCALARS data "+vtk_type<T>()+" "+to_string(memory.dimensions())+"\nLOOKUP_TABLE default\n"
	;
	T* data = new T[memory.range()];
	for(uint d=0u; d<memory.dimensions(); d++) {
		for(ulong i=0ull; i<memory.length(); i++) {
			data[i*(ulong)memory.dimensions()+(ulong)d] = reverse_bytes(memory(i, d)); // SoA <- AoS
		}
	}
	create_folder(path);
	const string filename = create_file_extension(path, ".vtk");
	std::ofstream file(filename, std::ios::out|std::ios::binary);
	file.write(header.c_str(), header.length()); // write non-binary file header
	file.write((char*)data, memory.capacity()); // write binary data
	file.close();
	delete[] data;
	info.allow_rendering = false; // temporarily disable interactive rendering
	print_info("File \""+path+"\" saved.");
	info.allow_rendering = true;
}

void LBM::rho_write_host_to_vtk(const string& path) {
	const string filename = default_filename(path, "rho", ".vtk");
	write_vtk(filename, rho, Nx, Ny, Nz);
}
void LBM::rho_write_device_to_vtk(const string& path) {
#ifndef UPDATE_FIELDS
	update_fields();
#endif // UPDATE_FIELDS
	rho.read_from_device();
	rho_write_host_to_vtk(path);
}
void LBM::u_write_host_to_vtk(const string& path) {
	const string filename = default_filename(path, "u", ".vtk");
	write_vtk(filename, u, Nx, Ny, Nz);
}
void LBM::u_write_device_to_vtk(const string& path) {
#ifndef UPDATE_FIELDS
	update_fields();
#endif // UPDATE_FIELDS
	u.read_from_device();
	u_write_host_to_vtk(path);
}
void LBM::flags_write_host_to_vtk(const string& path) {
	const string filename = default_filename(path, "flags", ".vtk");
	write_vtk(filename, flags, Nx, Ny, Nz);
}
void LBM::flags_write_device_to_vtk(const string& path) {
	flags.read_from_device();
	flags_write_host_to_vtk(path);
}

#ifdef FORCE_FIELD
void LBM::F_write_host_to_vtk(const string& path) {
	const string filename = default_filename(path, "F", ".vtk");
	write_vtk(filename, F, Nx, Ny, Nz);
}
void LBM::F_write_device_to_vtk(const string& path) {
	F.read_from_device();
	F_write_host_to_vtk(path);
}
#endif // FORCE_FIELD

#ifdef SURFACE
void LBM::phi_write_host_to_vtk(const string& path) {
	const string filename = default_filename(path, "phi", ".vtk");
	write_vtk(filename, phi, Nx, Ny, Nz);
}
void LBM::phi_write_device_to_vtk(const string& path) {
	phi.read_from_device();
	phi_write_host_to_vtk(path);
}
#endif // SURFACE

#ifdef TEMPERATURE
void LBM::T_write_host_to_vtk(const string& path) {
	const string filename = default_filename(path, "T", ".vtk");
	write_vtk(filename, T, Nx, Ny, Nz);
}
void LBM::T_write_device_to_vtk(const string& path) {
#ifndef UPDATE_FIELDS
	update_fields();
#endif // UPDATE_FIELDS
	T.read_from_device();
	T_write_host_to_vtk(path);
}
#endif // TEMPERATURE

void LBM::voxelize_mesh(const Mesh* mesh, const uchar flag) { // voxelize triangle mesh
	print_info("Voxelizing mesh. This may take a few minutes.");
	Memory<float3> p0(device, mesh->triangle_number, 1u, mesh->p0);
	Memory<float3> p1(device, mesh->triangle_number, 1u, mesh->p1);
	Memory<float3> p2(device, mesh->triangle_number, 1u, mesh->p2);
	const float x0=mesh->pmin.x, y0=mesh->pmin.y, z0=mesh->pmin.z, x1=mesh->pmax.x, y1=mesh->pmax.y, z1=mesh->pmax.z; // use bounding box of mesh to speed up voxelization
	Kernel kernel_voxelize_mesh(device, get_N(), "voxelize_mesh", flags, flag, p0, p1, p2, mesh->triangle_number, x0, y0, z0, x1, y1, z1);
	p0.write_to_device();
	p1.write_to_device();
	p2.write_to_device();
	flags.write_to_device();
	kernel_voxelize_mesh.run();
	flags.read_from_device();
}
void LBM::voxelize_stl(const string& path, const float3& center, const float3x3& rotation, const float size, const uchar flag) { // voxelize triangle mesh
	const Mesh* mesh = read_stl(path, float3(get_Nx(), get_Ny(), get_Nz()), center, rotation, size);
	voxelize_mesh(mesh, flag);
	delete mesh;
}
void LBM::voxelize_stl(const string& path, const float3x3& rotation, const float size, const uchar flag) { // read and voxelize binary .stl file (place in box center)
	voxelize_stl(path, center(), rotation, size, flag);
}
void LBM::voxelize_stl(const string& path, const float3& center, const float size, const uchar flag) { // read and voxelize binary .stl file (no rotation)
	voxelize_stl(path, center, float3x3(1.0f), size, flag);
}
void LBM::voxelize_stl(const string& path, const float size, const uchar flag) { // read and voxelize binary .stl file (place in box center, no rotation)
	voxelize_stl(path, center(), float3x3(1.0f), size, flag);
}

string LBM::device_defines() const { return
	"\n	#define def_Nx "+to_string(Nx)+"u"
	"\n	#define def_Ny "+to_string(Ny)+"u"
	"\n	#define def_Nz "+to_string(Nz)+"u"
	"\n	#define def_N "+to_string(get_N())+"ul"

	"\n	#define D"+to_string(dimensions)+"Q"+to_string(velocity_set)+"" // D2Q9/D3Q15/D3Q19/D3Q27
	"\n	#define def_velocity_set "+to_string(velocity_set)+"u" // LBM velocity set (D2Q9/D3Q15/D3Q19/D3Q27)
	"\n	#define def_dimensions "+to_string(dimensions)+"u" // number spatial dimensions (2D or 3D)

	"\n	#define def_c 0.57735027f" // lattice speed of sound c = 1/sqrt(3)*dt
	"\n	#define def_w " +to_string(1.0f/get_tau())+"f" // relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)
#if defined(D2Q9)
	"\n	#define def_w0 (1.0f/2.25f)" // center (0)
	"\n	#define def_ws (1.0f/9.0f)" // straight (1-4)
	"\n	#define def_we (1.0f/36.0f)" // edge (5-8)
#elif defined(D3Q15)
	"\n	#define def_w0 (1.0f/4.5f)" // center (0)
	"\n	#define def_ws (1.0f/9.0f)" // straight (1-6)
	"\n	#define def_wc (1.0f/72.0f)" // corner (7-14)
#elif defined(D3Q19)
	"\n	#define def_w0 (1.0f/3.0f)" // center (0)
	"\n	#define def_ws (1.0f/18.0f)" // straight (1-6)
	"\n	#define def_we (1.0f/36.0f)" // edge (7-18)
#elif defined(D3Q27)
	"\n	#define def_w0 (1.0f/3.375f)" // center (0)
	"\n	#define def_ws (1.0f/13.5f)" // straight (1-6)
	"\n	#define def_we (1.0f/54.0f)" // edge (7-18)
	"\n	#define def_wc (1.0f/216.0f)" // corner (19-26)
#endif // D3Q27

#if defined(SRT)
	"\n	#define SRT"
#elif defined(TRT)
	"\n	#define TRT"
#endif // TRT

	"\n	#define TYPE_S 0x01" // 0b00000001 // (stationary or moving) solid boundary
	"\n	#define TYPE_E 0x02" // 0b00000010 // equilibrium boundary (inflow/outflow)
	"\n	#define TYPE_T 0x04" // 0b00000100 // temperature boundary
	"\n	#define TYPE_F 0x08" // 0b00001000 // fluid
	"\n	#define TYPE_I 0x10" // 0b00010000 // interface
	"\n	#define TYPE_G 0x20" // 0b00100000 // gas
	"\n	#define TYPE_X 0x40" // 0b01000000 // reserved type X
	"\n	#define TYPE_Y 0x80" // 0b10000000 // reserved type Y

	"\n	#define TYPE_MS 0x03" // 0b00000011 // node next to moving solid boundary
	"\n	#define TYPE_BO 0x03" // 0b00000011 // any flag bit used for boundaries (temperature excluded)
	"\n	#define TYPE_IF 0x18" // 0b00011000 // change from interface to fluid
	"\n	#define TYPE_IG 0x30" // 0b00110000 // change from interface to gas
	"\n	#define TYPE_GI 0x38" // 0b00111000 // change from gas to interface
	"\n	#define TYPE_SU 0x38" // 0b00111000 // any flag bit used for SURFACE

#if defined(FP16S)
	"\n	#define fpxx half" // switchable data type (scaled IEEE-754 16-bit floating-point format: 1-5-10, exp-30, +-1.99902344, +-1.86446416E-9, +-1.81898936E-12, 3.311 digits)
	"\n	#define load(p,o) vload_half(o,p)*3.0517578E-5f" // special function for loading half
	"\n	#define store(p,o,x) vstore_half_rte((x)*32768.0f,o,p)" // special function for storing half
#elif defined(FP16C)
	"\n	#define fpxx ushort" // switchable data type (custom 16-bit floating-point format: 1-4-11, exp-15, +-1.99951168, +-6.10351562E-5, +-2.98023224E-8, 3.612 digits), 12.5% slower than IEEE-754 16-bit
	"\n	#define load(p,o) half_to_float_custom(p[o])" // special function for loading half
	"\n	#define store(p,o,x) p[o]=float_to_half_custom(x)" // special function for storing half
#else // FP32
	"\n	#define fpxx float" // switchable data type (regular 32-bit float)
	"\n	#define load(p,o) p[o]" // regular float read
	"\n	#define store(p,o,x) p[o]=x" // regular float write
#endif // FP32

#ifdef UPDATE_FIELDS
	"\n	#define UPDATE_FIELDS"
#endif // UPDATE_FIELDS

#ifdef VOLUME_FORCE
	"\n	#define VOLUME_FORCE"
#endif // VOLUME_FORCE

#ifdef MOVING_BOUNDARIES
	"\n	#define MOVING_BOUNDARIES"
#endif // MOVING_BOUNDARIES

#ifdef EQUILIBRIUM_BOUNDARIES
	"\n	#define EQUILIBRIUM_BOUNDARIES"
#endif // EQUILIBRIUM_BOUNDARIES

#ifdef FORCE_FIELD
	"\n	#define FORCE_FIELD"
#endif // FORCE_FIELD

#ifdef SURFACE
	"\n	#define SURFACE"
	"\n	#define def_6_sigma "+to_string(6.0f*sigma)+"f" // rho_laplace = 2*o*K, rho = 1-rho_laplace/c^2 = 1-(6*o)*K
#endif // SURFACE

#ifdef TEMPERATURE
	"\n	#define TEMPERATURE"
	"\n	#define def_w_T "+to_string(1.0f/(2.0f*alpha+0.5f))+"f" // wT = dt/tauT = 1/(2*alpha+1/2), alpha = thermal diffusion coefficient
	"\n	#define def_beta "+to_string(beta)+"f" // thermal expansion coefficient
	"\n	#define def_T_avg "+to_string(T_avg)+"f" // average temperature
#endif // TEMPERATURE

#ifdef SUBGRID
	"\n	#define SUBGRID"
#endif // SUBGRID
;}

#ifdef GRAPHICS
void LBM::Graphics::default_settings() {
	key_1 = true;
}

void LBM::Graphics::allocate(Device& device) {
	bitmap = Memory<int>(device, camera.width*camera.height);
	zbuffer = Memory<int>(device, camera.width*camera.height, 1u, false);
	camera_parameters = Memory<float>(device, 15u);
	kernel_clear = Kernel(device, bitmap.length(), "graphics_clear", bitmap, zbuffer);

	camera.set_zoom(0.5f*(float)fmax(fmax(lbm->get_Nx(), lbm->get_Ny()), lbm->get_Nz()));
	default_settings();

	kernel_graphics_flags = Kernel(device, lbm->get_N(), "graphics_flags", lbm->flags, camera_parameters, bitmap, zbuffer);
	kernel_graphics_field = Kernel(device, lbm->get_N(), "graphics_field", lbm->flags, lbm->u, camera_parameters, bitmap, zbuffer);
#ifndef D2Q9
	kernel_graphics_streamline = Kernel(device, (lbm->get_Nx()/GRAPHICS_STREAMLINE_SPARSE)*(lbm->get_Ny()/GRAPHICS_STREAMLINE_SPARSE)*(lbm->get_Nz()/GRAPHICS_STREAMLINE_SPARSE), "graphics_streamline", lbm->flags, lbm->u, camera_parameters, bitmap, zbuffer); // 3D
#else // D2Q9
	kernel_graphics_streamline = Kernel(device, (lbm->get_Nx()/GRAPHICS_STREAMLINE_SPARSE)*(lbm->get_Ny()/GRAPHICS_STREAMLINE_SPARSE), "graphics_streamline", lbm->flags, lbm->u, camera_parameters, bitmap, zbuffer); // 2D
#endif // D2Q9
	kernel_graphics_q = Kernel(device, lbm->get_N(), "graphics_q", lbm->flags, lbm->u, camera_parameters, bitmap, zbuffer);

#ifdef FORCE_FIELD
	kernel_graphics_flags.add_parameters(lbm->F);
#endif // FORCE_FIELD

#ifdef SURFACE
	skybox = Memory<int>(device, skybox_image->width()*skybox_image->height(), 1u, skybox_image->data());
	kernel_graphics_rasterize_phi = Kernel(device, lbm->get_N(), "graphics_rasterize_phi", lbm->phi, camera_parameters, bitmap, zbuffer);
	kernel_graphics_raytrace_phi = Kernel(device, bitmap.length(), "graphics_raytrace_phi", lbm->phi, lbm->flags, skybox, camera_parameters, bitmap);
#endif // SURFACE

#ifdef TEMPERATURE
	kernel_graphics_streamline.add_parameters(lbm->T);
#endif // TEMPERATURE
}
bool LBM::Graphics::update_camera() {
	camera.update_matrix();
	bool change = false;
	for(uint i=0u; i<15u; i++) {
		const float data = camera.data(i);
		change |= (camera_parameters[i]!=data);
		camera_parameters[i] = data;
	}
	return change; // return false if camera parameters remain unchanged
}
int* LBM::Graphics::draw_frame() {
	const bool camera_update = update_camera();
#if defined(INTERACTIVE_GRAPHICS)||defined(INTERACTIVE_GRAPHICS_ASCII)
	if(!camera_update&&!camera.key_update&&lbm->get_t()==t_last_frame) return bitmap.data(); // don't render a new frame if the scene hasn't changed since last frame
#endif // INTERACTIVE_GRAPHICS or INTERACTIVE_GRAPHICS_ASCII
	t_last_frame = lbm->get_t();
#ifndef UPDATE_FIELDS
	if(key_2||key_3||key_4) lbm->update_fields(); // only call update_fields() if the time step has changed since the last rendered frame
#endif // UPDATE_FIELDS
	camera.key_update = false;

	if(camera_update) camera_parameters.write_to_device(false); // camera_parameters PCIe transfer and kernel_clear execution can happen simulataneously
	kernel_clear.run();
#ifdef SURFACE
	if(key_6) kernel_graphics_raytrace_phi.run();
	if(key_5) kernel_graphics_rasterize_phi.run();
#endif // SURFACE
	if(key_1) kernel_graphics_flags.run();
	if(key_2) kernel_graphics_field.run();
	if(key_3) kernel_graphics_streamline.run();
	if(key_4) kernel_graphics_q.run();

	bitmap.read_from_device();
	return bitmap.data();
}
string LBM::Graphics::device_defines() const { return
	"\n	#define GRAPHICS"
	"\n	#define def_background_color " +to_string(GRAPHICS_BACKGROUND_COLOR)+""
	"\n	#define def_screen_width "     +to_string(camera.width)+"u"
	"\n	#define def_screen_height "    +to_string(camera.height)+"u"
	"\n	#define def_n "                +to_string(1.333f)+"f" // refractive index of water
	"\n	#define def_scale_u "          +to_string(0.57735027f/(GRAPHICS_U_MAX))+"f"
	"\n	#define def_scale_Q_min "      +to_string(GRAPHICS_Q_CRITERION)+"f"
	"\n	#define def_scale_F "          +to_string(GRAPHICS_BOUNDARY_FORCE_SCALE)+"f"
	"\n	#define def_streamline_sparse "+to_string(GRAPHICS_STREAMLINE_SPARSE)+"u"
	"\n	#define def_streamline_length "+to_string(GRAPHICS_STREAMLINE_LENGTH)+"u"

	"\n	#define COLOR_S (127<<16|127<<8|127)" // coloring scheme
	"\n	#define COLOR_T (255<<16|  0<<8|  0)"
	"\n	#define COLOR_E (  0<<16|255<<8|  0)"
	"\n	#define COLOR_M (255<<16|  0<<8|255)"
	"\n	#define COLOR_F (  0<<16|  0<<8|255)"
	"\n	#define COLOR_I (  0<<16|255<<8|255)"
	"\n	#define COLOR_X (255<<16|127<<8|  0)"
	"\n	#define COLOR_Y (255<<16|255<<8|  0)"
	"\n	#define COLOR_0 (127<<16|127<<8|127)"

#ifndef SURFACE
	"\n	#define def_skybox_width 1u"
	"\n	#define def_skybox_height 1u"
#else // SURFACE
	"\n	#define def_skybox_width " +to_string(skybox_image->width() )+"u"
	"\n	#define def_skybox_height "+to_string(skybox_image->height())+"u"
#endif // SURFACE

#ifdef TEMPERATURE
	"\n	#define GRAPHICS_TEMPERATURE"
#endif // TEMPERATURE
;}

void LBM::Graphics::set_camera_centered(const float rx, const float ry, const float fov, const float zoom) {
	camera.free = false;
	camera.rx = 0.5*pi+((double)rx*pi/180.0);
	camera.ry = pi-((double)ry*pi/180.0);
	camera.fov = clamp((float)fov, 1E-6f, 179.0f);
	camera.set_zoom(0.5f*(float)fmax(fmax(lbm->get_Nx(), lbm->get_Ny()), lbm->get_Nz())/zoom);
}
void LBM::Graphics::set_camera_free(const float3& p, const float rx, const float ry, const float fov) {
	camera.free = true;
	camera.rx = 0.5*pi+((double)rx*pi/180.0);
	camera.ry = pi-((double)ry*pi/180.0);
	camera.fov = clamp((float)fov, 1E-6f, 179.0f);
	camera.zoom = 1E16f;
	camera.pos = p;
}
void LBM::Graphics::print_frame() { // preview current frame in console
#ifndef INTERACTIVE_GRAPHICS_ASCII
	info.allow_rendering = false; // temporarily disable interactive rendering
	draw_frame(); // make sure the frame is fully rendered
	Image* image = new Image(camera.width, camera.height, bitmap.data());
	println();
	print_image(image);
	delete image;
	info.allow_rendering = true;
#endif // INTERACTIVE_GRAPHICS_ASCII
}
void encode_image(Image* image, const string& filename, const string& extension, std::atomic_int* running_encoders) {
	if(extension==".png") write_png(filename, image);
	if(extension==".qoi") write_qoi(filename, image);
	if(extension==".bmp") write_bmp(filename, image);
	delete image; // delete image when done
	(*running_encoders)--;
}
void LBM::Graphics::write_frame(const string& path, const string& name, const string& extension, bool print_preview) { // save current frame as .png file (smallest file size, but slow)
	write_frame(0u, 0u, camera.width, camera.height, path, name, extension, print_preview);
}
void LBM::Graphics::write_frame(const uint x1, const uint y1, const uint x2, const uint y2, const string& path, const string& name, const string& extension, bool print_preview) { // save a cropped current frame with two corner points (x1,y1) and (x2,y2)
	info.allow_rendering = false; // temporarily disable interactive rendering
	draw_frame(); // make sure the frame is fully rendered
	const string filename = lbm->default_filename(path, name, extension);
	const uint xa=max(min(x1, x2), 0u), xb=min(max(x1, x2), camera.width ); // sort coordinates if necessary
	const uint ya=max(min(y1, y2), 0u), yb=min(max(y1, y2), camera.height);
	Image* image = new Image(xb-xa, yb-ya); // create local copy of frame buffer
	for(uint y=0u; y<image->height(); y++) for(uint x=0u; x<image->width(); x++) image->set_color(x, y, bitmap[camera.width*(ya+y)+(xa+x)]);
#ifndef INTERACTIVE_GRAPHICS_ASCII
	if(print_preview) {
		println();
		print_image(image);
		print_info("Image \""+filename+"\" saved.");
	}
#endif // INTERACTIVE_GRAPHICS_ASCII
	running_encoders++;
	thread encoder(encode_image, image, filename, extension, &running_encoders); // the main bottleneck in rendering images to the hard disk is .png encoding, so encode image in new thread
	encoder.detach(); // detatch thread so it can run concurrently
	info.allow_rendering = true;
}
void LBM::Graphics::write_frame_png(const string& path, bool print_preview) { // save current frame as .png file (smallest file size, but slow)
	write_frame(path, "image", ".png", print_preview);
}
void LBM::Graphics::write_frame_qoi(const string& path, bool print_preview) { // save current frame as .qoi file (small file size, fast)
	write_frame(path, "image", ".qoi", print_preview);
}
void LBM::Graphics::write_frame_bmp(const string& path, bool print_preview) { // save current frame as .bmp file (large file size, fast)
	write_frame(path, "image", ".bmp", print_preview);
}
void LBM::Graphics::write_frame_png(const uint x1, const uint y1, const uint x2, const uint y2, const string& path, bool print_preview) { // save current frame as .png file (smallest file size, but slow)
	write_frame(x1, y1, x2, y2, path, "image", ".png", print_preview);
}
void LBM::Graphics::write_frame_qoi(const uint x1, const uint y1, const uint x2, const uint y2, const string& path, bool print_preview) { // save current frame as .qoi file (small file size, fast)
	write_frame(x1, y1, x2, y2, path, "image", ".qoi", print_preview);
}
void LBM::Graphics::write_frame_bmp(const uint x1, const uint y1, const uint x2, const uint y2, const string& path, bool print_preview) { // save current frame as .bmp file (large file size, fast)
	write_frame(x1, y1, x2, y2, path, "image", ".bmp", print_preview);
}
#endif // GRAPHICS