#include "lbm.hpp"
#include "graphics.hpp"



Units units; // for unit conversion

#if defined(D2Q9)
const uint velocity_set = 9u;
const uint dimensions = 2u;
const uint transfers = 3u;
#elif defined(D3Q15)
const uint velocity_set = 15u;
const uint dimensions = 3u;
const uint transfers = 5u;
#elif defined(D3Q19)
const uint velocity_set = 19u;
const uint dimensions = 3u;
const uint transfers = 5u;
#elif defined(D3Q27)
const uint velocity_set = 27u;
const uint dimensions = 3u;
const uint transfers = 9u;
#endif // D3Q27

string default_filename(const string& path, const string& name, const string& extension, const ulong t) { // generate a default filename with timestamp
	string time = "00000000"+to_string(t);
	time = substring(time, length(time)-9u, 9u);
	return create_file_extension((path=="" ? get_exe_path()+"export/" : path)+(name=="" ? "file" : name)+"-"+time, extension);
}
string default_filename(const string& name, const string& extension, const ulong t) { // generate a default filename with timestamp at exe_path/export/
	return default_filename("", name, extension, t);
}



LBM_Domain::LBM_Domain(const Device_Info& device_info, const uint Nx, const uint Ny, const uint Nz, const uint Dx, const uint Dy, const uint Dz, const int Ox, const int Oy, const int Oz, const float nu, const float fx, const float fy, const float fz, const float sigma, const float alpha, const float beta, const uint particles_N, const float particles_rho) { // constructor with manual device selection and domain offset
	this->Nx = Nx; this->Ny = Ny; this->Nz = Nz;
	this->Dx = Dx; this->Dy = Dy; this->Dz = Dz;
	this->Ox = Ox; this->Oy = Oy; this->Oz = Oz;
	this->nu = nu;
	this->fx = fx; this->fy = fy; this->fz = fz;
	this->sigma = sigma;
	this->alpha = alpha; this->beta = beta;
	this->particles_N = particles_N;
	this->particles_rho = particles_rho;
	string opencl_c_code;
#ifdef GRAPHICS
	graphics = Graphics(this);
	opencl_c_code = device_defines()+graphics.device_defines()+get_opencl_c_code();
#else // GRAPHICS
	opencl_c_code = device_defines()+get_opencl_c_code();
#endif // GRAPHICS
	this->device = Device(device_info, opencl_c_code);
	allocate(device); // lbm first
#ifdef GRAPHICS
	graphics.allocate(device); // graphics after lbm
#endif // GRAPHICS
}

void LBM_Domain::allocate(Device& device) {
	const ulong N = get_N();
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
	kernel_reset_force_field = Kernel(device, N, "reset_force_field", F);
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

#ifdef PARTICLES
	particles = Memory<float>(device, (ulong)particles_N, 3u);
	kernel_integrate_particles = Kernel(device, (ulong)particles_N, "integrate_particles", particles, u, flags, 1.0f);
#ifdef FORCE_FIELD
	kernel_integrate_particles.add_parameters(F, fx, fy, fz);
#endif // FORCE_FIELD
#endif // PARTICLES

	if(get_D()>1u) allocate_transfer(device);
}

void LBM_Domain::enqueue_initialize() { // call kernel_initialize
	kernel_initialize.enqueue_run();
}
void LBM_Domain::enqueue_stream_collide() { // call kernel_stream_collide to perform one LBM time step
	kernel_stream_collide.set_parameters(4u, t, fx, fy, fz).enqueue_run();
}
void LBM_Domain::enqueue_update_fields() { // update fields (rho, u, T) manually
#ifndef UPDATE_FIELDS
	if(t!=t_last_update_fields) { // only run kernel_update_fields if the time step has changed since last update
		kernel_update_fields.set_parameters(4u, t, fx, fy, fz).enqueue_run();
		t_last_update_fields = t;
	}
#endif // UPDATE_FIELDS
}
#ifdef SURFACE
void LBM_Domain::enqueue_surface_0() {
	kernel_surface_0.set_parameters(7u, t, fx, fy, fz).enqueue_run();
}
void LBM_Domain::enqueue_surface_1() {
	kernel_surface_1.enqueue_run();
}
void LBM_Domain::enqueue_surface_2() {
	kernel_surface_2.set_parameters(4u, t).enqueue_run();
}
void LBM_Domain::enqueue_surface_3() {
	kernel_surface_3.enqueue_run();
}
#endif // SURFACE
#ifdef FORCE_FIELD
void LBM_Domain::enqueue_calculate_force_on_boundaries() { // calculate forces from fluid on TYPE_S nodes
	kernel_calculate_force_on_boundaries.set_parameters(2u, t).enqueue_run();
}
#endif // FORCE_FIELD
#ifdef MOVING_BOUNDARIES
void LBM_Domain::enqueue_update_moving_boundaries() { // mark/unmark nodes next to TYPE_S nodes with velocity!=0 with TYPE_MS
	kernel_update_moving_boundaries.enqueue_run();
}
#endif // MOVING_BOUNDARIES
#ifdef PARTICLES
void LBM_Domain::enqueue_integrate_particles(const uint time_step_multiplicator) { // intgegrate particles forward in time and couple particles to fluid
#ifdef FORCE_FIELD
	if(particles_rho!=1.0f) kernel_reset_force_field.enqueue_run(); // only reset force field if particles have buoyancy and apply forces on fluid
	kernel_integrate_particles.set_parameters(5u, fx, fy, fz);
#endif // FORCE_FIELD
	kernel_integrate_particles.set_parameters(3u, (float)time_step_multiplicator).enqueue_run();
}
#endif // PARTICLES

void LBM_Domain::increment_time_step(const uint steps) {
	t += (ulong)steps; // increment time step
#ifdef UPDATE_FIELDS
	t_last_update_fields = t;
#endif // UPDATE_FIELDS
}
void LBM_Domain::reset_time_step() {
	t = 0ull; // increment time step
#ifdef UPDATE_FIELDS
	t_last_update_fields = t;
#endif // UPDATE_FIELDS
}
void LBM_Domain::finish_queue() {
	device.finish_queue();
}

uint LBM_Domain::get_velocity_set() const {
	return velocity_set;
}

void LBM_Domain::voxelize_mesh_on_device(const Mesh* mesh, const uchar flag, const float3& rotation_center, const float3& linear_velocity, const float3& rotational_velocity) { // voxelize triangle mesh
	Memory<float3> p0(device, mesh->triangle_number, 1u, mesh->p0);
	Memory<float3> p1(device, mesh->triangle_number, 1u, mesh->p1);
	Memory<float3> p2(device, mesh->triangle_number, 1u, mesh->p2);
	Memory<float> bounding_box_and_velocity(device, 16u);
	const float x0=mesh->pmin.x-2.0f, y0=mesh->pmin.y-2.0f, z0=mesh->pmin.z-2.0f, x1=mesh->pmax.x+2.0f, y1=mesh->pmax.y+2.0f, z1=mesh->pmax.z+2.0f; // use bounding box of mesh to speed up voxelization; add tolerance of 2 cells for re-voxelization of moving objects
	bounding_box_and_velocity[ 0] = as_float(mesh->triangle_number);
	bounding_box_and_velocity[ 1] = x0;
	bounding_box_and_velocity[ 2] = y0;
	bounding_box_and_velocity[ 3] = z0;
	bounding_box_and_velocity[ 4] = x1;
	bounding_box_and_velocity[ 5] = y1;
	bounding_box_and_velocity[ 6] = z1;
	bounding_box_and_velocity[ 7] = rotation_center.x;
	bounding_box_and_velocity[ 8] = rotation_center.y;
	bounding_box_and_velocity[ 9] = rotation_center.z;
	bounding_box_and_velocity[10] = linear_velocity.x;
	bounding_box_and_velocity[11] = linear_velocity.y;
	bounding_box_and_velocity[12] = linear_velocity.z;
	bounding_box_and_velocity[13] = rotational_velocity.x;
	bounding_box_and_velocity[14] = rotational_velocity.y;
	bounding_box_and_velocity[15] = rotational_velocity.z;
	uint direction = 0u;
	if(length(rotational_velocity)==0.0f) { // choose direction of minimum bounding-box cross-section area
		float v[3] = { (y1-y0)*(z1-z0), (z1-z0)*(x1-x0), (x1-x0)*(y1-y0) };
		float vmin = v[0];
		for(uint i=1u; i<3u; i++) {
			if(v[i]<vmin) {
				vmin = v[i];
				direction = i;
			}
		}
	} else { // choose direction closest to rotation axis
		float v[3] = { rotational_velocity.x, rotational_velocity.y, rotational_velocity.z };
		float vmax = v[0];
		for(uint i=1u; i<3u; i++) {
			if(v[i]>vmax) {
				vmax = v[i];
				direction = i; // find direction of minimum bounding-box cross-section area
			}
		}
	}
	const ulong A[3] = { (ulong)Ny*(ulong)Nz, (ulong)Nz*(ulong)Nx, (ulong)Nx*(ulong)Ny };
	Kernel kernel_voxelize_mesh(device, A[direction], "voxelize_mesh", direction, fi, u, flags, t+1ull, flag, p0, p1, p2, bounding_box_and_velocity);
	p0.write_to_device();
	p1.write_to_device();
	p2.write_to_device();
	bounding_box_and_velocity.write_to_device();
	kernel_voxelize_mesh.run();
}
void LBM_Domain::enqueue_unvoxelize_mesh_on_device(const Mesh* mesh, const uchar flag) { // remove voxelized triangle mesh from LBM grid
	const float x0=mesh->pmin.x, y0=mesh->pmin.y, z0=mesh->pmin.z, x1=mesh->pmax.x, y1=mesh->pmax.y, z1=mesh->pmax.z; // remove all flags in bounding box of mesh
	Kernel kernel_unvoxelize_mesh(device, get_N(), "unvoxelize_mesh", flags, flag, x0, y0, z0, x1, y1, z1);
	kernel_unvoxelize_mesh.run();
}

string LBM_Domain::device_defines() const { return
	"\n	#define def_Nx "+to_string(Nx)+"u"
	"\n	#define def_Ny "+to_string(Ny)+"u"
	"\n	#define def_Nz "+to_string(Nz)+"u"
	"\n	#define def_N "+to_string(get_N())+"ul"

	"\n	#define def_Dx "+to_string(Dx)+"u"
	"\n	#define def_Dy "+to_string(Dy)+"u"
	"\n	#define def_Dz "+to_string(Dz)+"u"

	"\n	#define def_Ox "+to_string(Ox)+"" // offsets are signed integer!
	"\n	#define def_Oy "+to_string(Oy)+""
	"\n	#define def_Oz "+to_string(Oz)+""

	"\n	#define def_Ax "+to_string(Ny*Nz)+"u"
	"\n	#define def_Ay "+to_string(Nz*Nx)+"u"
	"\n	#define def_Az "+to_string(Nx*Ny)+"u"

	"\n	#define def_domain_offset_x "+to_string((float)Ox+(float)(Dx>1u)-0.5f*((float)Dx-1.0f)*(float)(Nx-2u*(Dx>1u)))+"f"
	"\n	#define def_domain_offset_y "+to_string((float)Oy+(float)(Dy>1u)-0.5f*((float)Dy-1.0f)*(float)(Ny-2u*(Dy>1u)))+"f"
	"\n	#define def_domain_offset_z "+to_string((float)Oz+(float)(Dz>1u)-0.5f*((float)Dz-1.0f)*(float)(Nz-2u*(Dz>1u)))+"f"

	"\n	#define D"+to_string(dimensions)+"Q"+to_string(velocity_set)+"" // D2Q9/D3Q15/D3Q19/D3Q27
	"\n	#define def_velocity_set "+to_string(velocity_set)+"u" // LBM velocity set (D2Q9/D3Q15/D3Q19/D3Q27)
	"\n	#define def_dimensions "+to_string(dimensions)+"u" // number spatial dimensions (2D or 3D)
	"\n	#define def_transfers "+to_string(transfers)+"u" // number of DDFs that are transferred between multiple domains

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
	"\n	#define fpxx_copy ushort" // switchable data type for direct copying (scaled IEEE-754 16-bit floating-point format: 1-5-10, exp-30, +-1.99902344, +-1.86446416E-9, +-1.81898936E-12, 3.311 digits)
	"\n	#define load(p,o) vload_half(o,p)*3.0517578E-5f" // special function for loading half
	"\n	#define store(p,o,x) vstore_half_rte((x)*32768.0f,o,p)" // special function for storing half
#elif defined(FP16C)
	"\n	#define fpxx ushort" // switchable data type (custom 16-bit floating-point format: 1-4-11, exp-15, +-1.99951168, +-6.10351562E-5, +-2.98023224E-8, 3.612 digits), 12.5% slower than IEEE-754 16-bit
	"\n	#define fpxx_copy ushort" // switchable data type for direct copying (custom 16-bit floating-point format: 1-4-11, exp-15, +-1.99951168, +-6.10351562E-5, +-2.98023224E-8, 3.612 digits), 12.5% slower than IEEE-754 16-bit
	"\n	#define load(p,o) half_to_float_custom(p[o])" // special function for loading half
	"\n	#define store(p,o,x) p[o]=float_to_half_custom(x)" // special function for storing half
#else // FP32
	"\n	#define fpxx float" // switchable data type (regular 32-bit float)
	"\n	#define fpxx_copy float" // switchable data type for direct copying (regular 32-bit float)
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

#ifdef PARTICLES
	"\n	#define PARTICLES"
	"\n	#define def_particles_N "+to_string(particles_N)+"ul"
	"\n	#define def_particles_rho "+to_string(particles_rho)+"f"
#endif // PARTICLES
;}

#ifdef GRAPHICS
void LBM_Domain::Graphics::allocate(Device& device) {
	bitmap = Memory<int>(device, camera.width*camera.height);
	zbuffer = Memory<int>(device, camera.width*camera.height, 1u, lbm->get_D()>1u); // if there are multiple domains, allocate zbuffer also on host side
	camera_parameters = Memory<float>(device, 15u);
	kernel_clear = Kernel(device, bitmap.length(), "graphics_clear", bitmap, zbuffer);

	kernel_graphics_flags = Kernel(device, lbm->get_N(), "graphics_flags", lbm->flags, camera_parameters, bitmap, zbuffer);
	kernel_graphics_flags_mc = Kernel(device, lbm->get_N(), "graphics_flags_mc", lbm->flags, camera_parameters, bitmap, zbuffer);
	kernel_graphics_field = Kernel(device, lbm->get_N(), "graphics_field", lbm->flags, lbm->u, camera_parameters, bitmap, zbuffer, 0, 0, 0, 0);
#ifndef D2Q9
	kernel_graphics_streamline = Kernel(device, (lbm->get_Nx()/GRAPHICS_STREAMLINE_SPARSE)*(lbm->get_Ny()/GRAPHICS_STREAMLINE_SPARSE)*(lbm->get_Nz()/GRAPHICS_STREAMLINE_SPARSE), "graphics_streamline", lbm->flags, lbm->u, camera_parameters, bitmap, zbuffer, 0, 0, 0, 0); // 3D
#else // D2Q9
	kernel_graphics_streamline = Kernel(device, (lbm->get_Nx()/GRAPHICS_STREAMLINE_SPARSE)*(lbm->get_Ny()/GRAPHICS_STREAMLINE_SPARSE), "graphics_streamline", lbm->flags, lbm->u, camera_parameters, bitmap, zbuffer, 0, 0, 0, 0); // 2D
#endif // D2Q9
	kernel_graphics_q = Kernel(device, lbm->get_N(), "graphics_q", lbm->flags, lbm->u, camera_parameters, bitmap, zbuffer);

#ifdef FORCE_FIELD
	kernel_graphics_flags.add_parameters(lbm->F);
	kernel_graphics_flags_mc.add_parameters(lbm->F);
#endif // FORCE_FIELD

#ifdef SURFACE
	skybox = Memory<int>(device, skybox_image->width()*skybox_image->height(), 1u, skybox_image->data());
	kernel_graphics_rasterize_phi = Kernel(device, lbm->get_N(), "graphics_rasterize_phi", lbm->phi, camera_parameters, bitmap, zbuffer);
	kernel_graphics_raytrace_phi = Kernel(device, bitmap.length(), "graphics_raytrace_phi", lbm->phi, lbm->flags, skybox, camera_parameters, bitmap);
#endif // SURFACE

#ifdef TEMPERATURE
	kernel_graphics_streamline.add_parameters(lbm->T);
#endif // TEMPERATURE

#ifdef PARTICLES
	kernel_graphics_particles = Kernel(device, lbm->particles.length(), "graphics_particles", lbm->particles, camera_parameters, bitmap, zbuffer);
#endif // PARTICLES
}

bool LBM_Domain::Graphics::update_camera() {
	camera.update_matrix();
	bool change = false;
	for(uint i=0u; i<15u; i++) {
		const float data = camera.data(i);
		change |= (camera_parameters[i]!=data);
		camera_parameters[i] = data;
	}
	return change; // return false if camera parameters remain unchanged
}
void LBM_Domain::Graphics::enqueue_draw_frame(const int visualization_modes, const int slice_mode, const int slice_x, const int slice_y, const int slice_z) {
	const bool camera_update = update_camera();
#if defined(INTERACTIVE_GRAPHICS)||defined(INTERACTIVE_GRAPHICS_ASCII)
	if(!camera_update&&!camera.key_update&&lbm->get_t()==t_last_frame) return; // don't render a new frame if the scene hasn't changed since last frame
#endif // INTERACTIVE_GRAPHICS||INTERACTIVE_GRAPHICS_ASCII
	t_last_frame = lbm->get_t();
	camera.key_update = false;
	if(camera_update) camera_parameters.enqueue_write_to_device(); // camera_parameters PCIe transfer and kernel_clear execution can happen simulataneously
	kernel_clear.enqueue_run();
#ifdef SURFACE
	if((visualization_modes&0b01000000)&&lbm->get_D()==1u) kernel_graphics_raytrace_phi.enqueue_run(); // disable raytracing for multi-GPU (domain decomposition rendering doesn't work for raytracing)
	if(visualization_modes&0b00100000) kernel_graphics_rasterize_phi.enqueue_run();
#endif // SURFACE
	if((visualization_modes&0b11)==1||(visualization_modes&0b11)==2) kernel_graphics_flags.enqueue_run();
	if((visualization_modes&0b11)==2||(visualization_modes&0b11)==3) kernel_graphics_flags_mc.enqueue_run();
	if(visualization_modes&0b00000100) kernel_graphics_field.set_parameters(5u, slice_mode, slice_x-lbm->Ox, slice_y-lbm->Oy, slice_z-lbm->Oz).enqueue_run();
	if(visualization_modes&0b00001000) kernel_graphics_streamline.set_parameters(5u, slice_mode, slice_x-lbm->Ox, slice_y-lbm->Oy, slice_z-lbm->Oz).enqueue_run();
	if(visualization_modes&0b00010000) kernel_graphics_q.enqueue_run();
#ifdef PARTICLES
	if(visualization_modes&0b10000000) kernel_graphics_particles.enqueue_run();
#endif // PARTICLES
	bitmap.enqueue_read_from_device();
	if(lbm->get_D()>1u) zbuffer.enqueue_read_from_device();
}
int* LBM_Domain::Graphics::get_bitmap() { // returns pointer to zbuffer
	return bitmap.data();
}
int* LBM_Domain::Graphics::get_zbuffer() { // returns pointer to zbuffer
	return zbuffer.data();
}

string LBM_Domain::Graphics::device_defines() const { return
	"\n	#define GRAPHICS"
	"\n	#define def_background_color " +to_string(GRAPHICS_BACKGROUND_COLOR)+""
	"\n	#define def_screen_width "     +to_string(camera.width)+"u"
	"\n	#define def_screen_height "    +to_string(camera.height)+"u"
	"\n	#define def_scale_u "          +to_string(1.0f/(0.57735027f*(GRAPHICS_U_MAX)))+"f"
	"\n	#define def_scale_Q_min "      +to_string(GRAPHICS_Q_CRITERION)+"f"
	"\n	#define def_scale_F "          +to_string(1.0f/(GRAPHICS_F_MAX))+"f"
	"\n	#define def_streamline_sparse "+to_string(GRAPHICS_STREAMLINE_SPARSE)+"u"
	"\n	#define def_streamline_length "+to_string(GRAPHICS_STREAMLINE_LENGTH)+"u"
	"\n	#define def_n "                +to_string(1.333f)+"f" // refractive index of water for raytracing graphics
	"\n	#define def_attenuation "      +to_string(ln(GRAPHICS_RAYTRACING_TRANSMITTANCE)/(float)max(max(lbm->get_Nx(), lbm->get_Ny()), lbm->get_Nz()))+"f" // (negative) attenuation parameter for raytracing graphics
	"\n	#define def_absorption_color " +to_string(GRAPHICS_RAYTRACING_COLOR)+"" // absorption color of fluid for raytracing graphics

	"\n	#define COLOR_S (127<<16|127<<8|127)" // (stationary or moving) solid boundary
	"\n	#define COLOR_E (  0<<16|255<<8|  0)" // equilibrium boundary (inflow/outflow)
	"\n	#define COLOR_M (255<<16|  0<<8|255)" // cells next to moving solid boundary
	"\n	#define COLOR_T (255<<16|  0<<8|  0)" // temperature boundary
	"\n	#define COLOR_F (  0<<16|  0<<8|255)" // fluid
	"\n	#define COLOR_I (  0<<16|255<<8|255)" // interface
	"\n	#define COLOR_0 (127<<16|127<<8|127)" // regular cell or gas
	"\n	#define COLOR_X (255<<16|127<<8|  0)" // reserved type X
	"\n	#define COLOR_Y (255<<16|255<<8|  0)" // reserved type Y
	"\n	#define COLOR_P (255<<16|255<<8|191)" // particles

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
#endif // GRAPHICS



//{ thread* threads=new thread[N]; for(uint n=0u; n<N; n++) threads[n]=thread([=]() { ... }); for(uint n=0u; n<N; n++) threads[n].join(); delete[] threads; }
//#include <ppl.h> // concurrency::parallel_for(0, N, [&](int n) { ... });
//#include <omp.h> // #pragma omp parallel for \n for(int n=0; n<N; i++) { ... } // #pragma warning(disable:6993)

vector<Device_Info> smart_device_selection(const uint D) {
	const vector<Device_Info>& devices = get_devices(); // a vector of all available OpenCL devices
	vector<Device_Info> device_infos(D);
	const int user_specified_devices = (int)main_arguments.size();
	if(user_specified_devices>0) { // user has selevted specific devices as command line arguments
		if(user_specified_devices==D) { // as much specified devices as domains
			for(uint d=0; d<D; d++) device_infos[d] = select_device_with_id(to_uint(main_arguments[d]), devices); // use list of devices IDs specified by user
		} else {
			print_warning("Incorrect number of devices specified. Using single fastest device for all domains.");
			for(uint d=0; d<D; d++) device_infos[d] = select_device_with_most_flops(devices);
		}
	} else { // device auto-selection
		vector<vector<Device_Info>> device_type_ids; // a vector of all different devices, containing vectors of their device IDs
		for(uint i=0u; i<(uint)devices.size(); i++) {
			const string name_i = devices[i].name;
			bool already_exists = false;
			for(uint j=0u; j<(uint)device_type_ids.size(); j++) {
				const string name_j = device_type_ids[j][0].name;
				if(name_i==name_j) {
					device_type_ids[j].push_back(devices[i]);
					already_exists = true;
				}
			}
			if(!already_exists) device_type_ids.push_back(vector<Device_Info>(1, devices[i]));
		}
		float best_value = 0.0f;
		int best_j = -1;
		for(uint j=0u; j<(uint)device_type_ids.size(); j++) {
			const float value = device_type_ids[j][0].tflops;
			if((uint)device_type_ids[j].size()>=D && value>best_value) {
				best_value = value;
				best_j = j;
			}
		}
		if(best_j>=0) { // select all devices of fastest device type with at least D devices of the same type
			for(uint d=0; d<D; d++) device_infos[d] = device_type_ids[best_j][d];
		} else {
			print_warning("Not enough devices of the same type available. Using single fastest device for all domains.");
			for(uint d=0; d<D; d++) device_infos[d] = select_device_with_most_flops(devices);
		}
		//for(uint j=0u; j<(uint)device_type_ids.size(); j++) print_info("Device Type "+to_string(j)+" ("+device_type_ids[j][0].name+"): "+to_string((uint)device_type_ids[j].size())+"x");
	}
	return device_infos;
}

LBM::LBM(const uint Nx, const uint Ny, const uint Nz, const float nu, const float fx, const float fy, const float fz, const float sigma, const float alpha, const float beta, const uint particles_N, const float particles_rho) // single device
	:LBM(Nx, Ny, Nz, 1u, 1u, 1u, nu, fx, fy, fz, sigma, alpha, beta, particles_N, particles_rho) { // delegating constructor
}
LBM::LBM(const uint Nx, const uint Ny, const uint Nz, const float nu, const uint particles_N, const float particles_rho)
	:LBM(Nx, Ny, Nz, 1u, 1u, 1u, nu, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, particles_N, particles_rho) { // delegating constructor
}
LBM::LBM(const uint Nx, const uint Ny, const uint Nz, const float nu, const float fx, const float fy, const float fz, const uint particles_N, const float particles_rho)
	:LBM(Nx, Ny, Nz, 1u, 1u, 1u, nu, fx, fy, fz, 0.0f, 0.0f, 0.0f, particles_N, particles_rho) { // delegating constructor
}
LBM::LBM(const uint Nx, const uint Ny, const uint Nz, const uint Dx, const uint Dy, const uint Dz, const float nu, const float fx, const float fy, const float fz, const float sigma, const float alpha, const float beta, const uint particles_N, const float particles_rho) { // multiple devices
	this->Nx = Nx; this->Ny = Ny; this->Nz = Nz;
	this->Dx = Dx; this->Dy = Dy; this->Dz = Dz;
	const uint D = Dx*Dy*Dz;
	const uint Hx=Dx>1u, Hy=Dy>1u, Hz=Dz>1u; // halo offsets
	const vector<Device_Info>& device_infos = smart_device_selection(D);
	sanity_checks_constructor(device_infos, Nx, Ny, Nz, Dx, Dy, Dz, nu, fx, fy, fz, sigma, alpha, beta, particles_N, particles_rho);
	lbm = new LBM_Domain*[D];
	for(uint d=0u; d<D; d++) { // { thread* threads=new thread[D]; for(uint d=0u; d<D; d++) threads[d]=thread([=]() {
		const uint x=((uint)d%(Dx*Dy))%Dx, y=((uint)d%(Dx*Dy))/Dx, z=(uint)d/(Dx*Dy); // d = x+(y+z*Dy)*Dx
		lbm[d] = new LBM_Domain(device_infos[d], Nx/Dx+2u*Hx, Ny/Dy+2u*Hy, Nz/Dz+2u*Hz, Dx, Dy, Dz, (int)(x*Nx/Dx)-(int)Hx, (int)(y*Ny/Dy)-(int)Hy, (int)(z*Nz/Dz)-(int)Hz, nu, fx, fy, fz, sigma, alpha, beta, particles_N, particles_rho);
	} // }); for(uint d=0u; d<D; d++) threads[d].join(); delete[] threads; }
	{
		Memory<float>** buffers_rho = new Memory<float>*[D];
		for(uint d=0u; d<D; d++) buffers_rho[d] = &(lbm[d]->rho);
		rho = Memory_Container(this, buffers_rho, "rho");
	} {
		Memory<float>** buffers_u = new Memory<float>*[D];
		for(uint d=0u; d<D; d++) buffers_u[d] = &(lbm[d]->u);
		u = Memory_Container(this, buffers_u, "u");
	} {
		Memory<uchar>** buffers_flags = new Memory<uchar>*[D];
		for(uint d=0u; d<D; d++) buffers_flags[d] = &(lbm[d]->flags);
		flags = Memory_Container(this, buffers_flags, "flags");
	} {
#ifdef FORCE_FIELD
		Memory<float>** buffers_F = new Memory<float>*[D];
		for(uint d=0u; d<D; d++) buffers_F[d] = &(lbm[d]->F);
		F = Memory_Container(this, buffers_F, "F");
#endif // FORCE_FIELD
	} {
#ifdef SURFACE
		Memory<float>** buffers_phi = new Memory<float>*[D];
		for(uint d=0u; d<D; d++) buffers_phi[d] = &(lbm[d]->phi);
		phi = Memory_Container(this, buffers_phi, "phi");
#endif // SURFACE
	} {
#ifdef TEMPERATURE
		Memory<float>** buffers_T = new Memory<float>*[D];
		for(uint d=0u; d<D; d++) buffers_T[d] = &(lbm[d]->T);
		T = Memory_Container(this, buffers_T, "T");
#endif // TEMPERATURE
	} {
#ifdef PARTICLES
		particles = &(lbm[0]->particles);
#endif // PARTICLES
	}
#ifdef GRAPHICS
	graphics = Graphics(this);
#endif // GRAPHICS
	info.initialize(this);
}
LBM::~LBM() {
	info.print_finalize();
	for(uint d=0u; d<get_D(); d++) delete lbm[d];
	delete[] lbm;
}

void LBM::sanity_checks_constructor(const vector<Device_Info>& device_infos, const uint Nx, const uint Ny, const uint Nz, const uint Dx, const uint Dy, const uint Dz, const float nu, const float fx, const float fy, const float fz, const float sigma, const float alpha, const float beta, const uint particles_N, const float particles_rho) { // sanity checks on grid resolution and extension support
	if((ulong)Nx*(ulong)Ny*(ulong)Nz==0ull) print_error("Grid point number is 0: "+to_string(Nx)+"x"+to_string(Ny)+"x"+to_string(Nz)+" = 0.");
	if(Nx%Dx!=0u || Ny%Dy!=0u || Nz%Dz!=0u) print_error("LBM grid ("+to_string(Nx)+"x"+to_string(Ny)+"x"+to_string(Nz)+") is not equally divisible in domains ("+to_string(Dx)+"x"+to_string(Dy)+"x"+to_string(Dz)+").");
	if(Dx*Dy*Dz==0u) print_error("You specified 0 LBM grid domains ("+to_string(Dx)+"x"+to_string(Dy)+"x"+to_string(Dz)+"). There has to be at least 1 domain in every direction. Check your input in LBM constructor.");
	const uint local_Nx=Nx/Dx+2u*(Dx>1u), local_Ny=Ny/Dy+2u*(Dy>1u), local_Nz=Nz/Dz+2u*(Dz>1u);
	if((ulong)local_Nx*(ulong)local_Ny*(ulong)local_Nz>=(ulong)max_uint) print_error("Single domain grid resolution is too large: "+to_string(local_Nx)+"x"+to_string(local_Ny)+"x"+to_string(local_Nz)+" > 2^32.");
	uint memory_available = max_uint; // in MB
	for(Device_Info device_info : device_infos) memory_available = min(memory_available, device_info.memory);
	uint memory_required = (uint)((ulong)Nx*(ulong)Ny*(ulong)Nz/((ulong)(Dx*Dy*Dz))*((ulong)velocity_set*sizeof(fpxx)+17ull)/1048576ull); // in MB
	if(memory_required>memory_available) {
		float factor = cbrt((float)memory_available/(float)memory_required);
		const uint maxNx=(uint)(factor*(float)Nx), maxNy=(uint)(factor*(float)Ny), maxNz=(uint)(factor*(float)Nz);
		string message = "Grid resolution ("+to_string(Nx)+", "+to_string(Ny)+", "+to_string(Nz)+") is too large: "+to_string(Dx*Dy*Dz)+"x "+to_string(memory_required)+" MB required, "+to_string(Dx*Dy*Dz)+"x "+to_string(memory_available)+" MB available. Largest possible resolution is ("+to_string(maxNx)+", "+to_string(maxNy)+", "+to_string(maxNz)+"). Restart the simulation with lower resolution or on different device(s) with more memory.";
#if !defined(FP16S)&&!defined(FP16C)
		uint memory_required_fp16 = (uint)((ulong)Nx*(ulong)Ny*(ulong)Nz/((ulong)(Dx*Dy*Dz))*((ulong)velocity_set*2ull+17ull)/1048576ull); // in MB
		float factor_fp16 = cbrt((float)memory_available/(float)memory_required_fp16);
		const uint maxNx_fp16=(uint)(factor_fp16*(float)Nx), maxNy_fp16=(uint)(factor_fp16*(float)Ny), maxNz_fp16=(uint)(factor_fp16*(float)Nz);
		message += " Consider using FP16S/FP16C memory compression to double maximum grid resolution to a maximum of ("+to_string(maxNx_fp16)+", "+to_string(maxNy_fp16)+", "+to_string(maxNz_fp16)+"); for this, uncomment \"#define FP16S\" or \"#define FP16C\" in defines.hpp.";
#endif // !FP16S&&!FP16C
		print_error(message);
	}
	if(nu==0.0f) print_error("Viscosity cannot be 0. Change it in setup.cpp."); // sanity checks for viscosity
	else if(nu<0.0f) print_error("Viscosity cannot be negative. Remove the \"-\" in setup.cpp.");
#ifdef D2Q9
	if(Nz!=1u) print_error("D2Q9 is the 2D velocity set. You have to set Nz=1u in the LBM constructor! Currently you have set Nz="+to_string(Nz)+"u.");
#endif // D2Q9
#if !defined(SRT)&&!defined(TRT)
	print_error("No LBM collision operator selected. Uncomment either \"#define SRT\" or \"#define TRT\" in defines.hpp");
#elif defined(SRT)&&defined(TRT)
	print_error("Too many LBM collision operators selected. Comment out either \"#define SRT\" or \"#define TRT\" in defines.hpp");
#endif // SRT && TRT
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
#ifdef PARTICLES
	if(particles_N==0u) print_error("The PARTICLES extension is enabled but the number of particles is set to 0. Comment out \"#define PARTICLES\" in defines.hpp.");
	if(get_D()>1u) print_error("The PARTICLES extension is not supported in multi-GPU mode.");
#if !defined(VOLUME_FORCE)||!defined(FORCE_FIELD)
	if(particles_rho!=1.0f) print_error("Particle density is set unequal to 1, but particle-fluid 2-way-coupling is not enabled. Uncomment both \"#define VOLUME_FORCE\" and \"#define FORCE_FIELD\" in defines.hpp.");
#endif // !VOLUME_FORCE||!FORCE_FIELD
#ifdef FORCE_FIELD
	if(particles_rho==1.0f) print_warning("Particle density is set to 1, so particles behave as passive tracers without acting a force on the fluid, but particle-fluid 2-way-coupling is enabled. You may comment out \"#define FORCE_FIELD\" in defines.hpp.");
#endif // FORCE_FIELD
#else // PARTICLES
	if(particles_N>0u) print_error("The PARTICLES extension is disabled but the number of particles is set to "+to_string(particles_N)+">0. Uncomment \"#define PARTICLES\" in defines.hpp.");
#endif // PARTICLES
}

void LBM::sanity_checks_initialization() { // sanity checks during initialization on used extensions based on used flags
	bool moving_boundaries_used=false, equilibrium_boundaries_used=false, surface_used=false, temperature_used=false; // identify used extensions based used flags
	for(ulong n=0ull; n<get_N(); n++) {
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

void LBM::initialize() { // write all data fields to device and call kernel_initialize
	sanity_checks_initialization();

	for(uint d=0u; d<get_D(); d++) lbm[d]->rho.enqueue_write_to_device();
	for(uint d=0u; d<get_D(); d++) lbm[d]->u.enqueue_write_to_device();
	for(uint d=0u; d<get_D(); d++) lbm[d]->flags.enqueue_write_to_device();
#ifdef FORCE_FIELD
	for(uint d=0u; d<get_D(); d++) lbm[d]->F.enqueue_write_to_device();
#endif // FORCE_FIELD
#ifdef SURFACE
	for(uint d=0u; d<get_D(); d++) lbm[d]->phi.enqueue_write_to_device();
#endif // SURFACE
#ifdef TEMPERATURE
	for(uint d=0u; d<get_D(); d++) lbm[d]->T.enqueue_write_to_device();
#endif // TEMPERATURE
#ifdef PARTICLES
	for(uint d=0u; d<get_D(); d++) lbm[d]->particles.enqueue_write_to_device();
#endif // PARTICLES

	for(uint d=0u; d<get_D(); d++) lbm[d]->increment_time_step(); // the communicate calls at initialization need an odd time step
	communicate_rho_u_flags();
#ifdef SURFACE
	communicate_phi_massex_flags();
#endif // SURFACE
	for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_initialize(); // odd time step is baked-in the kernel
	communicate_rho_u_flags();
#ifdef SURFACE
	communicate_phi_massex_flags();
#endif // SURFACE
	communicate_fi(); // time step must be odd here
#ifdef TEMPERATURE
	communicate_gi(); // time step must be odd here
#endif // TEMPERATURE
	for(uint d=0u; d<get_D(); d++) lbm[d]->finish_queue();
	for(uint d=0u; d<get_D(); d++) lbm[d]->reset_time_step(); // set time step to 0 again
	initialized = true;
}

void LBM::do_time_step() { // call kernel_stream_collide to perform one LBM time step
#ifdef SURFACE
	for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_surface_0();
#endif // SURFACE
	for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_stream_collide(); // run LBM stream_collide kernel after domain communication
#ifdef SURFACE
	communicate_rho_u_flags();
	for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_surface_1();
	communicate_flags();
	for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_surface_2();
	communicate_flags();
	for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_surface_3();
	communicate_phi_massex_flags();
#endif // SURFACE
	communicate_fi();
#ifdef TEMPERATURE
	communicate_gi();
#endif // TEMPERATURE
#ifdef PARTICLES
	for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_integrate_particles(); // intgegrate particles forward in time and couple particles to fluid
#endif // PARTICLES
	if(get_D()==1u) for(uint d=0u; d<get_D(); d++) lbm[d]->finish_queue(); // this additional domain synchronization barrier is only required in single-GPU, as communication calls already provide all necessary synchronization barriers in multi-GPU
	for(uint d=0u; d<get_D(); d++) lbm[d]->increment_time_step();
}

void LBM::run(const ulong steps) { // initializes the LBM simulation (copies data to device and runs initialize kernel), then runs LBM
	info.append(steps, get_t());
	if(!initialized) {
		initialize();
		info.print_initialize(); // only print setup info if the setup is new (run() was not called before)
	}
	Clock clock;
	for(ulong i=1ull; i<=steps; i++) {
#if defined(INTERACTIVE_GRAPHICS)||defined(INTERACTIVE_GRAPHICS_ASCII)
		while(!key_P&&running) sleep(0.016);
		if(!running) break;
#endif // INTERACTIVE_GRAPHICS_ASCII || INTERACTIVE_GRAPHICS
		clock.start();
		do_time_step();
		info.update(clock.stop());
	}
	if(get_D()>1u) for(uint d=0u; d<get_D(); d++) lbm[d]->finish_queue(); // wait for everything to finish (multi-GPU only)
}

void LBM::update_fields() { // update fields (rho, u, T) manually
	for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_update_fields();
	for(uint d=0u; d<get_D(); d++) lbm[d]->finish_queue();
}

void LBM::reset() { // reset simulation (takes effect in following run() call)
	initialized = false;
}

#ifdef FORCE_FIELD
void LBM::calculate_force_on_boundaries() { // calculate forces from fluid on TYPE_S nodes
	for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_calculate_force_on_boundaries();
	for(uint d=0u; d<get_D(); d++) lbm[d]->finish_queue();
}
float3 LBM::calculate_force_on_object(const uchar flag_marker) { // add up force for all nodes flagged with flag_marker
	double3 force(0.0, 0.0, 0.0);
	for(ulong n=0ull; n<get_N(); n++) {
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
	ulong counter = 0ull;
	for(ulong n=0ull; n<get_N(); n++) {
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
	for(ulong n=0ull; n<get_N(); n++) {
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
	for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_update_moving_boundaries();
	communicate_rho_u_flags();
	for(uint d=0u; d<get_D(); d++) lbm[d]->finish_queue();
}
#endif // MOVING_BOUNDARIES

#if defined(PARTICLES)&&!defined(FORCE_FIELD)
void LBM::integrate_particles(const ulong steps, const uint time_step_multiplicator) { // intgegrate passive tracer particles forward in time in stationary flow field
	info.append(steps, get_t());
	Clock clock;
	for(ulong i=1ull; i<=steps; i+=(ulong)time_step_multiplicator) {
#if defined(INTERACTIVE_GRAPHICS)||defined(INTERACTIVE_GRAPHICS_ASCII)
		while(!key_P&&running) sleep(0.016);
		if(!running) break;
#endif // INTERACTIVE_GRAPHICS_ASCII || INTERACTIVE_GRAPHICS
		clock.start();
		for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_integrate_particles(time_step_multiplicator);
		for(uint d=0u; d<get_D(); d++) lbm[d]->finish_queue();
		for(uint d=0u; d<get_D(); d++) lbm[d]->increment_time_step(time_step_multiplicator);
		info.update(clock.stop());
	}
}
#endif // PARTICLES&&!FORCE_FIELD

void LBM::write_status(const string& path) { // write LBM status report to a .txt file
	string status = "";
	status += "Grid Resolution = ("+to_string(Nx)+", "+to_string(Ny)+", "+to_string(Nz)+")\n";
	status += "LBM type = D"+string(get_velocity_set()==9 ? "2" : "3")+"Q"+to_string(get_velocity_set())+" "+info.collision+"\n";
	status += "Memory Usage = "+to_string(info.cpu_mem_required)+" MB (CPU), "+to_string(info.gpu_mem_required)+" MB (GPU)\n";
	status += "Maximum Allocation Size = "+to_string((uint)(get_N()*(ulong)(get_velocity_set()*sizeof(fpxx))/1048576ull))+" MB\n";
	status += "Time Step = "+to_string(get_t())+" / "+(info.steps==max_ulong ? "infinite" : to_string(info.steps))+"\n";
	status += "Kinematic Viscosity = "+to_string(get_nu())+"\n";
	status += "Relaxation Time = "+to_string(get_tau())+"\n";
	status += "Maximum Reynolds Number = "+to_string(get_Re_max())+"\n";
#ifdef VOLUME_FORCE
	status += "Volume Force = ("+to_string(get_fx())+", "+to_string(get_fy())+", "+to_string(get_fz())+")\n";
#endif // VOLUME_FORCE
#ifdef SURFACE
	status += "Surface Tension Coefficient = "+to_string(get_sigma())+"\n";
#endif // SURFACE
#ifdef TEMPERATURE
	status += "Thermal Diffusion Coefficient = "+to_string(get_alpha())+"\n";
	status += "Thermal Expansion Coefficient = "+to_string(get_beta())+"\n";
#endif // TEMPERATURE
	const string filename = default_filename(path, "status", ".txt", get_t());
	write_file(filename, status);
}

void LBM::voxelize_mesh_on_device(const Mesh* mesh, const uchar flag, const float3& rotation_center, const float3& linear_velocity, const float3& rotational_velocity) { // voxelize triangle mesh
	if(get_D()==1u) {
		lbm[0]->voxelize_mesh_on_device(mesh, flag, rotation_center, linear_velocity, rotational_velocity); // if this crashes on Windows, create a TdrDelay 32-bit DWORD with decimal value 300 in Computer\HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\GraphicsDrivers
	} else {
		thread* threads=new thread[get_D()]; for(uint d=0u; d<get_D(); d++) threads[d]=thread([=]() {
			lbm[d]->voxelize_mesh_on_device(mesh, flag, rotation_center, linear_velocity, rotational_velocity);
		}); for(uint d=0u; d<get_D(); d++) threads[d].join(); delete[] threads;
	}
#ifdef MOVING_BOUNDARIES
	if(flag==TYPE_S&&(length(linear_velocity)>0.0f||length(rotational_velocity)>0.0f)) update_moving_boundaries();
#endif // MOVING_BOUNDARIES
	if(!initialized) {
		flags.read_from_device();
		u.read_from_device();
	}
}
void LBM::unvoxelize_mesh_on_device(const Mesh* mesh, const uchar flag) { // remove voxelized triangle mesh from LBM grid by removing all flags in mesh bounding box (only required when bounding box size changes during re-voxelization)
	for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_unvoxelize_mesh_on_device(mesh, flag);
	for(uint d=0u; d<get_D(); d++) lbm[d]->finish_queue();
}
void LBM::write_mesh_to_vtk(const Mesh* mesh, const string& path) { // write mesh to binary .vtk file
	const string header_1 = "# vtk DataFile Version 3.0\nData\nBINARY\nDATASET POLYDATA\nPOINTS "+to_string(3u*mesh->triangle_number)+" float\n";
	const string header_2 = "POLYGONS "+to_string(mesh->triangle_number)+" "+to_string(4u*mesh->triangle_number)+"\n";
	float* points = new float[9u*mesh->triangle_number];
	int* triangles = new int[4u*mesh->triangle_number];
	for(uint i=0u; i<mesh->triangle_number; i++) {
		points[9u*i   ] = reverse_bytes(mesh->p0[i].x-center().x);
		points[9u*i+1u] = reverse_bytes(mesh->p0[i].y-center().y);
		points[9u*i+2u] = reverse_bytes(mesh->p0[i].z-center().z);
		points[9u*i+3u] = reverse_bytes(mesh->p1[i].x-center().x);
		points[9u*i+4u] = reverse_bytes(mesh->p1[i].y-center().y);
		points[9u*i+5u] = reverse_bytes(mesh->p1[i].z-center().z);
		points[9u*i+6u] = reverse_bytes(mesh->p2[i].x-center().x);
		points[9u*i+7u] = reverse_bytes(mesh->p2[i].y-center().y);
		points[9u*i+8u] = reverse_bytes(mesh->p2[i].z-center().z);
		triangles[4u*i   ] = reverse_bytes(3); // 3 vertices per triangle
		triangles[4u*i+1u] = reverse_bytes(3*(int)i  ); // vertex 0
		triangles[4u*i+2u] = reverse_bytes(3*(int)i+1); // vertex 1
		triangles[4u*i+3u] = reverse_bytes(3*(int)i+2); // vertex 2
	}
	const string filename = default_filename(path, "mesh", ".vtk", get_t());
	create_folder(filename);
	std::ofstream file(filename, std::ios::out|std::ios::binary);
	file.write(header_1.c_str(), header_1.length()); // write non-binary file header
	file.write((char*)points, 4u*9u*mesh->triangle_number); // write binary data
	file.write(header_2.c_str(), header_2.length()); // write non-binary file header
	file.write((char*)triangles, 4u*4u*mesh->triangle_number); // write binary data
	file.close();
	delete[] points;
	delete[] triangles;
	info.allow_rendering = false; // temporarily disable interactive rendering
	print_info("File \""+filename+"\" saved.");
	info.allow_rendering = true;
}
void LBM::voxelize_stl(const string& path, const float3& center, const float3x3& rotation, const float size, const uchar flag) { // voxelize triangle mesh
	const Mesh* mesh = read_stl(path, this->size(), center, rotation, size);
	flags.write_to_device();
	voxelize_mesh_on_device(mesh, flag);
	delete mesh;
	flags.read_from_device();
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

#ifdef GRAPHICS
int* LBM::Graphics::draw_frame() {
#ifndef UPDATE_FIELDS
	if(visualization_modes&0b00011100) {
		for(uint d=0u; d<lbm->get_D(); d++) lbm->lbm[d]->enqueue_update_fields(); // only call update_fields() if the time step has changed since the last rendered frame
		//for(uint d=0u; d<lbm->get_D(); d++) lbm->communicate_rho_u_flags();
	}
#endif // UPDATE_FIELDS

	if(key_1) { visualization_modes = (visualization_modes&~0b11)|(((visualization_modes&0b11)+1)%4); key_1 = false; }
	if(key_2) { visualization_modes ^= 0b00000100; key_2 = false; }
	if(key_3) { visualization_modes ^= 0b00001000; key_3 = false; }
	if(key_4) { visualization_modes ^= 0b00010000; key_4 = false; }
	if(key_5) { visualization_modes ^= 0b00100000; key_5 = false; }
	if(key_6) { visualization_modes ^= 0b01000000; key_6 = false; }
	if(key_7) { visualization_modes ^= 0b10000000; key_7 = false; }

	if(key_T) {
		slice_mode = (slice_mode+1)%8; key_T = false;
	}
	if(slice_mode==1u) {
		if(key_Q) { slice_x = clamp(slice_x-1, 0, (int)lbm->get_Nx()-1); key_Q = false; }
		if(key_E) { slice_x = clamp(slice_x+1, 0, (int)lbm->get_Nx()-1); key_E = false; }
	}
	if(slice_mode==2u) {
		if(key_Q) { slice_y = clamp(slice_y-1, 0, (int)lbm->get_Ny()-1); key_Q = false; }
		if(key_E) { slice_y = clamp(slice_y+1, 0, (int)lbm->get_Ny()-1); key_E = false; }
	}
	if(slice_mode==3u) {
		if(key_Q) { slice_z = clamp(slice_z-1, 0, (int)lbm->get_Nz()-1); key_Q = false; }
		if(key_E) { slice_z = clamp(slice_z+1, 0, (int)lbm->get_Nz()-1); key_E = false; }
	}

	for(uint d=0u; d<lbm->get_D(); d++) lbm->lbm[d]->graphics.enqueue_draw_frame(visualization_modes, slice_mode, slice_x, slice_y, slice_z);
	for(uint d=0u; d<lbm->get_D(); d++) lbm->lbm[d]->finish_queue();
	int* bitmap = lbm->lbm[0]->graphics.get_bitmap();
	int* zbuffer = lbm->lbm[0]->graphics.get_zbuffer();
	for(uint d=1u; d<lbm->get_D(); d++) {
		const int* bitmap_d = lbm->lbm[d]->graphics.get_bitmap(); // each domain renders its own frame
		const int* zbuffer_d = lbm->lbm[d]->graphics.get_zbuffer();
		for(uint i=0u; i<camera.width*camera.height; i++) {
			const int zdi = zbuffer_d[i];
			if(zdi>zbuffer[i]) {
				bitmap[i] = bitmap_d[i]; // overlay frames using their z-buffers
				zbuffer[i] = zdi;
			}
		}
	}
	return bitmap;
}

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
	int* image_data = draw_frame(); // make sure the frame is fully rendered
	Image* image = new Image(camera.width, camera.height, image_data);
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
	int* image_data = draw_frame(); // make sure the frame is fully rendered
	const string filename = default_filename(path, name, extension, lbm->get_t());
	const uint xa=max(min(x1, x2), 0u), xb=min(max(x1, x2), camera.width ); // sort coordinates if necessary
	const uint ya=max(min(y1, y2), 0u), yb=min(max(y1, y2), camera.height);
	Image* image = new Image(xb-xa, yb-ya); // create local copy of frame buffer
	for(uint y=0u; y<image->height(); y++) for(uint x=0u; x<image->width(); x++) image->set_color(x, y, image_data[camera.width*(ya+y)+(xa+x)]);
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



void LBM_Domain::allocate_transfer(Device& device) { // allocate all memory for multi-device trqansfer
	ulong Amax = 0ull; // maximum domain side area of communicated directions
	if(Dx>1u) Amax = max(Amax, (ulong)Ny*(ulong)Nz); // Ax
	if(Dy>1u) Amax = max(Amax, (ulong)Nz*(ulong)Nx); // Ay
	if(Dz>1u) Amax = max(Amax, (ulong)Nx*(ulong)Ny); // Az

	transfer_buffer_p = Memory<char>(device, Amax, max(transfers*(uint)sizeof(fpxx), 17u)); // only allocate one set of transfer buffers in plus/minus directions, for all x/y/z transfers
	transfer_buffer_m = Memory<char>(device, Amax, max(transfers*(uint)sizeof(fpxx), 17u));

	kernel_transfer[enum_transfer_field::fi              ][0] = Kernel(device, 0u, "transfer_extract_fi"              , 0u, t, transfer_buffer_p, transfer_buffer_m, fi);
	kernel_transfer[enum_transfer_field::fi              ][1] = Kernel(device, 0u, "transfer__insert_fi"              , 0u, t, transfer_buffer_p, transfer_buffer_m, fi);
	kernel_transfer[enum_transfer_field::rho_u_flags     ][0] = Kernel(device, 0u, "transfer_extract_rho_u_flags"     , 0u, t, transfer_buffer_p, transfer_buffer_m, rho, u, flags);
	kernel_transfer[enum_transfer_field::rho_u_flags     ][1] = Kernel(device, 0u, "transfer__insert_rho_u_flags"     , 0u, t, transfer_buffer_p, transfer_buffer_m, rho, u, flags);
#ifdef SURFACE
	kernel_transfer[enum_transfer_field::flags           ][0] = Kernel(device, 0u, "transfer_extract_flags"           , 0u, t, transfer_buffer_p, transfer_buffer_m, flags);
	kernel_transfer[enum_transfer_field::flags           ][1] = Kernel(device, 0u, "transfer__insert_flags"           , 0u, t, transfer_buffer_p, transfer_buffer_m, flags);
	kernel_transfer[enum_transfer_field::phi_massex_flags][0] = Kernel(device, 0u, "transfer_extract_phi_massex_flags", 0u, t, transfer_buffer_p, transfer_buffer_m, phi, massex, flags);
	kernel_transfer[enum_transfer_field::phi_massex_flags][1] = Kernel(device, 0u, "transfer__insert_phi_massex_flags", 0u, t, transfer_buffer_p, transfer_buffer_m, phi, massex, flags);
#endif // SURFACE
#ifdef TEMPERATURE
	kernel_transfer[enum_transfer_field::gi              ][0] = Kernel(device, 0u, "transfer_extract_gi"              , 0u, t, transfer_buffer_p, transfer_buffer_m, gi);
	kernel_transfer[enum_transfer_field::gi              ][1] = Kernel(device, 0u, "transfer__insert_gi"              , 0u, t, transfer_buffer_p, transfer_buffer_m, gi);
#endif // TEMPERATURE
}

ulong LBM_Domain::get_area(const uint direction) {
	const ulong A[3] = { (ulong)Ny*(ulong)Nz, (ulong)Nz*(ulong)Nx, (ulong)Nx*(ulong)Ny };
	return A[direction];
}
void LBM_Domain::enqueue_transfer_extract_field(Kernel& kernel_transfer_extract_field, const uint direction, const uint bytes_per_cell) {
	kernel_transfer_extract_field.set_ranges(get_area(direction)); // direction: x=0, y=1, z=2
	kernel_transfer_extract_field.set_parameters(0u, direction, get_t()).enqueue_run(); // selective in-VRAM copy
	transfer_buffer_p.enqueue_read_from_device(0ull, kernel_transfer_extract_field.range()*(ulong)bytes_per_cell); // PCIe copy (+)
	transfer_buffer_m.enqueue_read_from_device(0ull, kernel_transfer_extract_field.range()*(ulong)bytes_per_cell); // PCIe copy (-)
}
void LBM_Domain::enqueue_transfer_insert_field(Kernel& kernel_transfer_insert_field, const uint direction, const uint bytes_per_cell) {
	kernel_transfer_insert_field.set_ranges(get_area(direction)); // direction: x=0, y=1, z=2
	transfer_buffer_p.enqueue_write_to_device(0ull, kernel_transfer_insert_field.range()*(ulong)bytes_per_cell); // PCIe copy (+)
	transfer_buffer_m.enqueue_write_to_device(0ull, kernel_transfer_insert_field.range()*(ulong)bytes_per_cell); // PCIe copy (-)
	kernel_transfer_insert_field.set_parameters(0u, direction, get_t()).enqueue_run(); // selective in-VRAM copy
}
void LBM::communicate_field(const enum_transfer_field field, const uint bytes_per_cell) {
	if(Dx>1u) { // communicate in x-direction
		for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_transfer_extract_field(lbm[d]->kernel_transfer[field][0], 0u, bytes_per_cell); // selective in-VRAM copy (x) + PCIe copy
		for(uint d=0u; d<get_D(); d++) lbm[d]->finish_queue(); // domain synchronization barrier
		for(uint d=0u; d<get_D(); d++) {
			const uint x=(d%(Dx*Dy))%Dx, y=(d%(Dx*Dy))/Dx, z=d/(Dx*Dy), dxp=((x+1u)%Dx)+(y+z*Dy)*Dx; // d = x+(y+z*Dy)*Dx
			lbm[d]->transfer_buffer_p.exchange_host_buffer(lbm[dxp]->transfer_buffer_m.exchange_host_buffer(lbm[d]->transfer_buffer_p.data())); // CPU pointer swaps
		}
		for(uint d=0u; d<get_D(); d++) lbm[d]-> enqueue_transfer_insert_field(lbm[d]->kernel_transfer[field][1], 0u, bytes_per_cell); // PCIe copy + selective in-VRAM copy (x)
	}
	if(Dy>1u) { // communicate in y-direction
		for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_transfer_extract_field(lbm[d]->kernel_transfer[field][0], 1u, bytes_per_cell); // selective in-VRAM copy (y) + PCIe copy
		for(uint d=0u; d<get_D(); d++) lbm[d]->finish_queue(); // domain synchronization barrier
		for(uint d=0u; d<get_D(); d++) {
			const uint x=(d%(Dx*Dy))%Dx, y=(d%(Dx*Dy))/Dx, z=d/(Dx*Dy), dyp=x+(((y+1u)%Dy)+z*Dy)*Dx; // d = x+(y+z*Dy)*Dx
			lbm[d]->transfer_buffer_p.exchange_host_buffer(lbm[dyp]->transfer_buffer_m.exchange_host_buffer(lbm[d]->transfer_buffer_p.data())); // CPU pointer swaps
		}
		for(uint d=0u; d<get_D(); d++) lbm[d]-> enqueue_transfer_insert_field(lbm[d]->kernel_transfer[field][1], 1u, bytes_per_cell); // PCIe copy + selective in-VRAM copy (y)
	}
	if(Dz>1u) { // communicate in z-direction
		for(uint d=0u; d<get_D(); d++) lbm[d]->enqueue_transfer_extract_field(lbm[d]->kernel_transfer[field][0], 2u, bytes_per_cell); // selective in-VRAM copy (z) + PCIe copy
		for(uint d=0u; d<get_D(); d++) lbm[d]->finish_queue(); // domain synchronization barrier
		for(uint d=0u; d<get_D(); d++) {
			const uint x=(d%(Dx*Dy))%Dx, y=(d%(Dx*Dy))/Dx, z=d/(Dx*Dy), dzp=x+(y+((z+1u)%Dz)*Dy)*Dx; // d = x+(y+z*Dy)*Dx
			lbm[d]->transfer_buffer_p.exchange_host_buffer(lbm[dzp]->transfer_buffer_m.exchange_host_buffer(lbm[d]->transfer_buffer_p.data())); // CPU pointer swaps
		}
		for(uint d=0u; d<get_D(); d++) lbm[d]-> enqueue_transfer_insert_field(lbm[d]->kernel_transfer[field][1], 2u, bytes_per_cell); // PCIe copy + selective in-VRAM copy (z)
	}
}

void LBM::communicate_fi() {
	communicate_field(enum_transfer_field::fi, transfers*sizeof(fpxx));
}
void LBM::communicate_rho_u_flags() {
	communicate_field(enum_transfer_field::rho_u_flags, 17u);
}
#ifdef SURFACE
void LBM::communicate_flags() {
	communicate_field(enum_transfer_field::flags, 1u);
}
void LBM::communicate_phi_massex_flags() {
	communicate_field(enum_transfer_field::phi_massex_flags, 9u);
}
#endif // SURFACE
#ifdef TEMPERATURE
void LBM::communicate_gi() {
	communicate_field(enum_transfer_field::gi, sizeof(fpxx));
}
#endif // TEMPERATURE