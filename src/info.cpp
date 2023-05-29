#include "info.hpp"
#include "lbm.hpp"

Info info;

void Info::initialize(LBM* lbm) {
	this->lbm = lbm;
	device_transfer = lbm->get_velocity_set()*(2u*sizeof(fpxx))+17u; // lattice.set()*(2*fi) + flags + rho + 3*u
#ifndef UPDATE_FIELDS
	device_transfer -= 16u; // rho, u
#endif // UPDATE_FIELDS
#if defined(SRT)
	collision = "SRT";
#elif defined(TRT)
	collision = "TRT";
#endif // TRT
#if defined(FP16S)
	collision += " (FP32/FP16S)";
#elif defined(FP16C)
	collision += " (FP32/FP16C)";
#else // FP32
	collision += " (FP32/FP32)";
#endif // FP32
#if defined(MOVING_BOUNDARIES)||defined(SURFACE)||defined(TEMPERATURE)
	device_transfer += (lbm->get_velocity_set()-1u)*1u; // neighbor flags have to be loaded
#endif // MOVING_BOUNDARIES, SURFACE or TEMPERATURE
#ifdef SURFACE
	host_allocation += 4u; // phi
	device_transfer += (1u+(2u*lbm->get_velocity_set()-1u)*sizeof(fpxx)+8u+(lbm->get_velocity_set()-1u)*4u) + 1u + 1u + (4u+lbm->get_velocity_set()+4u+4u+4u); // surface_0 (flags, fi, mass, massex), surface_1 (flags), surface_2 (flags), surface_3 (rho, flags, mass, massex, phi)
#endif // SURFACE
#ifdef TEMPERATURE
	host_allocation += 4u; // T
	device_transfer += 7u*2u*sizeof(fpxx)+4u; // 2*gi, T
#endif // TEMPERATURE
	cpu_mem_required = (uint)(lbm->get_N()*(ulong)host_allocation/1048576ull); // reset to get valid values for consecutive simulations
	gpu_mem_required = lbm->lbm[0]->get_device().info.memory_used;
}
void Info::append(const ulong steps, const ulong t) {
	this->steps = steps; // has to be executed before info.print_initialize()
	this->steps_last = t; // reset last step count if multiple run() commands are executed consecutively
	this->runtime_last = runtime; // reset last runtime if multiple run() commands are executed consecutively
}
void Info::update(const double dt) {
	this->dt = dt; // exact dt
	this->dt_smooth = (dt+0.3)/(0.3/dt_smooth+1.0); // smoothed dt
	this->runtime += dt; // skip first step since it is likely slower than average
}
double Info::time() const { // returns either elapsed time or remaining time
	return steps==max_ulong ? runtime : ((double)steps/(double)(lbm->get_t()-steps_last)-1.0)*(runtime-runtime_last); // time estimation on average so far
	//return steps==max_ulong ? runtime : ((double)steps-(double)(lbm->get_t()-steps_last))*dt_smooth; // instantaneous time estimation
}
void Info::print_logo() const {
	const int a=color_light_blue, b=color_orange, c=color_pink;
	print(".-----------------------------------------------------------------------------.\n");
	print("|                       "); print("______________   ", a);               print("______________", b);      print("                       |\n");
	print("|                       "); print("\\   ________  | ", a);               print("|  ________   /", b);     print("                       |\n");
	print("|                        "); print("\\  \\       | | ", a);              print("| |       /  /", b);     print("                        |\n");
	print("|                         "); print("\\  \\      | | ", a);              print("| |      /  /", b);     print("                         |\n");
	print("|                          "); print("\\  \\     | | ", a);              print("| |     /  /", b);     print("                          |\n");
	print("|                           "); print("\\  \\_.-\"  | ", a);             print("|  \"-._/  /", b);    print("                           |\n");
	print("|                            "); print("\\    _.-\" ", a);  print("_ ", c); print("\"-._    /", b);  print("                            |\n");
	print("|                             "); print("\\.-\" ", a); print("_.-\" \"-._ ", c); print("\"-./", b); print("                             |\n");
	print("|                               ");                print(".-\"  .-\"-.  \"-.", c);                print("                               |\n");
	print("|                               ");                print("\\  v\"     \"v  /", c);                print("                               |\n");
	print("|                                ");                print("\\  \\     /  /", c);                 print("                                |\n");
	print("|                                 ");                print("\\  \\   /  /", c);                 print("                                 |\n");
	print("|                                  ");                print("\\  \\ /  /", c);                 print("                                  |\n");
	print("|                                   ");                print("\\  '  /", c);                  print("                                   |\n");
	print("|                                    ");                print("\\   /", c);                  print("                                    |\n");
	print("|                                     ");                print("\\ /", c);                  print("                FluidX3D Version 2.7 |\n");
	print("|                                      ");                 print("'", c);                  print("     Copyright (c) Dr. Moritz Lehmann |\n");
}
void Info::print_initialize() {
	const float Re = lbm->get_Re_max();
	println("|-----------------.-----------------------------------------------------------|");
	println("| Grid Resolution | "+alignr(57u, to_string(lbm->get_Nx())+" x "+to_string(lbm->get_Ny())+" x "+to_string(lbm->get_Nz())+" = "+to_string(lbm->get_N()))+" |");
	println("| Grid Domains    | "+alignr(57u, to_string(lbm->get_Dx())+" x "+to_string(lbm->get_Dy())+" x "+to_string(lbm->get_Dz())+" = "+to_string(lbm->get_D()))+" |");
	println("| LBM Type        | "+alignr(57u, /***************/ "D"+to_string(lbm->get_velocity_set()==9?2:3)+"Q"+to_string(lbm->get_velocity_set())+" "+collision)+" |");
	println("| Memory Usage    | "+alignr(54u, /*******/ "CPU "+to_string(cpu_mem_required)+" MB, GPU "+to_string(lbm->get_D())+"x "+to_string(gpu_mem_required))+" MB |");
	println("| Max Alloc Size  | "+alignr(54u, /*************/ (uint)(lbm->get_N()/(ulong)lbm->get_D()*(ulong)(lbm->get_velocity_set()*sizeof(fpxx))/1048576ull))+" MB |");
	println("| Time Steps      | "+alignr(57u, /***************************************************************/ (steps==max_ulong ? "infinite" : to_string(steps)))+" |");
	println("| Kin. Viscosity  | "+alignr(57u, /*************************************************************************************/ to_string(lbm->get_nu(), 8u))+" |");
	println("| Relaxation Time | "+alignr(57u, /************************************************************************************/ to_string(lbm->get_tau(), 8u))+" |");
	println("| Reynolds Number | "+alignr(57u, /******************************************/ "Re < "+string(Re>=100.0f ? to_string(to_uint(Re)) : to_string(Re, 6u)))+" |");
#ifdef VOLUME_FORCE
	println("| Volume Force    | "+alignr(57u, alignr(15u, to_string(lbm->get_fx(), 8u))+","+alignr(15u, to_string(lbm->get_fy(), 8u))+","+alignr(15u, to_string(lbm->get_fz(), 8u)))+" |");
#endif // VOLUME_FORCE
#ifdef SURFACE
	println("| Surface Tension | "+alignr(57u, /**********************************************************************************/ to_string(lbm->get_sigma(), 8u))+" |");
#endif // SURFACE
#ifdef TEMPERATURE
	println("| Thermal Diff.   | "+alignr(57u, /**********************************************************************************/ to_string(lbm->get_alpha(), 8u))+" |");
	println("| Thermal Exp.    | "+alignr(57u, /***********************************************************************************/ to_string(lbm->get_beta(), 8u))+" |");
#endif // TEMPERATURE
#ifndef INTERACTIVE_GRAPHICS_ASCII
	println("|---------.-------'-----.-----------.-------------------.---------------------|");
	println("| MLUPs   | Bandwidth   | Steps/s   | Current Step      | "+string(steps==max_ulong?"Elapsed Time  ":"Time Remaining")+"      |");
#else // INTERACTIVE_GRAPHICS_ASCII
	println("'-----------------'-----------------------------------------------------------'");
#endif // INTERACTIVE_GRAPHICS_ASCII
	allow_rendering = true;
}
void Info::print_update() const {
	if(allow_rendering) reprint(
		"|"+alignr(8, to_uint((double)lbm->get_N()*1E-6/dt_smooth))+" |"+ // MLUPs
		alignr(7, to_uint((double)lbm->get_N()*(double)device_transfer*1E-9/dt_smooth))+" GB/s |"+ // memory bandwidth
		alignr(10, to_uint(1.0/dt_smooth))+" | "+ // steps/s
		(steps==max_ulong ? alignr(17, lbm->get_t()) : alignr(12, lbm->get_t())+" "+print_percentage((double)(lbm->get_t()-steps_last)/(double)steps))+" | "+ // current step
		alignr(19, print_time(time()))+" |" // either elapsed time or remaining time
	);
#ifdef GRAPHICS
	if(key_G) { // print camera settings
		const string camera_position = "float3("+to_string(camera.pos.x/(float)lbm->get_Nx(), 6u)+"f*(float)Nx, "+to_string(camera.pos.y/(float)lbm->get_Ny(), 6u)+"f*(float)Ny, "+to_string(camera.pos.z/(float)lbm->get_Nz(), 6u)+"f*(float)Nz)";
		const string camera_rx_ry_fov = to_string(degrees(camera.rx)-90.0, 1u)+"f, "+to_string(180.0-degrees(camera.ry), 1u)+"f, "+to_string(camera.fov, 1u)+"f";
		const string camera_zoom = to_string(camera.zoom*(float)fmax(fmax(lbm->get_Nx(), lbm->get_Ny()), lbm->get_Nz())/(float)min(camera.width, camera.height), 6u)+"f";
		if(camera.free) print_info("lbm.graphics.set_camera_free("+camera_position+", "+camera_rx_ry_fov+");");
		else print_info("lbm.graphics.set_camera_centered("+camera_rx_ry_fov+", "+camera_zoom+");");
		key_G = false;
	}
#endif // GRAPHICS
}
void Info::print_finalize() {
	allow_rendering = false;
	println("\n|---------'-------------'-----------'-------------------'---------------------|");
}