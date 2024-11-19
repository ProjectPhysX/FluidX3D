#include "info.hpp"
#include "lbm.hpp"

Info info;

void Info::append(const ulong steps, const ulong total_steps, const ulong t) {
	if(total_steps==max_ulong) { // total_steps is not provided/used
		this->steps = steps; // has to be executed before info.print_initialize()
		this->steps_last = t; // reset last step count if multiple run() commands are executed consecutively
		this->runtime_total_last = this->runtime_total; // reset last runtime if multiple run() commands are executed consecutively
		this->runtime_total = clock.stop();
	} else { // total_steps has been specified
		this->steps = total_steps; // has to be executed before info.print_initialize()
	}
}
void Info::update(const double dt) {
	this->runtime_lbm_timestep_last = dt; // exact dt
	this->runtime_lbm_timestep_smooth = (dt+0.3)/(0.3/runtime_lbm_timestep_smooth+1.0); // smoothed dt
	this->runtime_lbm += dt; // skip first step since it is likely slower than average
	this->runtime_total = clock.stop();
}
double Info::time() const { // returns either elapsed time or remaining time
	if(lbm==nullptr) return 0.0;
	return steps==max_ulong ? runtime_total : ((double)steps/(double)max(lbm->get_t()-steps_last, 1ull)-1.0)*(runtime_total-runtime_total_last); // time estimation on average so far
	//return steps==max_ulong ? runtime_lbm : ((double)steps-(double)(lbm->get_t()-steps_last))*runtime_lbm_timestep_smooth; // instantaneous time estimation
}
void Info::print_logo() const {
	const int a=color_light_blue, b=color_orange, c=color_pink;
	print(".-----------------------------------------------------------------------------.\n");
	print("|                      "); print(  " ______________  ", a);                  print(" ______________ ", b); print("                      |\n");
	print("|                       "); print( "\\   ________  | ", a);                  print("|  ________   /", b); print("                       |\n");
	print("|                        "); print("\\  \\       | | ", a);                  print("| |       /  /", b); print("                        |\n");
	print("|                         "); print("\\  \\      | | ", a);                  print("| |      /  /", b); print("                         |\n");
	print("|                          "); print("\\  \\     | | ", a);                  print("| |     /  /", b); print("                          |\n");
	print("|                           "); print("\\  \\_.-\"  | ", a);                print("|  \"-._/  /", b); print("                           |\n");
	print("|                            "); print("\\    _.-\" ", a);  print("_ ", c);  print("\"-._    /", b); print("                            |\n");
	print("|                             "); print("\\.-\" ", a); print("_.-\" \"-._ ", c); print("\"-./", b); print("                             |\n");
	print("|                              ");                 print(" .-\"  .-\"-.  \"-. ", c);               print("                              |\n");
	print("|                               ");                 print("\\  v\"     \"v  /", c);               print("                               |\n");
	print("|                                ");                 print("\\  \\     /  /", c);                print("                                |\n");
	print("|                                 ");                 print("\\  \\   /  /", c);                print("                                 |\n");
	print("|                                  ");                 print("\\  \\ /  /", c);                print("                                  |\n");
	print("|                                   ");                 print("\\  '  /", c);                 print("                                   |\n");
	print("|                                    ");                 print("\\   /", c);                 print("                                    |\n");
	print("|                                     ");                 print("\\ /", c);                 print("                FluidX3D Version 3.0 |\n");
	print("|                                      ");                 print( "'", c);                 print("     Copyright (c) Dr. Moritz Lehmann |\n");
	print("|-----------------------------------------------------------------------------|\n");
}
void Info::print_initialize(LBM* lbm) {
	info.allow_printing.lock(); // disable print_update() until print_initialize() has finished
	this->lbm = lbm;
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
	bool all_domains_use_ram = true; // reset cpu/gpu_mem_required to get valid values for consecutive simulations
	for(uint d=0u; d<lbm->get_D(); d++) {
		all_domains_use_ram = all_domains_use_ram&&lbm->lbm_domain[d]->get_device().info.uses_ram;
	}
	if(all_domains_use_ram) {
		cpu_mem_required = lbm->get_D()*lbm->lbm_domain[0]->get_device().info.memory_used;
		gpu_mem_required = 0u;
	} else {
		cpu_mem_required = (uint)(lbm->get_N()*(ulong)bytes_per_cell_host()/1048576ull);
		gpu_mem_required = lbm->lbm_domain[0]->get_device().info.memory_used;
	}
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
	clock.start();
	info.allow_printing.unlock();
}
void Info::print_update() const {
	if(lbm==nullptr) return;
	info.allow_printing.lock();
	reprint(
		"|"+alignr(8, to_uint((double)lbm->get_N()*1E-6/runtime_lbm_timestep_smooth))+" |"+ // MLUPs
		alignr(7, to_uint((double)lbm->get_N()*(double)bandwidth_bytes_per_cell_device()*1E-9/runtime_lbm_timestep_smooth))+" GB/s |"+ // memory bandwidth
		alignr(10, to_uint(1.0/runtime_lbm_timestep_smooth))+" | "+ // steps/s
		(steps==max_ulong ? alignr(17, lbm->get_t()) : alignr(12, lbm->get_t())+" "+print_percentage((float)(lbm->get_t()-steps_last)/(float)steps))+" | "+ // current step
		alignr(19, print_time(time()))+" |" // either elapsed time or remaining time
	);
#ifdef GRAPHICS
	if(key_G) { // print camera settings
		const string camera_position = "float3("+alignr(9u, to_string(camera.pos.x/(float)lbm->get_Nx(), 6u))+"f*(float)Nx, "+alignr(9u, to_string(camera.pos.y/(float)lbm->get_Ny(), 6u))+"f*(float)Ny, "+alignr(9u, to_string(camera.pos.z/(float)lbm->get_Nz(), 6u))+"f*(float)Nz)";
		const string camera_rx_ry_fov = alignr(6u, to_string(degrees(camera.rx)-90.0, 1u))+"f, "+alignr(5u, to_string(180.0-degrees(camera.ry), 1u))+"f, "+alignr(5u, to_string(camera.fov, 1u))+"f";
		const string camera_zoom = alignr(8u, to_string(camera.zoom*(float)fmax(fmax(lbm->get_Nx(), lbm->get_Ny()), lbm->get_Nz())/(float)min(camera.width, camera.height), 6u))+"f";
		if(camera.free) println("\rlbm.graphics.set_camera_free("+camera_position+", "+camera_rx_ry_fov+");");
		else println("\rlbm.graphics.set_camera_centered("+camera_rx_ry_fov+", "+camera_zoom+");          ");
		key_G = false;
	}
#endif // GRAPHICS
	info.allow_printing.unlock();
}
void Info::print_finalize() {
	lbm = nullptr;
	println("\n|---------'-------------'-----------'-------------------'---------------------|");
}