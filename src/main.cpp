#include "info.hpp"
#include "lbm.hpp"
#include "setup.hpp"

#ifdef GRAPHICS
void main_label(const double frametime) {
	if(info.allow_rendering) {
		info.print_update();
		const Color c = Color(255-red(GRAPHICS_BACKGROUND_COLOR), 255-green(GRAPHICS_BACKGROUND_COLOR), 255-blue(GRAPHICS_BACKGROUND_COLOR));
		const int ox=-36*FONT_WIDTH-2, oy=-11*FONT_HEIGHT-3; // character size 5x10
		int i = 0;
		const float Re = info.lbm->get_Re_max();
		const double pn=(double)info.lbm->get_N(), mt=(double)info.device_transfer;
		draw_label(c, "Resolution "     +alignr(25u, to_string(info.lbm->get_Nx())+"x"+to_string(info.lbm->get_Ny())+"x"+to_string(info.lbm->get_Nz())+" = "+to_string(info.lbm->get_N())), camera.width+ox, camera.height+oy+i); i+=FONT_HEIGHT;
		//draw_label(c, "Volume Force "   +alignr(15u,                        info.lbm->get_fx())+","+alignr(15, info.lbm->get_fy())+", "+alignr(15, info.lbm->get_fz()), camera.width+ox, camera.height+oy+i); i+=FONT_HEIGHT;
		draw_label(c, "Kin. Viscosity " +alignr(21u,                                                                                  to_string(info.lbm->get_nu(), 8u)), camera.width+ox, camera.height+oy+i); i+=FONT_HEIGHT;
		draw_label(c, "Relaxation Time "+alignr(20u,                                                                                 to_string(info.lbm->get_tau(), 8u)), camera.width+ox, camera.height+oy+i); i+=FONT_HEIGHT;
		draw_label(c, "Reynolds Number "+alignr(20u,                                            "Re < "+string(Re>=100.0f ? to_string(to_uint(Re)) : to_string(Re, 6u))), camera.width+ox, camera.height+oy+i); i+=FONT_HEIGHT;
		draw_label(c, "LBM Type "       +alignr(27u, "D"+to_string(info.lbm->get_velocity_set()==9u?2:3)+"Q"+to_string(info.lbm->get_velocity_set())+" "+info.collision), camera.width+ox, camera.height+oy+i); i+=FONT_HEIGHT;
		draw_label(c, "RAM Usage "      +alignr(26u,                         "CPU "+to_string(info.cpu_mem_required)+" MB, GPU "+to_string(info.gpu_mem_required)+" MB"), camera.width+ox, camera.height+oy+i); i+=FONT_HEIGHT;
		draw_label(c, (info.steps==max_ulong ? "Elapsed Time   " : "Remaining Time ")+alignr(21u,                                               print_time(info.time())), camera.width+ox, camera.height+oy+i); i+=FONT_HEIGHT;
		draw_label(c, "Simulation Time "+alignr(20u,             (units.si_t(1ull)==1.0f?to_string(info.lbm->get_t()):to_string(units.si_t(info.lbm->get_t()), 6u))+"s"), camera.width+ox, camera.height+oy+i); i+=FONT_HEIGHT;
		draw_label(c, "MLUPs "          +alignr(30u,        alignr(5u, to_uint(pn*1E-6/info.dt_smooth))+" ("+alignr(5u, to_uint(pn*mt*1E-9/info.dt_smooth))+"    GB/s)"), camera.width+ox, camera.height+oy+i); i+=FONT_HEIGHT;
		draw_label(c, "Steps "          +alignr(30u,                             alignr(10u, info.lbm->get_t())+" ("+alignr(5, to_uint(1.0/info.dt_smooth))+" Steps/s)"), camera.width+ox, camera.height+oy+i); i+=FONT_HEIGHT;
		draw_label(c, "FPS "            +alignr(32u,                                   alignr(4u, to_uint(1.0/frametime))+" ("+alignr(5u, camera.fps_limit)+" fps max)"), camera.width+ox, camera.height+oy+i);
	}
}

void main_graphics() {
	if(info.allow_rendering) draw_bitmap(info.lbm->graphics.draw_frame());
}
#endif // GRAPHICS

void main_physics() {
	info.print_logo();
	main_setup(); // execute setup
	running = false;
	exit(0); // make sure that the program stops
}

#ifndef GRAPHICS
int main(int argc, char* argv[]) {
	main_arguments = get_main_arguments(argc, argv);
	thread compute_thread(main_physics);
	do { // main console loop
		info.print_update();
		sleep(0.050);
	} while(running);
	compute_thread.join();
	return 0;
}
#endif // GRAPHICS