#include "info.hpp"
#include "lbm.hpp"
#include "setup.hpp"

#ifdef GRAPHICS
void main_label(const double frametime) {
	if(info.allow_rendering) {
		info.print_update();
		const int c = color(255-red(GRAPHICS_BACKGROUND_COLOR), 255-green(GRAPHICS_BACKGROUND_COLOR), 255-blue(GRAPHICS_BACKGROUND_COLOR));
		const int ox=-36*FONT_WIDTH-1, oy=-11*FONT_HEIGHT-1;
		int i = 0;
		const float Re = info.lbm->get_Re_max();
		const double pn=(double)info.lbm->get_N(), mt=(double)info.device_transfer;
		draw_label(camera.width+ox, camera.height+oy+i, "Resolution "     +alignr(25u, to_string(info.lbm->get_Nx())+"x"+to_string(info.lbm->get_Ny())+"x"+to_string(info.lbm->get_Nz())+" = "+to_string(info.lbm->get_N())), c); i+=FONT_HEIGHT;
		//draw_label(camera.width+ox, camera.height+oy+i, "Volume Force "   +alignr(15u,                          info.lbm->get_fx())+","+alignr(15, info.lbm->get_fy())+", "+alignr(15, info.lbm->get_fz()), c); i+=FONT_HEIGHT;
		draw_label(camera.width+ox, camera.height+oy+i, "Kin. Viscosity " +alignr(21u,                                                                                  to_string(info.lbm->get_nu(), 8u)), c); i+=FONT_HEIGHT;
		draw_label(camera.width+ox, camera.height+oy+i, "Relaxation Time "+alignr(20u,                                                                                 to_string(info.lbm->get_tau(), 8u)), c); i+=FONT_HEIGHT;
		draw_label(camera.width+ox, camera.height+oy+i, "Reynolds Number "+alignr(20u,                                            "Re < "+string(Re>=100.0f ? to_string(to_uint(Re)) : to_string(Re, 6u))), c); i+=FONT_HEIGHT;
		draw_label(camera.width+ox, camera.height+oy+i, "LBM Type "       +alignr(27u, "D"+to_string(info.lbm->get_velocity_set()==9u?2:3)+"Q"+to_string(info.lbm->get_velocity_set())+" "+info.collision), c); i+=FONT_HEIGHT;
		draw_label(camera.width+ox, camera.height+oy+i, "RAM Usage "      +alignr(26u,                         "CPU "+to_string(info.cpu_mem_required)+" MB, GPU "+to_string(info.gpu_mem_required)+" MB"), c); i+=FONT_HEIGHT;
		draw_label(camera.width+ox, camera.height+oy+i, (info.steps==max_ulong ? "Elapsed Time   " : "Remaining Time ")+alignr(21u,                                               print_time(info.time())), c); i+=FONT_HEIGHT;
		draw_label(camera.width+ox, camera.height+oy+i, "Simulation Time "+alignr(20u,             (units.si_t(1ull)==1.0f?to_string(info.lbm->get_t()):to_string(units.si_t(info.lbm->get_t()), 6u))+"s"), c); i+=FONT_HEIGHT;
		draw_label(camera.width+ox, camera.height+oy+i, "MLUPs "          +alignr(30u,        alignr(5u, to_uint(pn*1E-6/info.dt_smooth))+" ("+alignr(5u, to_uint(pn*mt*1E-9/info.dt_smooth))+"    GB/s)"), c); i+=FONT_HEIGHT;
		draw_label(camera.width+ox, camera.height+oy+i, "Steps "          +alignr(30u,                             alignr(10u, info.lbm->get_t())+" ("+alignr(5, to_uint(1.0/info.dt_smooth))+" Steps/s)"), c); i+=FONT_HEIGHT;
		draw_label(camera.width+ox, camera.height+oy+i, "FPS "            +alignr(32u,                                   alignr(4u, to_uint(1.0/frametime))+" ("+alignr(5u, camera.fps_limit)+" fps max)"), c);
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