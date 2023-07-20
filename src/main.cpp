#include "info.hpp"
#include "lbm.hpp"
#include "setup.hpp"

#ifdef GRAPHICS
void main_label(const double frametime) {
	if(info.allow_rendering) {
		info.print_update();
		const int c = color(255-red(GRAPHICS_BACKGROUND_COLOR), 255-green(GRAPHICS_BACKGROUND_COLOR), 255-blue(GRAPHICS_BACKGROUND_COLOR));
		{
			const int ox=camera.width-37*(FONT_WIDTH)-1, oy=camera.height-11*(FONT_HEIGHT)-1;
			int i = 0;
			const float Re = info.lbm->get_Re_max();
			const double pn=(double)info.lbm->get_N(), mt=(double)bandwidth_bytes_per_cell_device();
			draw_label(ox, oy+i, "Resolution "     +alignr(26u, /********/ to_string(info.lbm->get_Nx())+"x"+to_string(info.lbm->get_Ny())+"x"+to_string(info.lbm->get_Nz())+" = "+to_string(info.lbm->get_N())), c); i+=FONT_HEIGHT;
			//draw_label(ox, oy+i, "Volume Force "   +alignr(16u, /***************************************************/ info.lbm->get_fx())+","+alignr(15, info.lbm->get_fy())+", "+alignr(15, info.lbm->get_fz()), c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "Kin. Viscosity " +alignr(22u, /***********************************************************************************************************/ to_string(info.lbm->get_nu(), 8u)), c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "Relaxation Time "+alignr(21u, /**********************************************************************************************************/ to_string(info.lbm->get_tau(), 8u)), c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "Reynolds Number "+alignr(21u, /*********************************************************************/ "Re < "+string(Re>=100.0f ? to_string(to_uint(Re)) : to_string(Re, 6u))), c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "LBM Type "       +alignr(28u, /**************************/ "D"+to_string(info.lbm->get_velocity_set()==9u?2:3)+"Q"+to_string(info.lbm->get_velocity_set())+" "+info.collision), c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "Memory "         +alignr(30u, /****************/ "CPU "+to_string(info.cpu_mem_required)+" MB, GPU "+to_string(info.lbm->get_D())+"x "+to_string(info.gpu_mem_required)+" MB"), c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, (info.steps==max_ulong ? "Elapsed Time   " : "Remaining Time ")+alignr(22u, /************************************************************************/ print_time(info.time())), c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "Simulation Time "+alignr(21u, /**************************************/ (units.si_t(1ull)==1.0f?to_string(info.lbm->get_t()):to_string(units.si_t(info.lbm->get_t()), 6u))+"s"), c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "MLUPs "          +alignr(31u, alignr(5u, to_uint(pn*1E-6/info.runtime_lbm_timestep_smooth))+" ("+alignr(5u, to_uint(pn*mt*1E-9/info.runtime_lbm_timestep_smooth))+"    GB/s)"), c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "Steps "          +alignr(31u, /************************************/ alignr(10u, info.lbm->get_t())+" ("+alignr(5, to_uint(1.0/info.runtime_lbm_timestep_smooth))+" Steps/s)"), c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "FPS "            +alignr(33u, /************************************************************/ alignr(4u, to_uint(1.0/frametime))+" ("+alignr(5u, camera.fps_limit)+" fps max)"), c);
		}
		if(!key_H) {
			draw_label(camera.width-16*(FONT_WIDTH)-1, 2, "Press H for Help", c);
		} else {
#ifdef SURFACE
			const bool surface = true;
#else // SURFACE
			const bool surface = false;
#endif // SURFACE
#ifdef PARTICLES
			const bool particles = true;
#else // PARTICLES
			const bool particles = false;
#endif // PARTICLES
			const int ox=2, oy=2;
			int i = 0;

			const int mode = info.lbm->graphics.visualization_modes;
			string mode_1 = (mode&3)==0 ? "inactive" : (mode&3)==VIS_FLAG_LATTICE ? " flags  " : (mode&3)==VIS_FLAG_SURFACE ? " solid  " : "  both  ";
			string mode_2 = mode&VIS_FIELD ? " active " : "inactive";
			string mode_3 = mode&VIS_STREAMLINES ? " active " : "inactive";
			string mode_4 = mode&VIS_Q_CRITERION ? " active " : "inactive";
			string mode_5 = surface ? (mode&VIS_PHI_RASTERIZE ? " active " : "inactive") : "disabled";
			string mode_6 = surface&&info.lbm->get_D()==1u ? (mode&VIS_PHI_RAYTRACE ? " active " : "inactive") : "disabled";
			string mode_7 = particles ? (mode&VIS_PARTICLES ? " active " : "inactive") : "disabled";

			const int sl = info.lbm->graphics.slice_mode;
			const string sx = "x="+alignr(4u, info.lbm->graphics.slice_x);
			const string sy = "y="+alignr(4u, info.lbm->graphics.slice_y);
			const string sz = "z="+alignr(4u, info.lbm->graphics.slice_z);
			string slice = sl==0 ? "      disabled      " : sl==1 ? sx+"|      |      " : sl==2 ? "      |"+sy+"|      " : sl==3 ? "      |      |"+sz : sl==4 ? sx+"|      |"+sz : sl==5 ? sx+"|"+sy+"|"+sz : sl==6 ? "      |"+sy+"|"+sz : sx+"|"+sy+"|      ";

			draw_label(ox, oy+i, "Keyboard/Mouse Controls: ", c); i+=2*FONT_HEIGHT;
			draw_label(ox, oy+i, "P ("+string(key_P?"running ":" paused ")+"): start/pause simulation", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "H ("+string(key_H?" shown  ":" hidden ")+"): show/hide help", c); i+=2*FONT_HEIGHT;
			draw_label(ox, oy+i, "1 ("+mode_1+"): flag wireframe / solid surface (and force vectors on solid cells or surface pressure if the extension is used)", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "2 ("+mode_2+"): velocity field", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "3 ("+mode_3+"): streamlines", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "4 ("+mode_4+"): vorticity / velocity-colored Q-criterion isosurface", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "5 ("+mode_5+"): rasterized free surface", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "6 ("+mode_6+"): raytraced free surface", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "7 ("+mode_7+"): particles", c); i+=2*FONT_HEIGHT;
			draw_label(ox, oy+i, "T: ("+slice+"): toggle slice visualization mode", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "Q/E: move slice in slice visualization mode", c); i+=2*FONT_HEIGHT;
			draw_label(ox, oy+i, "Mouse or I/J/K/L (rx="+alignr(4u, to_int(fmod(degrees(camera.rx)+90.0+360.0, 360.0)-180.0))+", ry="+alignr(3u, to_int(180.0-degrees(camera.ry)))+"): rotate camera", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "Scrollwheel or +/- ("+to_string(camera.free ? (float)camera.free_camera_velocity : camera.zoom*(float)fmax(fmax(info.lbm->get_Nx(), info.lbm->get_Ny()), info.lbm->get_Nz())/(float)min(camera.width, camera.height), 3u)+"): zoom (centered camera mode) or camera movement speed (free camera mode)", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "Mouseclick or U: toggle rotation with Mouse and angle snap rotation with I/J/K/L", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "Y/X ("+alignr(3u, to_int(camera.fov))+"): adjust camera field of view", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "G: print current camera position/rotation in console as copy/paste command", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "R ("+string(camera.autorotation?" active ":"inactive")+"): toggle camera autorotation", c); i+=2*FONT_HEIGHT;
			draw_label(ox, oy+i, "F ("+string(camera.free?"  free  ":"centered")+"): toggle centered/free camera mode", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "W/A/S/D/Space/C ("+to_string(camera.pos.x/(float)info.lbm->get_Nx(), 2u)+"*Nx, "+to_string(camera.pos.y/(float)info.lbm->get_Ny(), 2u)+"*Ny, "+to_string(camera.pos.z/(float)info.lbm->get_Nz(), 2u)+"*Nz): move free camera", c); i+=2*FONT_HEIGHT;
			draw_label(ox, oy+i, "V ("+string(camera.vr?"VR mode":"regular")+"): toggle stereoscopic rendering for VR", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "B ("+string(camera.tv?"TV mode":"goggles")+"): toggle VR-goggles/3D-TV mode for stereoscopic rendering", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "N/M ("+to_string(camera.eye_distance, 1u)+"): adjust eye distance for stereoscopic rendering", c); i+=2*FONT_HEIGHT;
			draw_label(ox, oy+i, "Esc/Alt+F4: quit", c);
		}
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