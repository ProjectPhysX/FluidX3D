#include "info.hpp"
#include "lbm.hpp"
#include "setup.hpp"

#ifdef GRAPHICS
void draw_scale(const int field_mode, const int color) {
	float scale_min=0.0f, scale_max=1.0f;
	string title = "";
	const int label_count = 10;
	switch(field_mode) {
		case 0: // coloring by velocity
			scale_min = 0.0f;
			scale_max = units.si_u(1.0f)==1.0f ? (GRAPHICS_U_MAX) : units.si_u(0.57735027f*(GRAPHICS_U_MAX));
			title = "velocity u / "+string(units.si_u(1.0f)==1.0f ? "c" : "[m/s]");
			break;
		case 1: // coloring by density
			scale_min = units.si_rho(1.0f)-units.si_rho(GRAPHICS_RHO_DELTA);
			scale_max = units.si_rho(1.0f)+units.si_rho(GRAPHICS_RHO_DELTA);
			title = "density rho / "+string(units.si_u(1.0f)==1.0f ? "1" : "[kg/m^3]");
			break;
		case 2: // coloring by temperature
			scale_min = units.si_T(1.0f)-units.si_T(GRAPHICS_T_DELTA);
			scale_max = units.si_T(1.0f)+units.si_T(GRAPHICS_T_DELTA);
			title = "temperature T / "+string(units.si_u(1.0f)==1.0f ? "1" : "[K]");
			break;
	}
	const int margin_x=2*(FONT_WIDTH), margin_y=1*(FONT_HEIGHT); // margins in x and y
	const int ox=camera.width-16*(FONT_WIDTH)-margin_x-1, oy=(FONT_HEIGHT)*3/2+margin_y+8; // plot area offset x/y
	const int label_length_max = scale_max>=1000.0f ? (int)length(to_string(to_int(scale_max))) : 5;
	const int w = camera.width-(ox+margin_x+label_length_max*(FONT_WIDTH)+8); // width of color scale
	const int h = camera.height-(FONT_HEIGHT)*3/2-12*(FONT_HEIGHT)-2*margin_y-4; // height of color scale
	const int N = min(((h/(FONT_HEIGHT))/2)*2, label_count); // number of labels on the y-axis
	for(int i=0; i<h; i++) {
		const float v = (float)i/(float)h;
		int c = 0;
		switch(field_mode) {
			case 0: c = colorscale_rainbow(v); break; // coloring by velocity
			case 1: c = colorscale_twocolor(v); break; // coloring by density
			case 2: c = colorscale_iron(v); break; // coloring by temperature
		}
		draw_line_label(ox, oy+h-i, ox+w, oy+h-i, c);
	}
	draw_line_label(ox  , oy+h, ox+w+1, oy+h  , color); // x-axis
	draw_line_label(ox  , oy  , ox    , oy+h+1, color); // y-axis
	draw_line_label(ox  , oy  , ox+w+1, oy    , color); // x-axis mirror
	draw_line_label(ox+w, oy  , ox+w  , oy+h+1, color); // y-axis mirror
	for(int i=0; i<=N; i++) {
		const float f = (float)i/(float)N;
		const float v = scale_min+f*(scale_max-scale_min);
		const string ly = scale_max>=1000.0f ? to_string(to_int(v)) : scale_max>=100.0f ? to_string(v, 1u) : scale_max>=10.0f ? to_string(v, 2u) : to_string(v, 3u);
		const int y = h-to_int(h*f);
		draw_line_label(ox, oy+y, ox+4, oy+y, color); // y-axis tickmarks
		draw_line_label(ox+w-3, oy+y, ox+w+4, oy+y, color); // y-axis mirror tickmarks
		draw_label(ox+w+7, oy+y-(FONT_HEIGHT)/2, ly, color); // y-axis labels
	}
	draw_label(ox+min(0, w+7+label_length_max*(FONT_WIDTH)-(int)length(title)*(FONT_WIDTH)), oy-(FONT_HEIGHT)*3/2-6, title, color); // colorbar title
}
void main_label(const double frametime) {
	if(info.allow_rendering) {
		info.print_update();
		const int c = invert(GRAPHICS_BACKGROUND_COLOR);
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
		draw_label(2, camera.height-1*(FONT_HEIGHT)-1, "FluidX3D v2.14 Copyright (c) Dr. Moritz Lehmann", c);
		if(!key_H) {
			draw_label(camera.width-16*(FONT_WIDTH)-1, 2, "Press H for Help", c);
		} else {
			if(info.lbm->graphics.visualization_modes&(VIS_FIELD|VIS_STREAMLINES|VIS_Q_CRITERION)) draw_scale(info.lbm->graphics.field_mode, c);
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

			const int sl=info.lbm->graphics.slice_mode, fl=info.lbm->graphics.field_mode;
			const string sx="x="+alignr(4u, info.lbm->graphics.slice_x), sy="y="+alignr(4u, info.lbm->graphics.slice_y), sz="z="+alignr(4u, info.lbm->graphics.slice_z);
			string slice = sl==0 ? "      disabled      " : sl==1 ? sx+"|      |      " : sl==2 ? "      |"+sy+"|      " : sl==3 ? "      |      |"+sz : sl==4 ? sx+"|      |"+sz : sl==5 ? sx+"|"+sy+"|"+sz : sl==6 ? "      |"+sy+"|"+sz : sx+"|"+sy+"|      ";
			string field = fl==0 ? "     velocity u     " : fl==1 ? "     density rho    " : "    temperature T   ";

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
			draw_label(ox, oy+i, "Z: ("+field+"): toggle field visualization mode", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "Q/E: move slice in slice visualization mode", c); i+=2*FONT_HEIGHT;
			draw_label(ox, oy+i, "Mouse or I/J/K/L (rx="+alignr(4u, to_int(fmod(degrees(camera.rx)+90.0+360.0, 360.0)-180.0))+", ry="+alignr(3u, to_int(180.0-degrees(camera.ry)))+"): rotate camera", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "Scrollwheel or +/- ("+to_string(camera.free ? (float)camera.free_camera_velocity : camera.zoom*(float)fmax(fmax(info.lbm->get_Nx(), info.lbm->get_Ny()), info.lbm->get_Nz())/(float)min(camera.width, camera.height), 3u)+"): zoom (centered camera mode) or camera movement speed (free camera mode)", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "Mouseclick or U: toggle rotation with Mouse and angle snap rotation with I/J/K/L", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "Y/X ("+alignr(3u, to_int(camera.fov))+"): adjust camera field of view", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "F ("+string(camera.free?"  free  ":"centered")+"): toggle centered/free camera mode", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "W/A/S/D/Space/C ("+to_string(camera.pos.x/(float)info.lbm->get_Nx(), 2u)+"*Nx, "+to_string(camera.pos.y/(float)info.lbm->get_Ny(), 2u)+"*Ny, "+to_string(camera.pos.z/(float)info.lbm->get_Nz(), 2u)+"*Nz): move free camera", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "G: print current camera position/rotation in console as copy/paste command", c); i+=FONT_HEIGHT;
			draw_label(ox, oy+i, "R ("+string(camera.autorotation?" active ":"inactive")+"): toggle camera autorotation", c); i+=2*FONT_HEIGHT;
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