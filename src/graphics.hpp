#pragma once

#define WINDOW_NAME "FluidX3D"
//#define CONSOLE_GRAPHICS
//#define WINDOWS_GRAPHICS
//#define GRAPHICS

#include "defines.hpp"
#include "utilities.hpp"
#include <atomic>

extern vector<string> main_arguments; // console arguments
extern std::atomic_bool running;

#ifdef GRAPHICS
void main_label(const double frametime); // implement these three
void main_graphics();
void main_physics();

struct Camera {
#ifdef GRAPHICS_WIDTH
	uint width = (uint)(GRAPHICS_WIDTH);
#else // GRAPHICS_WIDTH
	uint width = 1920u; // screen width
#endif // GRAPHICS_WIDTH
#ifdef GRAPHICS_HEIGHT
	uint height = (uint)(GRAPHICS_HEIGHT);
#else // GRAPHICS_HEIGHT
	uint height = 1080u; // screen height
#endif // GRAPHICS_HEIGHT
	uint fps_limit = 60u; // default value for screen frames per second limit
	bool key_update = true; // a key variable has been updated
	float fov = 1E-6f; // field of view, default: orthogonal graphics (fov cannot be exactly 0)
	float zoom = 0.5f*(float)min(width, height); // zoom
	float dis = 0.5f*(float)width/tan(fov*pif/360.0f); // distance from camera to rotation center
	float3 pos = float3(0.0f); // free camera position
	float3x3 R = float3x3(1.0f); // camera rotation matrix
	bool free = false; // free camera mode
	bool vr=false, tv=false; // virtual reality mode (enables stereoscopic rendering), VR TV mode
	float eye_distance = 8.0f; // distance between cameras
	const double ms = 1.0; // mouse sensitivity
	double rx=0.5*pi, ry=pi; // rotation angles
	void update_matrix() {
		dis = 0.5f*(float)width/tan(fov*pif/360.0f);
		const float sinrx=sin((float)rx), cosrx=cos((float)rx), sinry=sin((float)ry), cosry=cos((float)ry);
		R.xx =  cosrx;       R.xy =  sinrx;       R.xz = 0.0f;
		R.yx =  sinrx*sinry; R.yy = -cosrx*sinry; R.yz = cosry;
		R.zx = -sinrx*cosry; R.zy =  cosrx*cosry; R.zz = sinry;
	}
	float data(const uint i) { // returns all camera data required for rendering
		switch(i) {
			case  0: return zoom   ; // camera zoom
			case  1: return dis    ; // distance from camera to rotation center
			case  2: return free ? pos.x : 0.0f; // camera position
			case  3: return free ? pos.y : 0.0f;
			case  4: return free ? pos.z : 0.0f;
			case  5: return R.xx; // camera rotation matrix
			case  6: return R.xy;
			case  7: return R.xz;
			case  8: return R.yx;
			case  9: return R.yy;
			case 10: return R.yz;
			case 11: return R.zx;
			case 12: return R.zy;
			case 13: return R.zz;
			case 14: return as_float((uint)vr<<31|(uint)tv<<30|((uint)float_to_half(eye_distance)&0xFFFF)); // stereoscopic rendering parameters
			default: return 0.0f;
		}
	}
};

extern Camera camera;
extern bool key_E, key_F, key_G, key_H, key_O, key_P, key_Q, key_T, key_Z; // defined in graphics.cpp
extern bool key_1, key_2, key_3, key_4, key_5, key_6, key_7, key_8, key_9, key_0; // defined in graphics.cpp

struct Color {
	uchar r, g, b;
	Color(const int r, const int g, const int b) {
		this->r = (uchar)r;
		this->g = (uchar)g;
		this->b = (uchar)b;
	}
	Color() = default;
};

#define FONT_HEIGHT 10 // default: 10
#define FONT_WIDTH 5 // default: 6

void set_zoom(const float rad);
void set_light(const uint i, const float3& p);

void draw_bitmap(const void* buffer);
void draw_label(const Color& c, const string& s, const int x, const int y);

#ifdef WINDOWS_GRAPHICS

#define GRAPHICS_CONSOLE
//#define SKIP_VISIBILITY_CHECKS // makes CPU graphics without polygons 40% faster

void draw_pixel(const Color& c, const int x, const int y); // 2D drawing functions
void draw_circle(const Color& c, const int x, const int y, const int r);
void draw_line(const Color& c, const int x0, const int y0, const int x1, const int y1);
void draw_triangle(const Color& c, const int x0, const int y0, const int x1, const int y1, const int x2, const int y2);
void draw_rectangle(const Color& c, const int x0, const int y0, const int x1, const int y1);
void draw_polygon(const Color& c, const int* const x, const int* const y, const int n);
void draw_text(const Color& c, const string& s, const int x, const int y);

void draw_pixel(const Color& c, const float3& p); // 3D drawing functions
void draw_circle(const Color& c, const float3& p, const float r);
void draw_line(const Color& c, const float3& p0, const float3& p1);
void draw_triangle(const Color& c, const float3& p0, const float3& p1, const float3& p2, const bool translucent=false);
void draw_quadrangle(const Color& c, const float3& p0, const float3& p1, const float3& p2, const float3& p3, const bool translucent=false);
void draw_text(const Color& c, const float3& p, const string& s, const float r);

#endif // WINDOWS_GRAPHICS
#endif // GRAPHICS