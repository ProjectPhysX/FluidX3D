#pragma once

#define WINDOW_NAME "FluidX3D"
//#define INTERACTIVE_GRAPHICS
//#define INTERACTIVE_GRAPHICS_ASCII
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

class Camera {
public:
	int* bitmap = nullptr;
	int* zbuffer = nullptr;
	uint width = 1920u; // screen width
	uint height = 1080u; // screen height
	uint fps_limit = 60u; // default value for screen frames per second limit
	float fov = 100.0f; // field of view, default: 100
	float zoom=0.5f*(float)min(width, height), dis=0.5f*(float)width/tan(fov*pif/360.0f); // zoom, distance from camera to rotation center
	float3x3 R = float3x3(1.0f); // camera rotation matrix
	double rx=0.5*pi, ry=pi; // rotation angles
	float3 pos = float3(0.0f); // free camera position
	bool free = false; // free camera mode
	double free_camera_velocity = 0.05; // free camera speed
	bool vr=false, tv=false; // virtual reality mode (enables stereoscopic rendering), VR TV mode
	float eye_distance = 8.0f; // distance between cameras in VR mode
	bool autorotation = false; // autorotation
	bool key_update = true; // a key variable has been updated

private:
	float log_zoom=4.0f*log(zoom), target_log_zoom=log_zoom;
	double mouse_x=0.0, mouse_y=0.0, target_mouse_x=0.0, target_mouse_y=0.0; // mouse position
	double mouse_sensitivity = 1.0; // mouse sensitivity
	bool lockmouse = false; // mouse movement won't change camera view when this is true
	bool key_state[512] = { 0 };

public:
	Camera(const uint width, const uint height, const uint fps_limit) {
		this->width = width;
		this->height = height;
		this->fps_limit = fps_limit;
		bitmap = new int[width*height];
		zbuffer = new int[width*height];
		set_zoom(1.0f); // set initial zoom
		update_matrix();
	}
	Camera() = default; // default constructor
	~Camera() {
		delete[] bitmap;
		delete[] zbuffer;
	}
	Camera& operator=(Camera&& camera) noexcept { // move assignment
		this->width = camera.width;
		this->height = camera.height;
		this->fps_limit = camera.fps_limit;
		std::swap(bitmap, camera.bitmap);
		std::swap(zbuffer, camera.zbuffer);
		set_zoom(1.0f); // set initial zoom
		update_matrix();
		return *this;
	}

	void set_zoom(const float rad) {
		zoom = 0.5f*(float)min(width, height)/rad;
		log_zoom = target_log_zoom = 4.0f*log(zoom);
	}
	void update_matrix() {
		dis = 0.5f*(float)width/tan(fov*pif/360.0f);
		const float sinrx=sin((float)rx), cosrx=cos((float)rx), sinry=sin((float)ry), cosry=cos((float)ry);
		R.xx =  cosrx;       R.xy =  sinrx;       R.xz = 0.0f;
		R.yx =  sinrx*sinry; R.yy = -cosrx*sinry; R.yz = cosry;
		R.zx = -sinrx*cosry; R.zy =  cosrx*cosry; R.zz = sinry;
	}
	void set_key_state(const int key, const bool state) {
		key_state[clamp(256+key, 0, 511)] = state;
	}
	bool get_key_state(const int key) {
		return key_state[clamp(256+key, 0, 511)];
	}
	void input_key(const int key) {
		switch(key) {
			case 'R': input_R(); break;
			case 'U': input_U(); break;
			case 'I': input_I(); break;
			case 'J': input_J(); break;
			case 'K': input_K(); break;
			case 'L': input_L(); break;
			case 'V': input_V(); break;
			case 'B': input_B(); break;
			case '+': input_scroll_down(); break;
			case '-': input_scroll_up(); break;
			case 'F': input_F(); break;
			case 27: running=false; exit(0);
		}
#ifdef INTERACTIVE_GRAPHICS_ASCII
		if(free) { // move free camera
			if(key=='W') input_W();
			if(key=='A') input_A();
			if(key=='S') input_S();
			if(key=='D') input_D();
			if(key==' ') input_Space();
			if(key=='C') input_C();
		}
		if(!lockmouse) {
			if(key=='I') input_I(); // rotating camera with keys
			if(key=='J') input_J();
			if(key=='K') input_K();
			if(key=='L') input_L();
		}
		if(key=='Y') input_Y(); // adjusting field of view
		if(key=='X') input_X();
		if(key=='N') input_N(); // adjust camera.vr eye distance
		if(key=='M') input_M();
#endif // INTERACTIVE_GRAPHICS_ASCII
	}
	void update_state() {
		if(!free) {
			log_zoom = 0.8f*log_zoom+0.2f*target_log_zoom; // continuous zoom
			zoom = exp(log_zoom*0.25f);
		} else { // move free camera
			if(get_key_state('W')) input_W();
			if(get_key_state('A')) input_A();
			if(get_key_state('S')) input_S();
			if(get_key_state('D')) input_D();
			if(get_key_state(' ')) input_Space();
			if(get_key_state('C')) input_C();
		}
		if(!lockmouse) {
			if(get_key_state('I')) input_I(); // rotate camera with keys
			if(get_key_state('J')) input_J();
			if(get_key_state('K')) input_K();
			if(get_key_state('L')) input_L();
		}
		if(autorotation) update_rotation(-1, 0);
		if(get_key_state('Y')) input_Y(); // adjust field of view
		if(get_key_state('X')) input_X();
		if(get_key_state('N')) input_N(); // adjust vr eye distance
		if(get_key_state('M')) input_M();
		mouse_x = 0.8*mouse_x+0.2*target_mouse_x; // continuous mouse movement
		mouse_y = 0.8*mouse_y+0.2*target_mouse_y;
		target_mouse_x = 0.0;
		target_mouse_y = 0.0;
		if(!lockmouse) update_rotation(mouse_x, mouse_y);
	}
	void clear_frame() {
		for(uint i=0u; i<width*height; i++) {
			bitmap[i] = GRAPHICS_BACKGROUND_COLOR;
			zbuffer[i] = min_int;
		}
	}
	float data(const uint i) const { // returns all camera data required for rendering
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

	void input_mouse_moved(const int x, const int y) {
		if(!lockmouse) {
			target_mouse_x = mouse_sensitivity*(double)((int)width /2-x);
			target_mouse_y = mouse_sensitivity*(double)((int)height/2-y);
		}
	}
	void input_mouse_dragged(const int dx, const int dy) {
		if(!lockmouse) {
			target_mouse_x -= mouse_sensitivity*(double)(dx);
			target_mouse_y -= mouse_sensitivity*(double)(dy);
		}
	}
	void input_scroll_up() {
		if(!free) { // zoom
			target_log_zoom -= 1.0f;
		} else if(!lockmouse) {
			free_camera_velocity *= 1.284;
		}
	}
	void input_scroll_down() {
		if(!free) { // zoom
			target_log_zoom += 1.0f;
		} else if(!lockmouse) {
			free_camera_velocity /= 1.284;
		}
	}

private:
	void input_F() {
		free = !free;
		if(!free) {
			zoom = exp(log_zoom*0.25f);
		} else {
			zoom = 1E16f;
		}
	}
	void input_V() {
		vr = !vr;
	}
	void input_B() {
		tv = !tv;
	}
	void input_W() {
		pos.x += R.xy*R.yz*(float)free_camera_velocity;
		pos.y -= R.xx*R.yz*(float)free_camera_velocity;
		pos.z -= R.zz*(float)free_camera_velocity;
	}
	void input_A() {
		pos.x -= R.xx*(float)free_camera_velocity;
		pos.y -= R.xy*(float)free_camera_velocity;
	}
	void input_S() {
		pos.x -= R.xy*R.yz*(float)free_camera_velocity;
		pos.y += R.xx*R.yz*(float)free_camera_velocity;
		pos.z += R.zz*(float)free_camera_velocity;
	}
	void input_D() {
		pos.x += R.xx*(float)free_camera_velocity;
		pos.y += R.xy*(float)free_camera_velocity;
	}
	void input_Space() {
		pos.x -= R.xy*R.zz*(float)free_camera_velocity;
		pos.y += R.xx*R.zz*(float)free_camera_velocity;
		pos.z -= R.yz*(float)free_camera_velocity;
	}
	void input_C() {
		pos.x += R.xy*R.zz*(float)free_camera_velocity;
		pos.y -= R.xx*R.zz*(float)free_camera_velocity;
		pos.z += R.yz*(float)free_camera_velocity;
	}
	void input_R() {
		autorotation = !autorotation;
	}
	void input_U() {
#if defined(INTERACTIVE_GRAPHICS) && defined(_WIN32)
		if(!lockmouse) {
			ShowCursor(true); // show cursor
		} else {
			ShowCursor(false); // hide cursor
			SetCursorPos(width/2, height/2); // reset mouse
		}
#endif // Windows
		lockmouse = !lockmouse;
	}
	void input_I() {
		if(lockmouse) {
			double d = (ry*18.0/pi)-(double)((int)(ry*18.0f/pi));
			d = d<1E-6 ? 1.0 : 1.0-d;
			update_rotation(0.0, 10.0*d);
		} else {
			target_mouse_y += mouse_sensitivity;
		}
	}
	void input_J() {
		if(lockmouse) {
			double d = (rx*18.0/pi)-(double)((int)(rx*18.0/pi));
			d = d<1E-6 ? 1.0 : 1.0-d;
			update_rotation(10.0*d, 0.0);
		} else {
			target_mouse_x += mouse_sensitivity;
		}
	}
	void input_K() {
		if(lockmouse) {
			double d = (ry*18.0/pi)-(double)((int)(ry*18.0/pi));
			d = d<1E-6 ? 1.0f : d;
			update_rotation(0.0, -10.0*d);
		} else {
			target_mouse_y -= mouse_sensitivity;
		}
	}
	void input_L() {
		if(lockmouse) {
			double d = (rx*18.0/pi)-(double)((int)(rx*18.0/pi));
			d = d<1E-6 ? 1.0 : d;
			update_rotation(-10.0*d, 0.0);
		} else {
			target_mouse_x -= mouse_sensitivity;
		}
	}
	void input_X() {
		fov = fmax(fov-1.0f, 1.0f);
		dis = 0.5f*(float)width/tan(fov*pif/360.0f);
	}
	void input_Y() {
		fov = fmin(fov<1.0f ? 1.0f : fov+1.0f, 179.0f);
		dis = 0.5f*(float)width/tan(fov*pif/360.0f);
	}
	void input_N() {
		eye_distance = fmax(eye_distance-0.2f, 0.0f);
	}
	void input_M() {
		eye_distance += 0.2f;
	}

	void update_rotation(const double arx, const double ary) {
#if defined(INTERACTIVE_GRAPHICS)&&defined(_WIN32)
		if(!lockmouse) SetCursorPos((int)width/2, (int)height/2);
#endif // INTERACTIVE_GRAPHICS && Windows
		rx += arx*pi/180.0;
		ry += ary*pi/180.0;
		rx = fmod(rx, 2.0*pi);
		ry = clamp(ry, 0.5*pi, 1.5*pi);
		update_matrix();
	}
};

extern Camera camera;
extern bool key_E, key_G, key_H, key_O, key_P, key_Q, key_T, key_Z; // defined in graphics.cpp
extern bool key_1, key_2, key_3, key_4, key_5, key_6, key_7, key_8, key_9, key_0; // defined in graphics.cpp

#define GRAPHICS_CONSOLE // open console additionally to graphics window
#define FONT_HEIGHT 11 // default: 11
#define FONT_WIDTH 6 // default: 6

void set_light(const uint i, const float3& p);

void draw_bitmap(const int* bitmap);
void draw_label(const int x, const int y, const string& s, const int color);

void draw_pixel(const int x, const int y, const int color); // 2D drawing functions
void draw_circle(const int x, const int y, const int r, const int color);
void draw_line(const int x0, const int y0, const int x1, const int y1, const int color);
void draw_triangle(const int x0, const int y0, const int x1, const int y1, const int x2, const int y2, const int color);
void draw_rectangle(const int x0, const int y0, const int x1, const int y1, const int color);
void draw_text(const int x, const int y, const string& s, const int color);

void draw_pixel(const float3& p, const int color); // 3D drawing functions
void draw_circle(const float3& p, const float r, const int color);
void draw_line(const float3& p0, const float3& p1, const int color);
void draw_triangle(const float3& p0, const float3& p1, const float3& p2, const int color, const bool translucent=false);
void draw_triangle(const float3& p0, const float3& p1, const float3& p2, const int c0, const int c1, const int c2, const bool translucent=false);
void draw_text(const float3& p, const float r, const string& s, const int color);

#endif // GRAPHICS