#include "graphics.hpp"

vector<string> main_arguments = vector<string>(); // console arguments
std::atomic_bool running = true;

#ifdef GRAPHICS
Camera camera;

// reserved keys for graphics: W,A,S,D, I,J,K,L, F, R,U, V,B, C,VK_SPACE, Y,X, N,M
//bool key_A=false, key_B=false, key_C=false, key_D=false, key_E=false, key_F=false, key_G=false, key_H=false, key_I=false, key_J=false, key_K=false, key_L=false, key_M=false;
//bool key_N=false, key_O=false, key_P=false, key_Q=false, key_R=false, key_S=false, key_T=false, key_U=false, key_V=false, key_W=false, key_X=false, key_Y=false, key_Z=false;
bool key_E=false, key_G=false, key_H=false, key_O=false, key_P=false, key_Q=false, key_T=false, key_Z=false;
bool key_1=false, key_2=false, key_3=false, key_4=false, key_5=false, key_6=false, key_7=false, key_8=false, key_9=false, key_0=false;

const uint max_light_sources = 100u; // maximal number of light sources
uint light_sources_n = 0u; // number of light sources
float3 light_sources[max_light_sources]; // coordinates of light sources
bool lockmouse = false; // mouse movement won't change camera view when this is true
bool autorotation = false; // autorotation
double mx=0.0, my=0.0, dmx=0.0, dmy=0.0; // mouse position
float zo=0.0f, dzo=0.0f; // zoom and perspective
float fl=0.0f, fvel=0.05f; // free camera speed

void update_rotation(const double arx, const double ary) {
#ifdef WINDOWS_GRAPHICS
	SetCursorPos(camera.width/2, camera.height/2);
#endif // WINDOWS_GRAPHICS
	camera.rx += arx*pi/180.0;
	camera.ry += ary*pi/180.0;
	camera.rx = fmod(camera.rx, 2.0*pi);
	camera.ry = clamp(camera.ry, 0.5*pi, 1.5*pi);
	camera.update_matrix();
}

void key_bindings(const int key) {
	camera.key_update = true;
	switch(key) {
		// reserved keys for graphics: W,A,S,D, I,J,K,L, R,U, V,B, C,VK_SPACE, Y,X, N,M
		//case 'A': key_A = !key_A; break;
		//case 'B': key_B = !key_B; break;
		//case 'C': key_C = !key_C; break;
		//case 'D': key_D = !key_D; break;
		case 'E': key_E = !key_E; break;
		//case 'F': key_F = !key_F; break;
		case 'G': key_G = !key_G; break;
		case 'H': key_H = !key_H; break;
		//case 'I': key_I = !key_I; break;
		//case 'J': key_J = !key_J; break;
		//case 'K': key_K = !key_K; break;
		//case 'L': key_L = !key_L; break;
		//case 'M': key_M = !key_M; break;
		//case 'N': key_N = !key_N; break;
		case 'O': key_O = !key_O; break;
		case 'P': key_P = !key_P; break;
		case 'Q': key_Q = !key_Q; break;
		//case 'R': key_R = !key_R; break;
		//case 'S': key_S = !key_S; break;
		case 'T': key_T = !key_T; break;
		//case 'U': key_U = !key_U; break;
		//case 'V': key_V = !key_V; break;
		//case 'W': key_W = !key_W; break;
		//case 'X': key_X = !key_X; break;
		//case 'Y': key_Y = !key_Y; break;
		case 'Z': key_Z = !key_Z; break;

		case '1': key_1 = !key_1; break;
		case '2': key_2 = !key_2; break;
		case '3': key_3 = !key_3; break;
		case '4': key_4 = !key_4; break;
		case '5': key_5 = !key_5; break;
		case '6': key_6 = !key_6; break;
		case '7': key_7 = !key_7; break;
		case '8': key_8 = !key_8; break;
		case '9': key_9 = !key_9; break;
		case '0': key_0 = !key_0; break;

		case 'R': autorotation = !autorotation; break;
		case 'U': lockmouse = !lockmouse; break;
		case 'I':
			if(lockmouse) {
				double d = (camera.ry*18.0/pi)-(double)((int)(camera.ry*18.0f/pi));
				d = d<1E-6 ? 1.0 : 1.0-d;
				update_rotation(0.0, 10.0*d);
			}
			break;
		case 'J':
			if(lockmouse) {
				double d = (camera.rx*18.0/pi)-(double)((int)(camera.rx*18.0/pi));
				d = d<1E-6 ? 1.0 : 1.0-d;
				update_rotation(10.0*d, 0.0);
			}
			break;
		case 'K':
			if(lockmouse) {
				double d = (camera.ry*18.0/pi)-(double)((int)(camera.ry*18.0/pi));
				d = d<1E-6 ? 1.0f : d;
				update_rotation(0.0, -10.0*d);
			}
			break;
		case 'L':
			if(lockmouse) {
				double d = (camera.rx*18.0/pi)-(double)((int)(camera.rx*18.0/pi));
				d = d<1E-6 ? 1.0 : d;
				update_rotation(-10.0*d, 0.0);
			}
			break;
		case 'V': camera.vr = !camera.vr; break;
		case 'B': camera.tv = !camera.tv; break;
		case ',': case 0xBC: // 0xBC is Windows Virtual-Key Code for ','
			if(!camera.free) { // centered camera zoom
				dzo -= 1.0f;
			} else if(!lockmouse) { // free camera speed
				fl -= 1.0f;
				fvel = 0.05f*exp(fl*0.25f);
			}
			break;
		case '.': case 0xBE: // 0xBE is Windows Virtual-Key Code for '.'
			if(!camera.free) { // centered camera zoom
				dzo += 1.0f;
			} else if(!lockmouse) { // free camera speed
				fl += 1.0f;
				fvel = 0.05f*exp(fl*0.25f);
			}
			break;
		case 'F':
			camera.free = !camera.free;
			if(!camera.free) {
				camera.zoom = exp(zo*0.25f);
			} else {
				camera.zoom = 1E16f;
			}
			break;
		case 27:
			running = false;
			exit(0);
	}

#ifndef WINDOWS_GRAPHICS
	if(camera.free) { // move free camera
		if(key=='W') {
			camera.pos.x += camera.R.xy*camera.R.yz*fvel;
			camera.pos.y -= camera.R.xx*camera.R.yz*fvel;
			camera.pos.z -= camera.R.zz*fvel;
		}
		if(key=='A') {
			camera.pos.x -= camera.R.xx*fvel;
			camera.pos.y -= camera.R.xy*fvel;
		}
		if(key=='S') {
			camera.pos.x -= camera.R.xy*camera.R.yz*fvel;
			camera.pos.y += camera.R.xx*camera.R.yz*fvel;
			camera.pos.z += camera.R.zz*fvel;
		}
		if(key=='D') {
			camera.pos.x += camera.R.xx*fvel;
			camera.pos.y += camera.R.xy*fvel;
		}
		if(key==' ') {
			camera.pos.x -= camera.R.xy*camera.R.zz*fvel;
			camera.pos.y += camera.R.xx*camera.R.zz*fvel;
			camera.pos.z -= camera.R.yz*fvel;
		}
		if(key=='C') {
			camera.pos.x += camera.R.xy*camera.R.zz*fvel;
			camera.pos.y -= camera.R.xx*camera.R.zz*fvel;
			camera.pos.z += camera.R.yz*fvel;
		}
	}
	if(!lockmouse) {
		if(key=='I') dmy += camera.ms; // rotating camera with keys
		if(key=='J') dmx += camera.ms;
		if(key=='K') dmy -= camera.ms;
		if(key=='L') dmx -= camera.ms;
	}
	if(key=='Y') { // adjusting field of view
		camera.fov = fmin(camera.fov<1.0f ? 1.0f : camera.fov+1.0f, 179.0f);
		camera.dis = 0.5f*(float)camera.width/tan(camera.fov*pif/360.0f);
	}
	if(key=='X') {
		camera.fov = fmax(camera.fov-1.0f, 1.0f);
		camera.dis = 0.5f*(float)camera.width/tan(camera.fov*pif/360.0f);
	}
	if(key=='N') { // adjust camera.vr eye distance
		camera.eye_distance -= 0.2f;
		if(camera.eye_distance<0.0f) camera.eye_distance = 0.0f;
	}
	if(key=='M') {
		camera.eye_distance += 0.2f;
	}
#endif // WINDOWS_GRAPHICS
}

void set_zoom(const float rad) {
	camera.zoom = 0.5f*(float)min(camera.width, camera.height)/rad;
	zo = 4.0f*log(camera.zoom);
	dzo = zo;
}
void set_light(const uint i, const float3& position) {
	if(i<max_light_sources) {
		light_sources[i] = position;
		light_sources_n = max(light_sources_n, i+1u);
	}
}

bool convert(int& rx, int& ry, float& rz, const float3& p, const int stereo) { // 3D -> 2D
	float3 t, r;
	t.x = p.x-(camera.free ? camera.pos.x : 0.0f)-(float)stereo*camera.eye_distance/camera.zoom*camera.R.xx; // transformation
	t.y = p.y-(camera.free ? camera.pos.y : 0.0f)-(float)stereo*camera.eye_distance/camera.zoom*camera.R.xy;
	t.z = p.z-(camera.free ? camera.pos.z : 0.0f);
	r.z = camera.R.zx*t.x+camera.R.zy*t.y+camera.R.zz*t.z; // z-position for z-buffer
	const float rs = camera.zoom*camera.dis/(camera.dis-r.z*camera.zoom); // perspective (reciprocal is more efficient)
	if(rs<=0.0f) return false; // point is behins camera
	const float tv = camera.tv&&stereo!=0 ? 0.5f : 1.0f;
	r.x = ((camera.R.xx*t.x+camera.R.xy*t.y+camera.R.xz*t.z)*rs+(float)stereo*camera.eye_distance)*tv+(0.5f+(float)stereo*0.25f)*(float)camera.width; // x position on screen
	r.y =  (camera.R.yx*t.x+camera.R.yy*t.y+camera.R.yz*t.z)*rs+0.5f*(float)camera.height; // y position on screen
	rx = (int)(r.x+0.5f);
	ry = (int)(r.y+0.5f);
	rz = r.z;
	return true;
}
bool is_off_screen(const int x, const int y, const int stereo) { // check if point is off screen
	switch(stereo) {
		default: return x<                  0||x>=(int)camera.width  ||y<0||y>=(int)camera.height; // entire screen
		case -1: return x<                  0||x>=(int)camera.width/2||y<0||y>=(int)camera.height; // left half
		case +1: return x<(int)camera.width/2||x>=(int)camera.width  ||y<0||y>=(int)camera.height; // right half
	}
}
bool intersect_lines(const int x0, const int y0, const int x1, const int y1, const int xA, const int yA, const int xB, const int yB) { // check if two lines intersect
	const float d = (float)((yB-yA)*(x1-x0)-(xB-xA)*(y1-y0));
	if(d==0.0f) return false; // lines are parallel
	const float ua = ((xB-xA)*(y0-yA)-(yB-yA)*(x0-xA))/d;
	const float ub = ((x1-x0)*(y0-yA)-(y1-y0)*(x0-xA))/d;
	return ua>=0.0f && ua<=1.0f && ub>=0.0f && ub<=1.0f;
}
bool intersect_line_rectangle(const int x0, const int y0, const int x1, const int y1, const int xA, const int yA, const int xB, const int yB) { // check if line intersects rectangle
	return intersect_lines(x0, y0, x1, y1, xA, yA, xB, yA) || intersect_lines(x0, y0, x1, y1, xA, yB, xB, yB) || intersect_lines(x0, y0, x1, y1, xA, yA, xA, yB) || intersect_lines(x0, y0, x1, y1, xB, yA, xB, yB);
}
bool intersects_screen(const int x0, const int y0, const int x1, const int y1, const int stereo) {
	switch(stereo) {
		case  0: return intersect_line_rectangle(x0, y0, x1, y1,              0, 0, camera.width  , camera.height);
		case -1: return intersect_line_rectangle(x0, y0, x1, y1,              0, 0, camera.width/2, camera.height);
		case +1: return intersect_line_rectangle(x0, y0, x1, y1, camera.width/2, 0, camera.width  , camera.height);
	}
	return false;
}
Color lighting(const Color& c, const float3& p, const float3& normal, const bool translucent=false) {
	const float snb = sq(normal.x)+sq(normal.y)+sq(normal.z); // only one sqrt instead of two
	float br = 0.0f;
	for(uint i=0u; i<light_sources_n; i++) {
		const float3 d = light_sources[i]-p; // direction of light source
		const float sdb = sq(d.x)+sq(d.y)+sq(d.z);
		const float nbr = dot(d, normal)/sqrt(snb*sdb);
		br = fmax(br, translucent ? abs(nbr) : nbr);
	}
	br = fmax(0.2f, br);
	return Color((int)(br*c.r), (int)(br*c.g), (int)(br*c.b));
}

#if defined(WINDOWS_GRAPHICS)

#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#include <windows.h>
#ifndef UTILITIES_REGEX
#include <algorithm> // included in <regex> in "utilities.hpp"
#endif // UTILITIES_REGEX
HWND hWnd;
HDC frontDC, backDC;
HBITMAP frame;

void key_hold() {
	if(!camera.free) {
		zo = 0.8f*zo+0.2f*dzo; // continuous camera.zoom
		camera.zoom = exp(zo*0.25f);
	} else { // move free camera
		if(GetAsyncKeyState('W')<0) {
			camera.pos.x += camera.R.xy*camera.R.yz*fvel;
			camera.pos.y -= camera.R.xx*camera.R.yz*fvel;
			camera.pos.z -= camera.R.zz*fvel;
		}
		if(GetAsyncKeyState('A')<0) {
			camera.pos.x -= camera.R.xx*fvel;
			camera.pos.y -= camera.R.xy*fvel;
		}
		if(GetAsyncKeyState('S')<0) {
			camera.pos.x -= camera.R.xy*camera.R.yz*fvel;
			camera.pos.y += camera.R.xx*camera.R.yz*fvel;
			camera.pos.z += camera.R.zz*fvel;
		}
		if(GetAsyncKeyState('D')<0) {
			camera.pos.x += camera.R.xx*fvel;
			camera.pos.y += camera.R.xy*fvel;
		}
		if(GetAsyncKeyState(VK_SPACE)<0) {
			camera.pos.x -= camera.R.xy*camera.R.zz*fvel;
			camera.pos.y += camera.R.xx*camera.R.zz*fvel;
			camera.pos.z -= camera.R.yz*fvel;
		}
		if(GetAsyncKeyState('C')<0) {
			camera.pos.x += camera.R.xy*camera.R.zz*fvel;
			camera.pos.y -= camera.R.xx*camera.R.zz*fvel;
			camera.pos.z += camera.R.yz*fvel;
		}
	}
	if(!lockmouse) {
		if(GetAsyncKeyState('I')<0) dmy += camera.ms; // rotate camera with keys
		if(GetAsyncKeyState('J')<0) dmx += camera.ms;
		if(GetAsyncKeyState('K')<0) dmy -= camera.ms;
		if(GetAsyncKeyState('L')<0) dmx -= camera.ms;
	}
	if(autorotation) update_rotation(-1, 0);
	if(GetAsyncKeyState('Y')<0) { // adjust field of view
		camera.fov = fmin(camera.fov<1.0f ? 1.0f : camera.fov+1.0f, 179.0f);
		camera.dis = 0.5f*(float)camera.width/tan(camera.fov*pif/360.0f);
	}
	if(GetAsyncKeyState('X')<0) {
		camera.fov = fmax(camera.fov-1.0f, 1.0f);
		camera.dis = 0.5f*(float)camera.width/tan(camera.fov*pif/360.0f);
	}
	if(GetAsyncKeyState('N')<0) { // adjust camera.vr eye distance
		camera.eye_distance = fmax(camera.eye_distance-0.2f, 0.0f);
	}
	if(GetAsyncKeyState('M')<0) {
		camera.eye_distance += 0.2f;
	}
	mx = 0.8*mx+0.2*dmx; // continuous mouse movement
	my = 0.8*my+0.2*dmy;
	dmx = 0.0;
	dmy = 0.0;
	if(!lockmouse) update_rotation(mx, my);
}

void draw_bitmap(const void* buffer) {
	SetBitmapBits(frame, 4*camera.width*camera.height, buffer);
}
void draw_label(const Color& c, const string& s, const int x, const int y) {
	SetTextColor(backDC, RGB(c.r, c.g, c.b));
	TextOut(backDC, x, y, s.c_str(), (int)s.length());
	if(camera.vr) {
		if(x-camera.width/2>0) {
			TextOut(backDC, x-camera.width/2, y, s.c_str(), (int)s.length());
		} else if(x+camera.width/2<camera.width) {
			TextOut(backDC, x+camera.width/2, y, s.c_str(), (int)s.length());
		}
	}
}

void draw_pixel(const Color& c, const int x, const int y) {
	SetPixel(backDC, x, y, RGB(c.r, c.g, c.b));
}
void draw_circle(const Color& c, const int x, const int y, const int r) {
	if(r<2) {
		SetPixel(backDC, x, y, RGB(c.r, c.g, c.b));
	} else {
		SelectObject(backDC, GetStockObject(NULL_BRUSH));
		SetDCPenColor(backDC, RGB(c.r, c.g, c.b));
		if(camera.vr&&camera.tv) Ellipse(backDC, x-r/2, y-r, x+r/2, y+r);
		else Ellipse(backDC, x-r, y-r, x+r, y+r);
		SelectObject(backDC, GetStockObject(DC_BRUSH));
	}
}
void draw_line(const Color& c, const int x0, const int y0, const int x1, const int y1) {
	SetDCPenColor(backDC, RGB(c.r, c.g, c.b));
	MoveToEx(backDC, x0, y0, NULL);
	LineTo(backDC, x1, y1);
}
void draw_triangle(const Color& c, const int x0, const int y0, const int x1, const int y1, const int x2, const int y2) {
	POINT p[3];
	p[0].x = x0; p[0].y = y0;
	p[1].x = x1; p[1].y = y1;
	p[2].x = x2; p[2].y = y2;
	SetDCPenColor(backDC, RGB(c.r, c.g, c.b));
	SetDCBrushColor(backDC, RGB(c.r, c.g, c.b));
	Polygon(backDC, p, 3);
}
void draw_rectangle(const Color& c, const int x0, const int y0, const int x1, const int y1) {
	SetDCPenColor(backDC, RGB(c.r, c.g, c.b));
	SetDCBrushColor(backDC, RGB(c.r, c.g, c.b));
	Rectangle(backDC, x0, y0, x1, y1);
}
void draw_polygon(const Color& c, const int* const x, const int* const y, const int n) {
	POINT* p = new POINT[n];
	for(int i=0; i<n; i++) {
		p[i].x = x[i];
		p[i].y = y[i];
	}
	SetDCPenColor(backDC, RGB(c.r, c.g, c.b));
	SetDCBrushColor(backDC, RGB(c.r, c.g, c.b));
	Polygon(backDC, p, n);
	delete[] p;
}
void draw_text(const Color& c, const string& s, const int x, const int y) {
	SetTextColor(backDC, RGB(c.r, c.g, c.b));
	TextOut(backDC, x, y, s.c_str(), (int)s.length());
}

class Shape {
public:
	Color c = Color(0, 0, 0);
	float z = 0.0f;
	virtual void draw() const = 0;
};
class Pixel: public Shape {
private:
	int x, y;
public:
	Pixel(const Color& c, const float z, const int x, const int y) {
		this->c = c;
		this->z = z;
		this->x = x+1;
		this->y = y+1;
	}
	void draw() const override {
		draw_pixel(c, x, y);
	}
};
class Circle: public Shape {
private:
	int x, y, r;
public:
	Circle(const Color& c, const float z, const int x, const int y, const int r) {
		this->c = c;
		this->z = z;
		this->x = x;
		this->y = y;
		this->r = r;
	}
	void draw() const override {
		draw_circle(c, x, y, r);
	}
};
class Line: public Shape {
private:
	int x0, y0, x1, y1;
public:
	Line(const Color& c, float z, const int x0, const int y0, const int x1, const int y1) {
		this->c = c;
		this->z = z;
		this->x0 = x0;
		this->y0 = y0;
		this->x1 = x1;
		this->y1 = y1;
	}
	void draw() const override {
		draw_line(c, x0, y0, x1, y1);
	}
};
class Triangle: public Shape {
private:
	POINT p[3];
public:
	Triangle(const Color& c, const float z, const int x[3], const int y[3]) {
		this->c = c;
		this->z = z;
		p[0].x = x[0]; p[0].y = y[0];
		p[1].x = x[1]; p[1].y = y[1];
		p[2].x = x[2]; p[2].y = y[2];
	}
	void draw() const override {
		SetDCPenColor(backDC, RGB(c.r, c.g, c.b));
		SetDCBrushColor(backDC, RGB(c.r, c.g, c.b));
		Polygon(backDC, p, 3);
	}
};
class Quadrangle: public Shape {
private:
	POINT p[4];
public:
	Quadrangle(const Color& c, const float z, const int x[4], const int y[4]) {
		this->c = c;
		this->z = z;
		p[0].x = x[0]; p[0].y = y[0];
		p[1].x = x[1]; p[1].y = y[1];
		p[2].x = x[2]; p[2].y = y[2];
		p[3].x = x[3]; p[3].y = y[3];
	}
	void draw() const override {
		SetDCPenColor(backDC, RGB(c.r, c.g, c.b));
		SetDCBrushColor(backDC, RGB(c.r, c.g, c.b));
		Polygon(backDC, p, 4);
	}
};
class Text: public Shape {
private:
	int x, y;
	string s;
public:
	Text(const Color& c, const float z, const string& s, const int x, const int y) {
		this->c = c;
		this->z = z;
		this->s = s;
		this->x = x;
		this->y = y;
	}
	void draw() const override {
		draw_text(c, s, x, y);
	}
};

const uint maxNum = 10000u; // maximum number of 3D shapes (lines, circles, triangles) that can be drawn per frame
uint nums=0u, numr=0u;
Shape* shapes[maxNum]; // shapes for main or left screen
Shape* rights[maxNum]; // shapes for right screen

void convert_pixel(const Color& c, const float3& p, const int stereo) {
	int rx, ry; float rz;
	if(nums<maxNum&&numr<maxNum && convert(rx, ry, rz, p, stereo) && !is_off_screen(rx, ry, stereo)) {
		switch(stereo) {
			default: shapes[nums++] = new Pixel(c, rz, rx, ry); break;
			case +1: rights[numr++] = new Pixel(c, rz, rx, ry); break;
		}
	}
}
void convert_circle(const Color& c, const float3& p, const float r, const int stereo) {
	int rx, ry; float rz;
	if(nums<maxNum&&numr<maxNum && convert(rx, ry, rz, p, stereo)) {
		const float rs = camera.zoom*camera.dis/(camera.dis-rz*camera.zoom);
		const int radius = (int)(rs*r+0.5f);
		switch(stereo) {
			default: if((rx+radius>=                  0 && rx-radius<(int)camera.width   || ry+radius>=0 || ry-radius<(int)camera.height)) shapes[nums++] = new Circle(c, rz, rx, ry, radius); break; // cancel drawing if circle is off screen
			case -1: if((rx+radius>=                  0 && rx-radius<(int)camera.width/2 || ry+radius>=0 || ry-radius<(int)camera.height)) shapes[nums++] = new Circle(c, rz, rx, ry, radius); break;
			case +1: if((rx+radius>=(int)camera.width/2 && rx-radius<(int)camera.width   || ry+radius>=0 || ry-radius<(int)camera.height)) rights[numr++] = new Circle(c, rz, rx, ry, radius); break;
		}
	}
}
void convert_line(const Color& c, const float3& p0, const float3& p1, const int stereo) {
	int r0x, r0y, r1x, r1y; float r0z, r1z;
	if(nums<maxNum&&numr<maxNum && convert(r0x, r0y, r0z, p0, stereo) && convert(r1x, r1y, r1z, p1, stereo)
		&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo) && !intersects_screen(r0x, r0y, r1x, r1y, stereo))) {
		const float z = (r0z+r1z)*0.5f;
		switch(stereo) {
			default: shapes[nums++] = new Line(c, z, r0x, r0y, r1x, r1y); break;
			case +1: rights[numr++] = new Line(c, z, r0x, r0y, r1x, r1y); break;
		}
	}
}
void convert_triangle(const Color& c, const float3& p0, const float3& p1, const float3& p2, const int stereo) {
	int r0x, r0y, r1x, r1y, r2x, r2y; float r0z, r1z, r2z;
	if(nums<maxNum&&numr<maxNum && convert(r0x, r0y, r0z, p0, stereo) && convert(r1x, r1y, r1z, p1, stereo) && convert(r2x, r2y, r2z, p2, stereo)
		&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo) && is_off_screen(r2x, r2y, stereo)
		&& !intersects_screen(r0x, r0y, r1x, r1y, stereo) && !intersects_screen(r1x, r1y, r2x, r2y, stereo) && !intersects_screen(r2x, r2y, r0x, r0y, stereo))) {
		const int x[3] = { r0x, r1x, r2x };
		const int y[3] = { r0y, r1y, r2y };
		const float z = (r0z+r1z+r2z)/3.0f;
		switch(stereo) {
			default: shapes[nums++] = new Triangle(c, z, x, y); break;
			case +1: rights[numr++] = new Triangle(c, z, x, y); break;
		}
	}
}
void convert_quadrangle(const Color& c, const float3& p0, const float3& p1, const float3& p2, const float3& p3, const int stereo) {
	int r0x, r0y, r1x, r1y, r2x, r2y, r3x, r3y; float r0z, r1z, r2z, r3z;
	if(nums<maxNum&&numr<maxNum && convert(r0x, r0y, r0z, p0, stereo) && convert(r1x, r1y, r1z, p1, stereo) && convert(r2x, r2y, r2z, p2, stereo) && convert(r3x, r3y, r3z, p3, stereo)
		&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo) && is_off_screen(r2x, r2y, stereo) && is_off_screen(r3x, r3y, stereo)
		&& !intersects_screen(r0x, r0y, r1x, r1y, stereo) && !intersects_screen(r1x, r1y, r2x, r2y, stereo) && !intersects_screen(r2x, r2y, r3x, r3y, stereo) && !intersects_screen(r3x, r3y, r0x, r0y, stereo))) {
		const int x[4] = { r0x, r1x, r2x, r3x} ;
		const int y[4] = { r0y, r1y, r2y, r3y} ;
		const float z = (r0z+r1z+r2z+r3z)*0.25f;
		switch(stereo) {
			default: shapes[nums++] = new Quadrangle(c, z, x, y); break;
			case +1: rights[numr++] = new Quadrangle(c, z, x, y); break;
		}
	}
}
void convert_text(const Color& c, const float3& p, const string& s, const float r, const int stereo) {
	int rx, ry; float rz;
	if(nums<maxNum&&numr<maxNum && convert(rx, ry, rz, p, stereo)) {
		const float rs = camera.zoom*camera.dis/(camera.dis-rz*camera.zoom);
		const int radius = (int)(rs*r+0.5f);
		const float tr = fmax(0.85f*radius, 2.0f);
		switch(stereo) {
			default: if((rx+radius>=                  0 && rx-radius<(int)camera.width   || ry+radius>=0 || ry-radius<(int)camera.height)) shapes[nums++] = new Text(c, rz, s, (int)(rx+4.0f+tr), (int)(ry+3.0f+tr)); break; // cancel drawing if circle is off screen
			case -1: if((rx+radius>=                  0 && rx-radius<(int)camera.width/2 || ry+radius>=0 || ry-radius<(int)camera.height)) shapes[nums++] = new Text(c, rz, s, (int)(rx+4.0f+tr), (int)(ry+3.0f+tr)); break;
			case +1: if((rx+radius>=(int)camera.width/2 && rx-radius<(int)camera.width   || ry+radius>=0 || ry-radius<(int)camera.height)) rights[numr++] = new Text(c, rz, s, (int)(rx+4.0f+tr), (int)(ry+3.0f+tr)); break;
		}
	}
}

void draw_pixel(const Color& c, const float3& p) {
	if(!camera.vr) {
		convert_pixel(c, p,  0);
	} else {
		convert_pixel(c, p, -1);
		convert_pixel(c, p, +1);
	}
}
void draw_circle(const Color& c, const float3& p, const float r) {
	if(!camera.vr) {
		convert_circle(c, p, r,  0);
	} else {
		convert_circle(c, p, r, -1);
		convert_circle(c, p, r, +1);
	}
}
void draw_line(const Color& c, const float3& p0, const float3& p1) {
	if(!camera.vr) {
		convert_line(c, p0, p1,  0);
	} else {
		convert_line(c, p0, p1, -1);
		convert_line(c, p0, p1, +1);
	}
}
void draw_triangle(const Color& c, const float3& p0, const float3& p1, const float3& p2, const bool translucent) { // points clockwise from above
	const Color cl = lighting(c, (p0+p1+p2)/3.0f, cross(p1-p0, p2-p0), translucent);
	if(!camera.vr) {
		convert_triangle(cl, p0, p1, p2,  0);
	} else {
		convert_triangle(cl, p0, p1, p2, -1);
		convert_triangle(cl, p0, p1, p2, +1);
	}
}
void draw_quadrangle(const Color& c, const float3& p0, const float3& p1, const float3& p2, const float3& p3, const bool translucent) { // points clockwise from above, only planar points
	const Color cl = lighting(c, (p0+p1+p2+p3)*0.25f, cross(p1-p0, p2-p0), translucent);
	if(!camera.vr) {
		convert_quadrangle(cl, p0, p1, p2, p3,  0);
	} else {
		convert_quadrangle(cl, p0, p1, p2, p3, -1);
		convert_quadrangle(cl, p0, p1, p2, p3, +1);
	}
}
void draw_text(const Color& c, const float3& p, const string& s, const float r) {
	if(!camera.vr) {
		convert_text(c, p, s, r,  0);
	} else {
		convert_text(c, p, s, r, -1);
		convert_text(c, p, s, r, +1);
	}
}

void set_mouse(const int x, const int y) {
	SetCursorPos(x, y);
}
void set_clip(const int x, const int y, const int w, const int h) {
	SelectClipRgn(backDC, CreateRectRgn(x, y, x+w, y+h));
}
inline bool z_order(const Shape* const lhs, const Shape* const rhs) {
	return lhs->z<rhs->z;
}
void draw_frame() {
#ifndef SKIP_VISIBILITY_CHECKS
	if(!camera.vr) {
		std::sort(shapes, shapes+nums, z_order); // sort data array (--> visibility <--)
	} else if(numr>0u) {
		std::sort(shapes, shapes+nums, z_order);
		std::sort(rights, rights+numr, z_order);
	}
#endif // SKIP_VISIBILITY_CHECKS
	if(!camera.vr) {
		for(uint i=0u; i<nums; i++) {
			shapes[i]->draw(); // draw in right order on frame
			delete shapes[i];
			shapes[i] = nullptr;
		}
		nums = 0u;
	} else {
		set_clip(0, 0, camera.width/2, camera.height); // draw on left image only
		for(uint i=0u; i<nums; i++) {
			shapes[i]->draw();
			delete shapes[i];
			shapes[i] = nullptr;
		}
		nums = 0u;
		set_clip(camera.width/2, 0, camera.width/2, camera.height); // draw on right image only
		for(uint i=0u; i<numr; i++) {
			rights[i]->draw();
			delete rights[i];
			rights[i] = nullptr;
		}
		numr = 0u;
		set_clip(0, 0, camera.width, camera.height); // enable full drawing area again
	}
}

void update_frame(const double frametime) {
	main_label(frametime);
	BitBlt(frontDC, 0, 0, camera.width, camera.height, backDC, 0, 0, SRCCOPY); // copy back buffer to front buffer
	HPEN   oldPen   = (HPEN  )SelectObject(backDC, GetStockObject(BLACK_PEN  ));
	HBRUSH oldBrush = (HBRUSH)SelectObject(backDC, GetStockObject(BLACK_BRUSH));
	Rectangle(backDC, 0, 0, camera.width, camera.height); // clear back buffer
	SelectObject(backDC, oldPen);
	SelectObject(backDC, oldBrush);
}
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam) {
	switch(message) {
		case WM_DESTROY:
			running = false;
			PostQuitMessage(0);
			exit(0);
			return 0;
		case WM_MOUSEMOVE:
			if(!lockmouse) {
				dmx = camera.ms*(double)((int)camera.width/2 -(int)LOWORD(lParam));
				dmy = camera.ms*(double)((int)camera.height/2-(int)HIWORD(lParam));
			}
			break;
		case WM_MOUSEWHEEL:
			if(!camera.free) { // camera.zoom
				if((short)HIWORD(wParam)<0) dzo += 1.0f;
				else                        dzo -= 1.0f;
			} else if(!lockmouse) {
				if((short)HIWORD(wParam)<0) fl -= 1.0f;
				else                        fl += 1.0f;
				fvel = 0.05f*exp(fl*0.25f);
			}
			break;
		case WM_KEYDOWN:
			int key = (int)wParam;
			switch(key) {
				case VK_ESCAPE: key =  27; break; // escape
				case VK_UP    : key = -38; break; // up arrow
				case VK_DOWN  : key =  40; break; // down arrow
				case VK_LEFT  : key = -37; break; // left arrow
				case VK_RIGHT : key = -39; break; // right arrow
				case VK_PRIOR : key = -33; break; // page up
				case VK_NEXT  : key = -34; break; // page down
				case VK_HOME  : key = -36; break; // pos1
				case VK_END   : key = -35; break; // end
			}
			if(key=='U') {
				if(!lockmouse) {
					ShowCursor(true); // show cursor
				} else {
					ShowCursor(false); // hide cursor
					set_mouse(camera.width/2, camera.height/2); // reset mouse
				}
			}
			key_bindings(key);
			break;
	}
	return DefWindowProc(hWnd, message, wParam, lParam);
}
#ifdef GRAPHICS_CONSOLE
int main(int argc, char* argv[]) { // call WinMain from dummy main function in order to have an additional console window
	main_arguments = get_main_arguments(argc, argv);
	return WinMain(GetModuleHandle(NULL), NULL, GetCommandLineA(), SW_SHOWMINIMIZED);
}
#endif // GRAPHICS_CONSOLE
INT WINAPI WinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE, _In_ PSTR, _In_ INT) {
	WNDCLASS wndClass;
	wndClass.style = CS_HREDRAW|CS_VREDRAW;
	wndClass.lpfnWndProc = WndProc;
	wndClass.cbClsExtra = 0;
	wndClass.cbWndExtra = 0;
	wndClass.hInstance = hInstance;
	wndClass.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	wndClass.hCursor = LoadCursor(NULL, IDC_ARROW);
	wndClass.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
	wndClass.lpszMenuName = NULL;
	wndClass.lpszClassName = TEXT("WindowClass");
	RegisterClass(&wndClass);
	HMONITOR hMon = MonitorFromWindow(hWnd, MONITOR_DEFAULTTONEAREST);
	MONITORINFO mi = {sizeof(mi)};
	if(!GetMonitorInfo(hMon, &mi)) return 1;

	DEVMODE lpDevMode; // get monitor fps
	memset(&lpDevMode, 0, sizeof(DEVMODE));
	camera.fps_limit = (uint)EnumDisplaySettings(NULL, ENUM_CURRENT_SETTINGS, &lpDevMode)!=0 ? (uint)lpDevMode.dmDisplayFrequency : 60u; // find out screen refresh rate
	camera.width  = (uint)(mi.rcMonitor.right-mi.rcMonitor.left); // get screen size, initialize variables
	camera.height = (uint)(mi.rcMonitor.bottom-mi.rcMonitor.top);
	camera.fov = 100.0f;
	set_zoom(1.0f); // set initial zoom (weird 1.0f/1.296875f value required for backwards compatibility with previous work)
	ShowCursor(false); // hide cursor
	set_mouse(camera.width/2, camera.height/2);
	camera.update_matrix();

	hWnd = CreateWindow("WindowClass", WINDOW_NAME, WS_POPUP|WS_VISIBLE, mi.rcMonitor.left, mi.rcMonitor.top, camera.width, camera.height, NULL, NULL, hInstance, 0); // create fullscreen window
	frontDC = GetDC(hWnd);
	frame = CreateCompatibleBitmap(frontDC, camera.width, camera.height); // initialize back buffer
	backDC = CreateCompatibleDC(frontDC);
	HBITMAP oldBMP = (HBITMAP)SelectObject(backDC, frame);
	DeleteObject(oldBMP);
	SelectObject(backDC, GetStockObject(DC_PEN  ));
	SelectObject(backDC, GetStockObject(DC_BRUSH));
	//SelectObject(backDC, GetStockObject(NULL_BRUSH)); for no filling of circles and polygons
	SetDCPenColor(backDC, RGB(255, 255, 255)); // define drawing properties
	SetDCBrushColor(backDC, RGB(0, 0, 0));
	SetTextAlign(backDC, TA_TOP);
	SetBkMode(backDC, TRANSPARENT);
	SetPolyFillMode(backDC, ALTERNATE);
	HFONT hFont = CreateFont(FONT_HEIGHT+5, FONT_WIDTH+1, 0, 0, 500, 0, 0, 0, ANSI_CHARSET, 0, 0, 0, 0, "Courier New"); // (HFONT)GetStockObject(ANSI_FIXED_FONT);
	SelectObject(backDC, hFont);

	thread compute_thread(main_physics); // start main_physics() in a new thread

	MSG msg = {0};
	Clock clock;
	double frametime = 1.0;
	while(msg.message!=WM_QUIT) {
		while(PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
			if(msg.message==WM_QUIT) break;
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
		// main loop ################################################################
		key_hold();
		main_graphics();
		draw_frame();
		update_frame(frametime);
		frametime = clock.stop();
		sleep(1.0/(double)camera.fps_limit-frametime);
		clock.start();
		// ##########################################################################
	}
	ReleaseDC(hWnd, frontDC);
	compute_thread.join();
	return 0;
}

#elif defined(CONSOLE_GRAPHICS)

Image* frame = nullptr;

void draw_bitmap(const void* buffer) {
	std::copy((int*)buffer, (int*)buffer+camera.width*camera.height, frame->data());
}
void draw_label(const Color& c, const string& s, const int x, const int y) {}

uint ltw=0u, lth=0u;
uint fontwidth=8u, fontheight=16u;
void update_frame(const double frametime) {
	if(ltw==0u&&lth==0u) get_console_font_size(fontwidth, fontheight);
	uint textwidth=0u, textheight=0u;
	get_console_size(textwidth, textheight);
	textwidth = max(textwidth, 2u)-1u;
	textwidth = min(textwidth, textheight*frame->width()*fontheight/(frame->height()*fontwidth));
	textheight = min(textheight, textwidth*frame->height()*fontwidth/(frame->width()*fontheight));
	textwidth = max(textwidth, 1u);
	textheight = max(textheight, 1u);
	if(textwidth!=ltw||textheight!=lth) {
		clear_console();
		ltw = textwidth;
		lth = textheight;
	}
	show_console_cursor(false);
	print_video_dither(frame, textwidth, textheight);
	print(alignr(textwidth, to_string(textwidth)+"x"+to_string(textheight)+" "+alignr(4, to_int(1.0/frametime))+"fps"));
	show_console_cursor(true);
}

void key_detection() {
	while(running) {
		int key = key_press();
		key -= (key>96&&key<123)*32; // convert lower case to upper case
		key_bindings(key);
	}
}

void key_hold() {
	if(autorotation) update_rotation(-1, 0);
	if(!camera.free) {
		zo = 0.8f*zo+0.2f*dzo; // continuous camera.zoom
		camera.zoom = exp(zo*0.25f);
	}
	mx = 0.8*mx+0.2*dmx; // continuous mouse movement
	my = 0.8*my+0.2*dmy; // continuous mouse movement
	dmx = 0.0;
	dmy = 0.0;
	if(!lockmouse) update_rotation(mx, my);
}

int main(int argc, char* argv[]) {
	main_arguments = get_main_arguments(argc, argv);
	camera.fps_limit = 60u; // find out screen refresh rate
	camera.width  = 384u; // must be divisible by 8
	camera.height = 216u; // must be divisible by 8
	camera.fov = 100.0f;
	set_zoom(1.0f); // set initial zoom
	camera.update_matrix();

	frame = new Image(camera.width, camera.height);

	thread compute_thread(main_physics); // start main_physics() in a new thread
	thread key_thread(key_detection);

	Clock clock;
	double frametime = 1.0;
#ifdef UTILITIES_CONSOLE_DITHER_LOOKUP
	print_image_dither_initialize_lookup();
#endif // UTILITIES_CONSOLE_DITHER_LOOKUP
	clear_console();
	while(running) {
		// main loop ################################################################
		key_hold();
		main_graphics();
		update_frame(frametime);
		frametime = clock.stop();
		sleep(1.0/(double)camera.fps_limit-frametime);
		clock.start();
		// ##########################################################################
	}
	compute_thread.join();
	key_thread.join();
	return 0;
}

#else // GRAPHICS

void draw_bitmap(const void* buffer) {}
void draw_label(const Color& c, const string& s, const int x, const int y) {}

int main(int argc, char* argv[]) {
	main_arguments = get_main_arguments(argc, argv);
	camera.fps_limit = 60u; // find out screen refresh rate
	camera.width  = GRAPGICS_FRAME_WIDTH; // must be divisible by 8
	camera.height = GRAPGICS_FRAME_HEIGHT; // must be divisible by 8
	camera.fov = 100.0f;
	set_zoom(1.0f); // set initial zoom
	camera.update_matrix();

	thread compute_thread(main_physics); // start main_physics() in a new thread

	while(running) {
		// main loop ################################################################
		main_label(1.0);
		sleep(0.050);
		// ##########################################################################
	}
	compute_thread.join();
	return 0;
}

#endif // GRAPHICS
#endif // GRAPHICS