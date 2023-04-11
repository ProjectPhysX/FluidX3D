#include "kernel.hpp" // note: unbalanced round brackets () are not allowed and string literals can't be arbitrarily long, so periodically interrupt with )+R(
string opencl_c_container() { return R( // ########################## begin of OpenCL C code ####################################################################



// ################################################## utility functions ##################################################

)+R(float sq(const float x) {
	return x*x;
}
)+R(float cb(const float x) {
	return x*x*x;
}
)+R(float angle(const float3 v1, const float3 v2) {
	return acos(dot(v1, v2)/(length(v1)*length(v2)));
}
)+R(float fast_rsqrt(const float x) { // slightly fastwer approximation
	return as_float(0x5F37642F-(as_int(x)>>1));
}
)+R(float fast_asin(const float x) { // slightly fastwer approximation
	return x*fma(0.5702f, sq(sq(sq(x))), 1.0f); // 0.5707964f = (pi-2)/2
}
)+R(float fast_acos(const float x) { // slightly fastwer approximation
	return fma(fma(-0.5702f, sq(sq(sq(x))), -1.0f), x, 1.5712963f); // 0.5707964f = (pi-2)/2
}
)+R(void lu_solve(float* M, float* x, float* b, const int N, const int Nsol) { // solves system of N linear equations M*x=b within dimensionality Nsol<=N
	for(int i=0; i<Nsol; i++) { // decompose M in M=L*U
		for(int j=i+1; j<Nsol; j++) {
			M[N*j+i] /= M[N*i+i];
			for(int k=i+1; k<Nsol; k++) M[N*j+k] -= M[N*j+i]*M[N*i+k];
		}
	}
	for(int i=0; i<Nsol; i++) { // find solution of L*y=b
		x[i] = b[i];
		for(int k=0; k<i; k++) x[i] -= M[N*i+k]*x[k];
	}
	for(int i=Nsol-1; i>=0; i--) { // find solution of U*x=y
		for(int k=i+1; k<Nsol; k++) x[i] -= M[N*i+k]*x[k];
		x[i] /= M[N*i+i];
	}
}
//bool workgroup_any(const bool condition) { // returns true if any thread within the workgroup enters true
//	volatile local uint workgroup_condition; // does not work on AMD GPUs (error: non-kernel function variable cannot be declared in local address space)
//	workgroup_condition = 0u;
//	barrier(CLK_LOCAL_MEM_FENCE);
//	atomic_or(&workgroup_condition, (uint)condition);
//	barrier(CLK_LOCAL_MEM_FENCE);
//	return (bool)workgroup_condition;
//}
//bool workgroup_all(const bool condition) { // returns true if all threads within the workgroup enter true
//	volatile local uint workgroup_condition; // does not work on AMD GPUs (error: non-kernel function variable cannot be declared in local address space)
//	workgroup_condition = 1u;
//	barrier(CLK_LOCAL_MEM_FENCE);
//	atomic_and(&workgroup_condition, (uint)condition);
//	barrier(CLK_LOCAL_MEM_FENCE);
//	return (bool)workgroup_condition;
//}
)+R(void atomic_add_f(volatile global float* addr, const float val) { // not deterministic because the order of addition can vary: (a+b)+c is rounded differently than a+(b+c)
	union { // https://streamhpc.com/blog/2016-02-09/atomic-operations-for-floats-in-opencl-improved/
		uint  u32;
		float f32;
	} next, expected, current;
	current.f32 = *addr;
	do {
		next.f32 = (expected.f32=current.f32)+val; // ...*val for atomic_mul_f()
		current.u32 = atomic_cmpxchg((volatile global uint*)addr, expected.u32, next.u32);
	} while(current.u32!=expected.u32);
}
//)+"#ifdef cl_khr_int64_base_atomics"+R( // OpenCL C defines don't work in R() stringification macro
//)+R(void atomic_add_d(volatile global double* addr, const double val) { // not deterministic because the order of addition can vary: (a+b)+c is rounded differently than a+(b+c)
//	union { // https://streamhpc.com/blog/2016-02-09/atomic-operations-for-floats-in-opencl-improved/
//		ulong  u64;
//		double f64;
//	} next, expected, current;
//	current.f64 = *addr;
//	do {
//		next.f64 = (expected.f64=current.f64)+val; // ...*val for atomic_mul_d()
//		current.u64 = atom_cmpxchg((volatile global ulong*)addr, expected.u64, next.u64); // does not work on some older GPUs
//	} while(current.u64!=expected.u64);
//}
//)+"#endif"+R( // cl_khr_int64_base_atomics



// ################################################## Line3D code ##################################################

// Line3D OpenCL C version (c) Moritz Lehmann
// draw_point(...)    : draw 3D pixel
// draw_circle(...)   : draw 3D circle
// draw_line(...)     : draw 3D line
// draw_triangle(...) : draw 3D triangle
// iron_color(...)    : convert float in [0,255] to iron spectrum int color
// graphics_clear()   : kernel to reset bitmap and zbuffer

)+"#ifdef GRAPHICS"+R(
)+R(int iron_color(float x) { // coloring scheme (float 0-255 -> int color)
	x = clamp(360.0f-x*360.0f/255.0f, 0.0f, 360.0f);
	float r=255.0f, g=0.0f, b=0.0f;
	if(x<60.0f) { // white - yellow
		g = 255.0f;
		b = 255.0f-255.0f*x/60.0f;
	} else if(x<180.0f) { // yellow - red
		g = 255.0f-255.0f*(x-60.0f)/120.0f;
	} else if(x<270.0f) { // red - violet
		r = 255.0f-255.0f*(x-180.0f)/180.0f;
		b = 255.0f*(x-180.0f)/90.0f;
	} else { // violet - black
		r = 255.0f-255.0f*(x-180.0f)/180.0f;
		b = 255.0f-255.0f*(x-270.0f)/90.0f;
	}
	return (((int)r)<<16)|(((int)g)<<8)|((int)b);
}
)+R(int rainbow_color(float x) { // coloring scheme (float 0-255 -> int color)
	x = clamp(360.0f-x*360.0f/255.0f, 0.0f, 360.0f);
	float r=0.0f, g=0.0f, b=0.0f; // black
	if(x<60.0f) { // red - yellow
		r = 255.0f;
		g = 255.0f*x/60.0f;
	} else if(x>=60.0f&&x<120.0f) { // yellow - green
		r = 255.0f-255.0f*(x-60.0f)/60.0f;
		g = 255.0f;
	} else if(x>=120.0f&&x<180.0f) { // green - cyan
		g = 255.0f;
		b = 255.0f*(x-120.0f)/60.0f;
	} else if(x>=180.0f&&x<240.0f) { // cyan - blue
		g = 255.0f-255.0f*(x-180.0f)/60.0f;
		b = 255.0f;
	} else if(x>=240.0f&&x<300.0f) { // blue - violet
		r = (255.0f*(x-240.0f)/60.0f)/2.0f;
		b = 255.0f;
	} else { // violet - black
		r = (255.0f-255.0f*(x-300.0f)/60.0f)/2.0f;
		b = 255.0f-255.0f*(x-300.0f)/60.0f;
	}
	return (((int)r)<<16)|(((int)g)<<8)|((int)b);
}
)+R(int color_dim(const int c, const float x) {
	const int r = clamp((int)fma((float)((c>>16)&255), x, 0.5f), 0, 255);
	const int g = clamp((int)fma((float)((c>> 8)&255), x, 0.5f), 0, 255);
	const int b = clamp((int)fma((float)( c     &255), x, 0.5f), 0, 255);
	return (r&255)<<16|(g&255)<<8|(b&255);
}
)+R(int color_mix(const int c1, const int c2, const float w) {
	const uchar4 cc1=as_uchar4(c1), cc2=as_uchar4(c2);
	const float3 fc1=(float3)((float)cc1.x, (float)cc1.y, (float)cc1.z), fc2=(float3)((float)cc2.x, (float)cc2.y, (float)cc2.z);
	const float3 fcm = fma(w, fc1, fma(1.0f-w, fc2, (float3)(0.5f, 0.5f, 0.5f)));
	return as_int((uchar4)((uchar)fcm.x, (uchar)fcm.y, (uchar)fcm.z, (uchar)0u));
}
)+R(int color_mix_3(const int c0, const int c1, const int c2, const float w0, const float w1, const float w2) { // w0+w1+w2 = 1
	const uchar4 cc0=as_uchar4(c0), cc1=as_uchar4(c1), cc2=as_uchar4(c2);
	const float3 fc0=(float3)((float)cc0.x, (float)cc0.y, (float)cc0.z),  fc1=(float3)((float)cc1.x, (float)cc1.y, (float)cc1.z), fc2=(float3)((float)cc2.x, (float)cc2.y, (float)cc2.z);
	const float3 fcm = fma(w0, fc0, fma(w1, fc1, fma(w2, fc2, (float3)(0.5f, 0.5f, 0.5f))));
	return as_int((uchar4)((uchar)fcm.x, (uchar)fcm.y, (uchar)fcm.z, (uchar)0u));
}
)+R(int hsv_to_rgb(const float h, const float s, const float v) {
	const float c = v*s;
	const float x = c*(1.0f-fabs(fmod(h/60.0f, 2.0f)-1.0f));
	const float m = v-c;
	float r=0.0f, g=0.0f, b=0.0f;
	if(0.0f<=h&&h<60.0f) { r = c; g = x; }
	else if(h<120.0f) { r = x; g = c; }
	else if(h<180.0f) { g = c; b = x; }
	else if(h<240.0f) { g = x; b = c; }
	else if(h<300.0f) { r = x; b = c; }
	else if(h<360.0f) { r = c; b = x; }
	return (int)((r+m)*255.0f)<<16|(int)((g+m)*255.0f)<<8|(int)((b+m)*255.0f);
}
)+R(int lighting(const int c, const float3 p, const float3 normal, const float* camera_cache) { // calculate lighting of triangle
	const float dis  = camera_cache[ 1]; // fetch camera parameters (rotation matrix, camera position, etc.)
	const float posx = camera_cache[ 2]-def_domain_offset_x;
	const float posy = camera_cache[ 3]-def_domain_offset_y;
	const float posz = camera_cache[ 4]-def_domain_offset_z;
	const float Rzx  = camera_cache[11];
	const float Rzy  = camera_cache[12];
	const float Rzz  = camera_cache[13];
	const uchar cr=c>>16&255, cg=c>>8&255, cb=c&255;
	const float nl2 = sq(normal.x)+sq(normal.y)+sq(normal.z); // only one native_sqrt instead of two
	const float dx = p.x-fma(Rzx, dis, posx); // direction of light source is viewing direction
	const float dy = p.y-fma(Rzy, dis, posy);
	const float dz = p.z-fma(Rzz, dis, posz);
	const float dl2 = sq(dx)+sq(dy)+sq(dz);
	const float br = max(1.5f*fabs(normal.x*dx+normal.y*dy+normal.z*dz)*rsqrt(nl2*dl2), 0.3f);
	return min((int)(br*cr), 255)<<16|min((int)(br*cg), 255)<<8|min((int)(br*cb), 255);
}
)+R(bool is_off_screen(const int x, const int y, const int stereo) {
	switch(stereo) {
		default: return x<                 0||x>=def_screen_width  ||y<0||y>=def_screen_height; // entire screen
		case -1: return x<                 0||x>=def_screen_width/2||y<0||y>=def_screen_height; // left half
		case +1: return x<def_screen_width/2||x>=def_screen_width  ||y<0||y>=def_screen_height; // right half
	}
}
)+R(void draw(const int x, const int y, const float z, const int color, global int* bitmap, volatile global int* zbuffer, const int stereo) {
	const int index=x+y*def_screen_width, iz=(int)(z*(2147483647.0f/10000.0f)); // use int z-buffer and atomic_max to minimize noise in image
	if(!is_off_screen(x, y, stereo)&&iz>atomic_max(&zbuffer[index], iz)) bitmap[index] = color; // only draw if point is on screen and first in zbuffer
}
)+R(bool convert(int* rx, int* ry, float* rz, const float3 p, const float* camera_cache, const int stereo) { // 3D -> 2D
	const float zoom = camera_cache[0]; // fetch camera parameters (rotation matrix, camera position, etc.)
	const float dis  = camera_cache[1];
	const float posx = camera_cache[2];
	const float posy = camera_cache[3];
	const float posz = camera_cache[4];
	const float Rxx  = camera_cache[5];
	const float Rxy  = camera_cache[6];
	const float Rxz  = camera_cache[7];
	const float Ryx  = camera_cache[8];
	const float Ryy  = camera_cache[9];
	const float Ryz  = camera_cache[10];
	const float Rzx  = camera_cache[11];
	const float Rzy  = camera_cache[12];
	const float Rzz  = camera_cache[13];
	const float eye_distance = vload_half(28, (half*)camera_cache);
	float3 t, r;
	t.x = p.x+def_domain_offset_x-posx-(float)stereo*eye_distance/zoom*Rxx; // transformation
	t.y = p.y+def_domain_offset_y-posy-(float)stereo*eye_distance/zoom*Rxy;
	t.z = p.z+def_domain_offset_z-posz;
	r.z = Rzx*t.x+Rzy*t.y+Rzz*t.z; // z-position for z-buffer
	const float rs = zoom*dis/(dis-r.z*zoom); // perspective (reciprocal is more efficient)
	if(rs<=0.0f) return false; // point is behins camera
	const float tv = ((as_int(camera_cache[14])>>30)&0x1)&&stereo!=0 ? 0.5f : 1.0f;
	r.x = ((Rxx*t.x+Rxy*t.y+Rxz*t.z)*rs+(float)stereo*eye_distance)*tv+(0.5f+(float)stereo*0.25f)*(float)def_screen_width; // x position on screen
	r.y =  (Ryx*t.x+Ryy*t.y+Ryz*t.z)*rs+0.5f*(float)def_screen_height; // y position on screen
	*rx = (int)(r.x+0.5f);
	*ry = (int)(r.y+0.5f);
	*rz = r.z;
	return true;
}
)+R(void convert_circle(float3 p, const float r, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer, const int stereo) { // 3D -> 2D
	int rx, ry; float rz;
	if(convert(&rx, &ry, &rz, p, camera_cache, stereo)) {
		const float zoom = camera_cache[0];
		const float dis  = camera_cache[1];
		const float rs = zoom*dis/(dis-rz*zoom);
		const int radius = (int)(rs*r+0.5f);
		switch(stereo) {
			default: if(rx<                       -radius||rx>=(int)def_screen_width  +radius || ry<-radius||ry>=(int)def_screen_height+radius) return; break; // cancel drawing if circle is off screen
			case -1: if(rx<                       -radius||rx>=(int)def_screen_width/2+radius || ry<-radius||ry>=(int)def_screen_height+radius) return; break;
			case +1: if(rx<(int)def_screen_width/2-radius||rx>=(int)def_screen_width  +radius || ry<-radius||ry>=(int)def_screen_height+radius) return; break;
		}
		int d=-radius, x=radius, y=0; // Bresenham algorithm for circle
		while(x>=y) {
			draw(rx+x, ry+y, rz, color, bitmap, zbuffer, stereo);
			draw(rx-x, ry+y, rz, color, bitmap, zbuffer, stereo);
			draw(rx+x, ry-y, rz, color, bitmap, zbuffer, stereo);
			draw(rx-x, ry-y, rz, color, bitmap, zbuffer, stereo);
			draw(rx+y, ry+x, rz, color, bitmap, zbuffer, stereo);
			draw(rx-y, ry+x, rz, color, bitmap, zbuffer, stereo);
			draw(rx+y, ry-x, rz, color, bitmap, zbuffer, stereo);
			draw(rx-y, ry-x, rz, color, bitmap, zbuffer, stereo);
			d += 2*y+1;
			y++;
			if(d>0) d-=2*(--x);
		}
	}
}
)+R(void convert_line(const float3 p0, const float3 p1, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer, const int stereo) { // 3D -> 2D
	int r0x, r0y, r1x, r1y; float r0z, r1z;
	if(convert(&r0x, &r0y, &r0z, p0, camera_cache, stereo) && convert(&r1x, &r1y, &r1z, p1, camera_cache, stereo)
		&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo))) { // cancel drawing if both points are off screen
		int x=r0x, y=r0y; // Bresenham algorithm
		const float z = 0.5f*(r0z+r1z); // approximate line z position for each pixel to be equal
		const int dx= abs(r1x-r0x), sx=2*(r0x<r1x)-1;
		const int dy=-abs(r1y-r0y), sy=2*(r0y<r1y)-1;
		int err = dx+dy;
		while(x!=r1x||y!=r1y) {
			draw(x, y, z, color, bitmap, zbuffer, stereo);
			const int e2 = 2*err;
			if(e2>dy) { err+=dy; x+=sx; }
			if(e2<dx) { err+=dx; y+=sy; }
		}
	}
}
)+R(void convert_triangle(float3 p0, float3 p1, float3 p2, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer, const int stereo) { // 3D -> 2D
	int r0x, r0y, r1x, r1y, r2x, r2y; float r0z, r1z, r2z;
	if(convert(&r0x, &r0y, &r0z, p0, camera_cache, stereo) && convert(&r1x, &r1y, &r1z, p1, camera_cache, stereo) && convert(&r2x, &r2y, &r2z, p2, camera_cache, stereo)
		&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo) && is_off_screen(r2x, r2y, stereo))) { // cancel drawing if all points are off screen
		if(r0x*(r1y-r2y)+r1x*(r2y-r0y)+r2x*(r0y-r1y)>100000 || (r0y==r1y&&r0y==r2y)) return; // return for large triangle area or degenerate triangles
		if(r0y>r1y) { const int xt = r0x; const int yt = r0y; r0x = r1x; r0y = r1y; r1x = xt; r1y = yt; } // sort vertices ascending by y
		if(r0y>r2y) { const int xt = r0x; const int yt = r0y; r0x = r2x; r0y = r2y; r2x = xt; r2y = yt; }
		if(r1y>r2y) { const int xt = r1x; const int yt = r1y; r1x = r2x; r1y = r2y; r2x = xt; r2y = yt; }
		const float z = (r0z+r1z+r2z)/3.0f; // approximate triangle z position for each pixel to be equal
		for(int y=r0y; y<r1y; y++) { // Bresenham algorithm (lower triangle half)
			const int xA = r0x+(r2x-r0x)*(y-r0y)/(r2y-r0y);
			const int xB = r0x+(r1x-r0x)*(y-r0y)/(r1y-r0y);
			for(int x=min(xA, xB); x<max(xA, xB); x++) {
				draw(x, y, z, color, bitmap, zbuffer, stereo);
			}
		}
		for(int y=r1y; y<r2y; y++) { // Bresenham algorithm (upper triangle half)
			const int xA = r0x+(r2x-r0x)*(y-r0y)/(r2y-r0y);
			const int xB = r1x+(r2x-r1x)*(y-r1y)/(r2y-r1y);
			for(int x=min(xA, xB); x<max(xA, xB); x++) {
				draw(x, y, z, color, bitmap, zbuffer, stereo);
			}
		}
	}
}
)+R(void convert_triangle_interpolated(float3 p0, float3 p1, float3 p2, int c0, int c1, int c2, const float* camera_cache, global int* bitmap, global int* zbuffer, const int stereo) { // 3D -> 2D
	int r0x, r0y, r1x, r1y, r2x, r2y; float r0z, r1z, r2z;
	if(convert(&r0x, &r0y, &r0z, p0, camera_cache, stereo) && convert(&r1x, &r1y, &r1z, p1, camera_cache, stereo) && convert(&r2x, &r2y, &r2z, p2, camera_cache, stereo)
		&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo) && is_off_screen(r2x, r2y, stereo))) { // cancel drawing if all points are off screen
		if(r0x*(r1y-r2y)+r1x*(r2y-r0y)+r2x*(r0y-r1y)>100000 || (r0y==r1y&&r0y==r2y)) return; // return for large triangle area or degenerate triangles
		if(r0y>r1y) { const int xt = r0x; const int yt = r0y; r0x = r1x; r0y = r1y; r1x = xt; r1y = yt; const int ct = c0; c0 = c1; c1 = ct; } // sort vertices ascending by y
		if(r0y>r2y) { const int xt = r0x; const int yt = r0y; r0x = r2x; r0y = r2y; r2x = xt; r2y = yt; const int ct = c0; c0 = c2; c2 = ct; }
		if(r1y>r2y) { const int xt = r1x; const int yt = r1y; r1x = r2x; r1y = r2y; r2x = xt; r2y = yt; const int ct = c1; c1 = c2; c2 = ct; }
		const float z = (r0z+r1z+r2z)/3.0f; // approximate triangle z position for each pixel to be equal
		const float d = (float)((r1y-r2y)*(r0x-r2x)+(r2x-r1x)*(r0y-r2y));
		for(int y=r0y; y<r1y; y++) { // Bresenham algorithm (lower triangle half)
			const int xA = r0x+(r2x-r0x)*(y-r0y)/(r2y-r0y);
			const int xB = r0x+(r1x-r0x)*(y-r0y)/(r1y-r0y);
			for(int x=min(xA, xB); x<max(xA, xB); x++) {
				const float w0 = (float)((r1y-r2y)*(x-r2x)+(r2x-r1x)*(y-r2y))/d; // barycentric coordinates
				const float w1 = (float)((r2y-r0y)*(x-r2x)+(r0x-r2x)*(y-r2y))/d;
				const float w2 = 1.0f-w0-w1;
				const int color = color_mix_3(c0, c1, c2, w0, w1, w2); // interpolate color
				draw(x, y, z, color, bitmap, zbuffer, stereo);
			}
		}
		for(int y=r1y; y<r2y; y++) { // Bresenham algorithm (upper triangle half)
			const int xA = r0x+(r2x-r0x)*(y-r0y)/(r2y-r0y);
			const int xB = r1x+(r2x-r1x)*(y-r1y)/(r2y-r1y);
			for(int x=min(xA, xB); x<max(xA, xB); x++) {
				const float w0 = (float)((r1y-r2y)*(x-r2x)+(r2x-r1x)*(y-r2y))/d; // barycentric coordinates
				const float w1 = (float)((r2y-r0y)*(x-r2x)+(r0x-r2x)*(y-r2y))/d;
				const float w2 = 1.0f-w0-w1;
				const int color = color_mix_3(c0, c1, c2, w0, w1, w2); // interpolate color
				draw(x, y, z, color, bitmap, zbuffer, stereo);
			}
		}
	}
}
)+R(void draw_point(const float3 p, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer) { // 3D -> 2D
	const bool vr = (as_int(camera_cache[14])>>31)&0x1;
	int rx, ry; float rz;
	if(!vr) {
		if(convert(&rx, &ry, &rz, p, camera_cache,  0)) draw(rx, ry, rz, color, bitmap, zbuffer,  0);
	} else {
		if(convert(&rx, &ry, &rz, p, camera_cache, -1)) draw(rx, ry, rz, color, bitmap, zbuffer, -1); // left eye
		if(convert(&rx, &ry, &rz, p, camera_cache, +1)) draw(rx, ry, rz, color, bitmap, zbuffer, +1); // right eye
	}
}
)+R(void draw_circle(const float3 p, const float r, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer) { // 3D -> 2D
	const bool vr = (as_int(camera_cache[14])>>31)&0x1;
	if(!vr) {
		convert_circle(p, r, color, camera_cache, bitmap, zbuffer,  0);
	} else {
		convert_circle(p, r, color, camera_cache, bitmap, zbuffer, -1); // left eye
		convert_circle(p, r, color, camera_cache, bitmap, zbuffer, +1); // right eye
	}
}
)+R(void draw_line(const float3 p0, const float3 p1, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer) { // 3D -> 2D
	const bool vr = (as_int(camera_cache[14])>>31)&0x1;
	if(!vr) {
		convert_line(p0, p1, color, camera_cache, bitmap, zbuffer,  0);
	} else {
		convert_line(p0, p1, color, camera_cache, bitmap, zbuffer, -1); // left eye
		convert_line(p0, p1, color, camera_cache, bitmap, zbuffer, +1); // right eye
	}
}
)+R(void draw_triangle(const float3 p0, const float3 p1, const float3 p2, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer) { // 3D -> 2D
	const bool vr = (as_int(camera_cache[14])>>31)&0x1;
	if(!vr) {
		convert_triangle(p0, p1, p2, color, camera_cache, bitmap, zbuffer,  0);
	} else {
		convert_triangle(p0, p1, p2, color, camera_cache, bitmap, zbuffer, -1); // left eye
		convert_triangle(p0, p1, p2, color, camera_cache, bitmap, zbuffer, +1); // right eye
	}
}
)+R(void draw_triangle_interpolated(const float3 p0, const float3 p1, const float3 p2, const int c0, const int c1, const int c2, const float* camera_cache, global int* bitmap, global int* zbuffer) { // 3D -> 2D
	const bool vr = (as_int(camera_cache[14])>>31)&0x1;
	if(!vr) {
		convert_triangle_interpolated(p0, p1, p2, c0, c1, c2, camera_cache, bitmap, zbuffer,  0);
	} else {
		convert_triangle_interpolated(p0, p1, p2, c0, c1, c2, camera_cache, bitmap, zbuffer, -1); // left eye
		convert_triangle_interpolated(p0, p1, p2, c0, c1, c2, camera_cache, bitmap, zbuffer, +1); // right eye
	}
}
)+R(kernel void graphics_clear(global int* bitmap, global int* zbuffer) {
	const uint n = get_global_id(0);
	bitmap[n] = def_background_color; // black background = 0x000000, use 0xFFFFFF for white background
	zbuffer[n] = -2147483648;
}
)+R(constant ushort edge_table_data[128] = { // source: Paul Bourke, http://paulbourke.net/geometry/polygonise/, mirror symmetry applied, makes marching-cubes 31% faster
	0x000, 0x109, 0x203, 0x30A, 0x406, 0x50F, 0x605, 0x70C, 0x80C, 0x905, 0xA0F, 0xB06, 0xC0A, 0xD03, 0xE09, 0xF00,
	0x190, 0x099, 0x393, 0x29A, 0x596, 0x49F, 0x795, 0x69C, 0x99C, 0x895, 0xB9F, 0xA96, 0xD9A, 0xC93, 0xF99, 0xE90,
	0x230, 0x339, 0x033, 0x13A, 0x636, 0x73F, 0x435, 0x53C, 0xA3C, 0xB35, 0x83F, 0x936, 0xE3A, 0xF33, 0xC39, 0xD30,
	0x3A0, 0x2A9, 0x1A3, 0x0AA, 0x7A6, 0x6AF, 0x5A5, 0x4AC, 0xBAC, 0xAA5, 0x9AF, 0x8A6, 0xFAA, 0xEA3, 0xDA9, 0xCA0,
	0x460, 0x569, 0x663, 0x76A, 0x066, 0x16F, 0x265, 0x36C, 0xC6C, 0xD65, 0xE6F, 0xF66, 0x86A, 0x963, 0xA69, 0xB60,
	0x5F0, 0x4F9, 0x7F3, 0x6FA, 0x1F6, 0x0FF, 0x3F5, 0x2FC, 0xDFC, 0xCF5, 0xFFF, 0xEF6, 0x9FA, 0x8F3, 0xBF9, 0xAF0,
	0x650, 0x759, 0x453, 0x55A, 0x256, 0x35F, 0x055, 0x15C, 0xE5C, 0xF55, 0xC5F, 0xD56, 0xA5A, 0xB53, 0x859, 0x950,
	0x7C0, 0x6C9, 0x5C3, 0x4CA, 0x3C6, 0x2CF, 0x1C5, 0x0CC, 0xFCC, 0xEC5, 0xDCF, 0xCC6, 0xBCA, 0xAC3, 0x9C9, 0x8C0
};
)+R(constant uchar triangle_table_data[1920] = { // source: Paul Bourke, http://paulbourke.net/geometry/polygonise/, termination value 15, bit packed
	255,255,255,255,255,255,255, 15, 56,255,255,255,255,255,255, 16,249,255,255,255,255,255, 31, 56,137,241,255,255,255,255, 33,250,255,255,255,255,255, 15, 56, 33,250,255,255,255,255, 41, 10,146,
	255,255,255,255, 47, 56,162,168,137,255,255,255,179,242,255,255,255,255,255, 15, 43,184,240,255,255,255,255,145, 32,179,255,255,255,255, 31, 43,145,155,184,255,255,255,163,177, 58,255,255,255,
	255, 15, 26,128,138,171,255,255,255,147, 48,155,171,249,255,255,159,168,138,251,255,255,255,255,116,248,255,255,255,255,255, 79,  3, 55,244,255,255,255,255, 16,137,116,255,255,255,255, 79,145,
	116,113, 19,255,255,255, 33,138,116,255,255,255,255, 63,116,  3, 20,162,255,255,255, 41,154, 32, 72,247,255,255, 47,154,146, 39, 55,151,244,255, 72, 55, 43,255,255,255,255,191,116, 43, 36, 64,
	255,255,255,  9,129,116, 50,251,255,255, 79,183, 73,155, 43, 41,241,255,163, 49,171,135,244,255,255, 31,171, 65, 27, 64,183,244,255,116,152,176,185,186, 48,255, 79,183,180,153,171,255,255,255,
	 89,244,255,255,255,255,255,159, 69,128,243,255,255,255,255, 80, 20,  5,255,255,255,255,143, 69, 56, 53, 81,255,255,255, 33,154, 69,255,255,255,255, 63,128, 33, 74, 89,255,255,255, 37, 90, 36,
	  4,242,255,255, 47, 90, 35, 53, 69, 67,248,255, 89, 36,179,255,255,255,255, 15, 43,128, 75, 89,255,255,255, 80,  4, 81, 50,251,255,255, 47, 81, 82, 40,184,132,245,255, 58,171, 49, 89,244,255,
	255, 79, 89,128,129, 26,184,250,255, 69, 80,176,181,186, 48,255, 95,132,133,170,184,255,255,255,121, 88,151,255,255,255,255,159,  3, 89, 83, 55,255,255,255,112,  8,113, 81,247,255,255, 31, 53,
	 83,247,255,255,255,255,121,152,117, 26,242,255,255,175, 33, 89, 80,  3,117,243,255,  8,130, 82, 88,167, 37,255, 47, 90, 82, 51,117,255,255,255,151,117,152,179,242,255,255,159,117,121,146,  2,
	114,251,255, 50, 11,129,113, 24,117,255,191, 18, 27,119, 81,255,255,255, 89,136,117, 26,163,179,255, 95,  7,  5,121, 11,  1,186, 10,171,176, 48, 90,128,112,117,176, 90,183,245,255,255,255,255,
	106,245,255,255,255,255,255, 15, 56,165,246,255,255,255,255,  9, 81,106,255,255,255,255, 31, 56,145, 88,106,255,255,255, 97, 37, 22,255,255,255,255, 31, 86, 33, 54,128,255,255,255,105,149, 96,
	 32,246,255,255, 95,137,133, 82, 98, 35,248,255, 50,171, 86,255,255,255,255,191,128, 43,160, 86,255,255,255, 16, 41,179,165,246,255,255, 95,106,145,146, 43,137,251,255, 54,107, 53, 21,243,255,
	255, 15,184,176,  5, 21,181,246,255,179,  6, 99, 96,  5,149,255,111,149,150,187,137,255,255,255,165, 70,135,255,255,255,255, 79,  3,116, 99,165,255,255,255,145, 80,106, 72,247,255,255,175, 86,
	145, 23, 55,151,244,255, 22, 98, 21,116,248,255,255, 31, 82, 37, 54, 64, 67,247,255, 72,151, 80, 96,  5, 98,255,127,147,151, 52,146,149, 38,150,179,114, 72,106,245,255,255, 95,106,116, 66,  2,
	114,251,255, 16, 73,135, 50, 91,106,255,159, 18,185,146,180,183, 84,106, 72, 55, 91, 83, 81,107,255, 95,177,181, 22,176,183,  4,180, 80,  9, 86, 48,182, 54, 72,103,149,150, 75,151,183,249,255,
	 74,105,164,255,255,255,255, 79,106,148, 10, 56,255,255,255, 10,161,  6, 70,240,255,255,143, 19, 24,134, 70, 22,250,255, 65, 25, 66, 98,244,255,255, 63,128, 33, 41,148, 98,244,255, 32, 68, 98,
	255,255,255,255,143, 35, 40, 68, 98,255,255,255, 74,169, 70, 43,243,255,255, 15, 40,130, 75,169,164,246,255,179,  2, 97, 96,100,161,255,111, 20, 22, 74, 24, 18,139, 27,105,148, 99, 25,179, 54,
	255,143, 27, 24,176, 22, 25,100, 20,179, 54,  6, 96,244,255,255,111,132,107,248,255,255,255,255,167,118,168,152,250,255,255, 15, 55,160,  7,169,118,250,255,106, 23,122,113, 24,  8,255,175,118,
	122, 17, 55,255,255,255, 33, 22,134,129,137,118,255, 47,150,146, 97,151,144,115,147,135,112, 96,  6,242,255,255,127, 35,118,242,255,255,255,255, 50,171,134,138,137,118,255, 47,112,114, 11,121,
	118,154,122,129, 16,135,161,103,167, 50,187, 18, 27,167, 22,118,241,255,152,134,118, 25,182, 54, 49,  6, 25,107,247,255,255,255,255,135,112, 96,179,176,  6,255,127,107,255,255,255,255,255,255,
	103,251,255,255,255,255,255, 63,128,123,246,255,255,255,255, 16,185,103,255,255,255,255,143,145, 56,177,103,255,255,255, 26, 98,123,255,255,255,255, 31,162,  3,104,123,255,255,255,146, 32,154,
	182,247,255,255,111,123,162,163, 56,154,248,255, 39, 99,114,255,255,255,255,127,128,103, 96,  2,255,255,255,114, 38,115, 16,249,255,255, 31, 38,129, 22,137,120,246,255,122,166,113, 49,247,255,
	255,175,103,113, 26,120,  1,248,255, 48,  7,167,160,105,122,255,127,166,167,136,154,255,255,255,134,180,104,255,255,255,255, 63,182,  3,  6,100,255,255,255,104,139,100,  9,241,255,255,159,100,
	105,147, 19, 59,246,255,134,100,139,162,241,255,255, 31,162,  3, 11,182, 64,246,255,180, 72,182, 32, 41,154,255,175, 57, 58,146, 52, 59, 70, 54, 40,131, 36,100,242,255,255, 15, 36,100,242,255,
	255,255,255,145, 32, 67, 66, 70,131,255, 31, 73, 65, 34,100,255,255,255, 24,131, 22, 72,102, 26,255,175,  1, 10,102, 64,255,255,255,100, 67,131,166,  3,147,154,163, 73,166,244,255,255,255,255,
	148,117,182,255,255,255,255, 15, 56,148,181,103,255,255,255,  5, 81,  4,103,251,255,255,191,103, 56, 52, 69, 19,245,255, 89,164, 33,103,251,255,255,111,123, 33, 10, 56,148,245,255,103, 91,164,
	 36, 74, 32,255, 63,132, 83, 52, 82, 90,178,103, 39,115, 38, 69,249,255,255,159, 69,128,  6, 38,134,247,255, 99, 50,103, 81, 80,  4,255,111,130,134, 39,129,132, 21,133, 89,164, 97,113, 22,115,
	255, 31,166,113, 22,112,120,144, 69,  4, 74, 90, 48,106,122,115,122,166,167, 88,164,132,250,255,150,101,155,139,249,255,255, 63,182, 96,  3,101,144,245,255,176,  8,181, 16, 85,182,255,111, 59,
	 54, 85, 19,255,255,255, 33,154,181,185,184,101,255, 15, 59, 96, 11,105,101, 25,162,139,181,101,  8,165, 37, 32,101, 59, 54, 37, 58, 90,243,255,133, 89,130,101, 50, 40,255,159,101,105,  0, 38,
	255,255,255, 81, 24,  8,101, 56, 40, 38, 24,101, 18,246,255,255,255,255, 49, 22,166,131, 86,150,152,166,  1, 10,150,  5,101,240,255, 48, 88,166,255,255,255,255,175,101,255,255,255,255,255,255,
	 91,122,181,255,255,255,255,191,165,123,133,  3,255,255,255,181, 87,186,145,240,255,255,175, 87,186,151, 24, 56,241,255, 27,178, 23, 87,241,255,255, 15, 56, 33, 23, 87, 39,251,255,121,149,114,
	  9, 34,123,255,127, 37, 39, 91, 41, 35,152, 40, 82, 42, 83,115,245,255,255,143,  2, 88,130, 87, 42,245,255,  9, 81, 58, 53, 55, 42,255,159, 40, 41,129, 39, 42,117, 37, 49, 53, 87,255,255,255,
	255, 15,120,112, 17, 87,255,255,255,  9,147, 83, 53,247,255,255,159,120,149,247,255,255,255,255,133, 84,138,186,248,255,255, 95, 64,181, 80,186, 59,240,255, 16,137,164,168,171, 84,255,175, 75,
	 74,181, 67, 73, 49, 65, 82, 33, 88,178, 72,133,255, 15,180,176, 67,181,178, 81,177, 32,  5,149,178, 69,133,139,149, 84,178,243,255,255,255,255, 82, 58, 37, 67, 53, 72,255, 95, 42, 37, 68,  2,
	255,255,255,163, 50,165,131, 69,133, 16, 89, 42, 37, 20, 41, 73,242,255, 72,133, 53, 83,241,255,255, 15, 84,  1,245,255,255,255,255, 72,133, 53,  9,  5, 83,255,159, 84,255,255,255,255,255,255,
	180, 71,185,169,251,255,255, 15, 56,148,151,123,169,251,255,161, 27, 75, 65,112,180,255, 63, 65, 67, 24, 74, 71,171, 75,180,151, 75, 41,155, 33,255,159, 71,185,151,177,178,  1, 56,123,180, 36,
	 66,240,255,255,191, 71, 75,130, 67, 35,244,255,146, 42,151, 50,119,148,255,159,122,121,164,114,120, 32,112,115, 58, 42, 71, 26, 10,  4, 26, 42,120,244,255,255,255,255,148, 65,113, 23,243,255,
	255, 79, 25, 20,  7, 24,120,241,255,  4,115, 52,255,255,255,255, 79,120,255,255,255,255,255,255,169,168,139,255,255,255,255, 63,144,147,187,169,255,255,255, 16, 10,138,168,251,255,255, 63,161,
	 59,250,255,255,255,255, 33, 27,155,185,248,255,255, 63,144,147, 27,146,178,249,255, 32,139,176,255,255,255,255, 63,178,255,255,255,255,255,255, 50, 40,168,138,249,255,255,159, 42,144,242,255,
	255,255,255, 50, 40,168, 16, 24,138,255, 31, 42,255,255,255,255,255,255, 49,152,129,255,255,255,255, 15, 25,255,255,255,255,255,255, 48,248,255,255,255,255,255,255,255,255,255,255,255,255,255
};
)+R(ushort edge_table(const uint i) {
	return edge_table_data[i<128u?i:255u-i];
}
)+R(uchar triangle_table(const uint i) {
	return (triangle_table_data[i/2u]>>(4u*(i%2u)))&0xF;
}
)+R(float3 interpolate_vertex(const float3 p1, const float3 p2, const float v1, const float v2, const float iso) { // linearly interpolate position where isosurface cuts an edge between 2 vertices
	const float w = (iso-v1)/(v2-v1);
	return (1.0f-w)*p1+w*p2;
}
)+R(uint marching_cubes(const float* v, const float iso, float3* triangles) { // input: 8 values v, isovalue; output: returns number of triangles, 15 triangle vertices t
	uint cube = 0u; // determine index of which vertices are inside of the isosurface
	for(uint i=0u; i<8u; i++) cube |= (v[i]<iso)<<i;
	if(cube==0u || cube==255u) return 0u; // cube is entirely inside/outside of the isosurface
	float3 p[8]; // definition of unit cube corners
	p[0] = (float3)(0.0f, 0.0f, 0.0f);
	p[1] = (float3)(1.0f, 0.0f, 0.0f);
	p[2] = (float3)(1.0f, 0.0f, 1.0f);
	p[3] = (float3)(0.0f, 0.0f, 1.0f);
	p[4] = (float3)(0.0f, 1.0f, 0.0f);
	p[5] = (float3)(1.0f, 1.0f, 0.0f);
	p[6] = (float3)(1.0f, 1.0f, 1.0f);
	p[7] = (float3)(0.0f, 1.0f, 1.0f);
	const uint edges = edge_table(cube);
	float3 vertex[12]; // find the vertices where the surface intersects the cube
	if(edges&   1u) vertex[ 0] = interpolate_vertex(p[0], p[1], v[0], v[1], iso); // calculate vertices on all 12 edges
	if(edges&   2u) vertex[ 1] = interpolate_vertex(p[1], p[2], v[1], v[2], iso);
	if(edges&   4u) vertex[ 2] = interpolate_vertex(p[2], p[3], v[2], v[3], iso);
	if(edges&   8u) vertex[ 3] = interpolate_vertex(p[3], p[0], v[3], v[0], iso);
	if(edges&  16u) vertex[ 4] = interpolate_vertex(p[4], p[5], v[4], v[5], iso);
	if(edges&  32u) vertex[ 5] = interpolate_vertex(p[5], p[6], v[5], v[6], iso);
	if(edges&  64u) vertex[ 6] = interpolate_vertex(p[6], p[7], v[6], v[7], iso);
	if(edges& 128u) vertex[ 7] = interpolate_vertex(p[7], p[4], v[7], v[4], iso);
	if(edges& 256u) vertex[ 8] = interpolate_vertex(p[0], p[4], v[0], v[4], iso);
	if(edges& 512u) vertex[ 9] = interpolate_vertex(p[1], p[5], v[1], v[5], iso);
	if(edges&1024u) vertex[10] = interpolate_vertex(p[2], p[6], v[2], v[6], iso);
	if(edges&2048u) vertex[11] = interpolate_vertex(p[3], p[7], v[3], v[7], iso);
	cube *= 15u;
	uint i; // number of triangle vertices
	for(i=0u; i<15u&&triangle_table(cube+i)!=15u; i+=3u) { // create the triangles
		triangles[i   ] = vertex[triangle_table(cube+i   )];
		triangles[i+1u] = vertex[triangle_table(cube+i+1u)];
		triangles[i+2u] = vertex[triangle_table(cube+i+2u)];
	}
	return i/3u; // return number of triangles
}
)+R(typedef struct __attribute__((packed)) struct_ray {
	float3 origin;
	float3 direction;
} ray;
)+R(float intersect_sphere(const ray r, const float3 center, const float radius) {
	const float3 oc = center-r.origin;
	const float b=dot(oc, r.direction), c=sq(b)-dot(oc, oc)+sq(radius);
	return c<0.0f ? -1.0f : b-sqrt(c);
}
)+R(float intersect_sphere_inside(const ray r, const float3 center, const float radius) {
	const float3 oc = center-r.origin;
	const float b=dot(oc, r.direction), c=sq(b)-dot(oc, oc)+sq(radius);
	return c<0.0f ? -1.0f : b+sqrt(c);
}
)+R(float intersect_triangle(const ray r, const float3 p0, const float3 p1, const float3 p2) { // Moeller-Trumbore algorithm
	const float3 u=p1-p0, v=p2-p0, w=r.origin-p0, h=cross(r.direction, v), q=cross(w, u);
	const float f=1.0f/dot(u, h), s=f*dot(w, h), t=f*dot(r.direction, q);
	return (f<0.0f||s<0.0f||s>1.0f||t<0.0f||s+t>1.0f) ? -1.0f : f*dot(v, q);
}
)+R(float intersect_triangle_bidirectional(const ray r, const float3 p0, const float3 p1, const float3 p2) { // Moeller-Trumbore algorithm
	const float3 u=p1-p0, v=p2-p0, w=r.origin-p0, h=cross(r.direction, v), q=cross(w, u);
	const float f=1.0f/dot(u, h), s=f*dot(w, h), t=f*dot(r.direction, q);
	return (s<0.0f||s>1.0f||t<0.0f||s+t>1.0f) ? -1.0f : f*dot(v, q);
}
)+R(float intersect_rhombus(const ray r, const float3 p0, const float3 p1, const float3 p2) { // Moeller-Trumbore algorithm
	const float3 u=p1-p0, v=p2-p0, w=r.origin-p0, h=cross(r.direction, v), q=cross(w, u);
	const float f=1.0f/dot(u, h), s=f*dot(w, h), t=f*dot(r.direction, q);
	return (f<0.0f||s<0.0f||s>1.0f||t<0.0f||t>1.0f) ? -1.0f : f*dot(v, q);
}
)+R(float intersect_plane(const ray r, const float3 p0, const float3 p1, const float3 p2) { // ray-triangle intersection, but skip barycentric coordinates
	const float3 u=p1-p0, v=p2-p0, w=r.origin-p0, h=cross(r.direction, v);
	const float f = 1.0f/dot(u, h);
	return f<0.0f ? -1.0f : f*dot(v, cross(w, u));
}
)+R(float intersect_plane_always(const ray r, const float3 p0, const float3 p1, const float3 p2) { // ray-triangle intersection, but skip barycentric coordinates and visibility check
	const float3 u=p1-p0, v=p2-p0, w=r.origin-p0, h=cross(r.direction, v);
	return dot(v, cross(w, u))/dot(u, h);
}
)+R(bool intersect_cuboid_bool(const ray r, const float3 center, const float Lx, const float Ly, const float Lz) {
	const float3 bmin = center-0.5f*(float3)(Lx, Ly, Lz);
	const float3 bmax = center+0.5f*(float3)(Lx, Ly, Lz);
	const float txa = (bmin.x-r.origin.x)/r.direction.x;
	const float txb = (bmax.x-r.origin.x)/r.direction.x;
	const float txmin = fmin(txa, txb);
	const float txmax = fmax(txa, txb);
	const float tya = (bmin.y-r.origin.y)/r.direction.y;
	const float tyb = (bmax.y-r.origin.y)/r.direction.y;
	const float tymin = fmin(tya, tyb);
	const float tymax = fmax(tya, tyb);
	if(txmin>tymax||tymin>txmax) return false;
	const float tza = (bmin.z-r.origin.z)/r.direction.z;
	const float tzb = (bmax.z-r.origin.z)/r.direction.z;
	const float tzmin = fmin(tza, tzb);
	const float tzmax = fmax(tza, tzb);
	return fmax(txmin, tymin)<=tzmax&&tzmin<=fmin(txmax, tymax);
}
)+R(float intersect_cuboid(const ray r, const float3 center, const float Lx, const float Ly, const float Lz) {
	const float3 bmin = center-0.5f*(float3)(Lx, Ly, Lz);
	const float3 bmax = center+0.5f*(float3)(Lx, Ly, Lz);
	if(r.origin.x>=bmin.x&&r.origin.y>=bmin.y&&r.origin.z>=bmin.z&&r.origin.x<=bmax.x&&r.origin.y<=bmax.y&&r.origin.z<=bmax.z) return 0.0f; // ray origin is within cuboid
	float3 p[8]; // 8 cuboid vertices
	p[0] = (float3)(bmin.x, bmin.y, bmin.z);
	p[1] = (float3)(bmax.x, bmin.y, bmin.z);
	p[2] = (float3)(bmax.x, bmin.y, bmax.z);
	p[3] = (float3)(bmin.x, bmin.y, bmax.z);
	p[4] = (float3)(bmin.x, bmax.y, bmin.z);
	p[5] = (float3)(bmax.x, bmax.y, bmin.z);
	p[6] = (float3)(bmax.x, bmax.y, bmax.z);
	p[7] = (float3)(bmin.x, bmax.y, bmax.z);
	float intersect = -1.0f;
	intersect = fmax(intersect, intersect_rhombus(r, p[0], p[3], p[4])); // test for intersections with the 6 cuboid faces
	intersect = fmax(intersect, intersect_rhombus(r, p[3], p[2], p[7])); // ray will intersect with either 0 or 1 rhombuses
	intersect = fmax(intersect, intersect_rhombus(r, p[2], p[1], p[6]));
	intersect = fmax(intersect, intersect_rhombus(r, p[1], p[0], p[5]));
	intersect = fmax(intersect, intersect_rhombus(r, p[7], p[6], p[4]));
	intersect = fmax(intersect, intersect_rhombus(r, p[1], p[2], p[0]));
	return intersect;
}
)+R(float3 reflect(const float3 direction, const float3 normal) {
	return direction-2.0f*dot(direction, normal)*normal;
}
)+R(float3 refract(const float3 direction, const float3 normal, const float n) {
	const float direction_normal = dot(direction, normal);
	const float sqrt_part = sq(n)-1.0f+sq(direction_normal);
	return sqrt_part>=0.0f ? (direction-(direction_normal+sqrt(sqrt_part))*normal)/n : direction-2.0f*direction_normal*normal; // refraction : total internal reflection
}
)+R(ray get_camray(const int x, const int y, const float* camera_cache) {
	const float zoom = camera_cache[0]; // fetch camera parameters (rotation matrix, camera position, etc.)
	const float dis  = camera_cache[1];
	const float posx = camera_cache[2];
	const float posy = camera_cache[3];
	const float posz = camera_cache[4];
	const float Rxx  = camera_cache[5];
	const float Rxy  = camera_cache[6];
	const float Rxz  = camera_cache[7];
	const float Ryx  = camera_cache[8];
	const float Ryy  = camera_cache[9];
	const float Ryz  = camera_cache[10];
	const float Rzx  = camera_cache[11];
	const float Rzy  = camera_cache[12];
	const float Rzz  = camera_cache[13];
	const bool  vr  = (as_int(camera_cache[14])>>31)&0x1;
	const float rtv = (as_int(camera_cache[14])>>30)&0x1 ? 2.0f : 1.0f;
	const float eye_distance = vload_half(28, (half*)camera_cache);
	const float stereo = (x<(int)def_screen_width/2 ? -1.0f : 1.0f);
	float3 p0 = (float3)(!vr ? 0.0f : stereo*eye_distance/zoom, 0.0f, dis/zoom);
	float3 p1 = p0+normalize((float3)(!vr ? (float)(x-(int)def_screen_width/2) : ((float)(x-(int)def_screen_width/2)-stereo*(float)(def_screen_width/4u))*rtv-stereo*eye_distance, (float)(y-(int)def_screen_height/2), -dis));
	const float x0 = Rxx*p0.x+Ryx*p0.y+Rzx*p0.z; // reverse rotate p0
	const float y0 = Rxy*p0.x+Ryy*p0.y+Rzy*p0.z;
	const float z0 = Rxz*p0.x+Ryz*p0.y+Rzz*p0.z;
	const float x1 = Rxx*p1.x+Ryx*p1.y+Rzx*p1.z; // reverse rotate p1
	const float y1 = Rxy*p1.x+Ryy*p1.y+Rzy*p1.z;
	const float z1 = Rxz*p1.x+Ryz*p1.y+Rzz*p1.z;
	p0 = (float3)(x0, y0, z0)-(float3)(def_domain_offset_x, def_domain_offset_y, def_domain_offset_z);
	p1 = (float3)(x1, y1, z1)-(float3)(def_domain_offset_x, def_domain_offset_y, def_domain_offset_z);
	p0.x=p0.x+posx; p0.y=p0.y+posy; p0.z=p0.z+posz; // reverse transformation of p0
	p1.x=p1.x+posx; p1.y=p1.y+posy; p1.z=p1.z+posz; // reverse transformation of p1
	ray camray;
	camray.origin = p0;
	camray.direction = p1-p0;
	return camray;
}
)+R(uint skybox_bottom(const ray r, const int c1, const int c2, const uint skybox_color) {
	const float3 p0=(float3)(0.0f, 0.0f, -0.5f*(float)def_Nz), p1=(float3)(1.0f, 0.0f, -0.5f*(float)def_Nz), p2=(float3)(0.0f, 1.0f, -0.5f*(float)def_Nz);
	const float distance = intersect_plane(r, p0, p1, p2);
	if(distance>0.0f) { // ray intersects with bottom
		const float3 normal = normalize(cross(p1-p0, p2-p0));
		float3 intersection = r.origin+distance*r.direction;
		const float scale = 2.0f/fmin((float)def_Nx, (float)def_Ny);
		int a = abs((int)floor(scale*intersection.x));
		int b = abs((int)floor(scale*intersection.y));
		const float r = scale*sqrt(sq(intersection.x)+sq(intersection.y));
		const int w = (a%2==b%2);
		return color_mix(w*c1+(1-w)*c2, color_mix(c1, c2, 0.5f), clamp(10.0f/r, 0.0f, 1.0f));
	} else {
		return skybox_color;
	}
}
)+R(uint skybox_color_bw(const float x, const float y) {
	return color_dim(0xFFFFFF, 1.0f-y);
}
)+R(uint skybox_color_hsv(const float x, const float y) {
	const float h = fmod(x*360.0f+120.0f, 360.0f);
	const float s = y>0.5f ? 1.0f : 2.0f*y;
	const float v = y>0.5f ? 2.0f-2.0f*y : 1.0f;
	return hsv_to_rgb(h, s, v);
}
)+R(uint skybox_color_sunset(const float x, const float y) {
	return color_mix(255<<16|175<<8|55, y<0.5f ? 55<<16|111<<8|255 : 0, 2.0f*(0.5f-fabs(y-0.5f)));
}
)+R(uint skybox_color_grid(const float x, const float y, const int c1, const int c2) {
	int a = (int)(72.0f*x);
	int b = (int)(36.0f*y);
	const int w = (a%2==b%2);
	return w*c1+(1-w)*c2;
}
)+R(uint skybox_color(const ray r, const global int* skybox) {
	const float3 direction = normalize(r.direction); // to avoid artifacts from asin(direction.z)
	//const float x = fma(atan2(direction.x, direction.y),  0.5f/3.1415927f, 0.5f);
	//const float y = fma(asin (direction.z             ), -1.0f/3.1415927f, 0.5f);
	//return skybox_color_bw(x, y);
	//return color_mix(skybox_color_hsv(x, y), skybox_color_grid(x, y, 0xFFFFFF, 0x000000), 0.95f-0.33f*(2.0f*(0.5f-fabs(y-0.5f))));
	//return skybox_bottom(r, 0xFFFFFF, 0xF0F0F0, skybox_color_grid(x, y, 0xFFFFFF, 0xF0F0F0));
	const float fu = (float)def_skybox_width *fma(atan2(direction.x, direction.y),  0.5f/3.1415927f, 0.5f);
	const float fv = (float)def_skybox_height*fma(asin (direction.z             ), -1.0f/3.1415927f, 0.5f);
	const int ua=clamp((int)fu, 0, (int)def_skybox_width-1), va=clamp((int)fv, 0, (int)def_skybox_height-1), ub=(ua+1)%def_skybox_width, vb=min(va+1, (int)def_skybox_height-1); // bilinear interpolation positions
	const uint s00=skybox[ua+va*def_skybox_width], s01=skybox[ua+vb*def_skybox_width], s10=skybox[ub+va*def_skybox_width], s11=skybox[ub+vb*def_skybox_width];
	const float u1=fu-(float)ua, v1=fv-(float)va, u0=1.0f-u1, v0=1.0f-v1; // interpolation factors
	return color_mix(color_mix(s00, s01, v0), color_mix(s10, s11, v0), u0); // perform bilinear interpolation
}
)+R(uint last_ray_reflectivity(const ray reflection, const ray transmission, const float reflectivity, const float transmissivity, const global int* skybox) {
	return color_mix(skybox_color(reflection, skybox), color_mix(skybox_color(transmission, skybox), def_absorption_color, transmissivity), reflectivity);
}
)+R(float ray_grid_traverse(const ray r, const global float* phi, const global uchar* flags, float3* normal, const uint Nx, const uint Ny, const uint Nz) {
	const float3 pa = r.origin;
	const float3 pb = r.origin+r.direction;
	const float xa=pa.x-0.5f+0.5f*(float)Nx, ya=pa.y-0.5f+0.5f*(float)Ny, za=pa.z-0.5f+0.5f*(float)Nz; // start point
	const float xb=pb.x-0.5f+0.5f*(float)Nx, yb=pb.y-0.5f+0.5f*(float)Ny, zb=pb.z-0.5f+0.5f*(float)Nz; // end point
	const int dx=(int)sign(xb-xa), dy=(int)sign(yb-ya), dz=(int)sign(zb-za); // fast ray-grid-traversal
	const float fxa=xa-floor(xa), fya=ya-floor(ya), fza=za-floor(za);
	int3 xyz = (int3)(floor(xa), floor(ya), floor(za));
	const float tdx = fmin((float)dx/(xb-xa), 1E7f);
	const float tdy = fmin((float)dy/(yb-ya), 1E7f);
	const float tdz = fmin((float)dz/(zb-za), 1E7f);
	float tmx = tdx*(dx>0 ? 1.0f-fxa : fxa);
	float tmy = tdy*(dy>0 ? 1.0f-fya : fya);
	float tmz = tdz*(dz>0 ? 1.0f-fza : fza);
	while(true) {
		if(tmx<tmy) {
			if(tmx<tmz) { xyz.x += dx; tmx += tdx; } else { xyz.z += dz; tmz += tdz; }
		} else {
			if(tmy<tmz) { xyz.y += dy; tmy += tdy; } else { xyz.z += dz; tmz += tdz; }
		}
		if(xyz.x<-1 || xyz.y<-1 || xyz.z<-1 || xyz.x>=(int)Nx || xyz.y>=(int)Ny || xyz.z>=(int)Nz) break;
		else if(xyz.x<0 || xyz.y<0 || xyz.z<0 || xyz.x>=(int)Nx-1 || xyz.y>=(int)Ny-1 || xyz.z>=(int)Nz-1) continue;
		const uint x0 =  (uint)xyz.x; // cube stencil
		const uint xp =  (uint)xyz.x+1u;
		const uint y0 =  (uint)xyz.y    *Nx;
		const uint yp = ((uint)xyz.y+1u)*Nx;
		const uint z0 =  (uint)xyz.z    *Ny*Nx;
		const uint zp = ((uint)xyz.z+1u)*Ny*Nx;
		uint j[8];
		j[0] = x0+y0+z0; // 000
		j[1] = xp+y0+z0; // +00
		j[2] = xp+y0+zp; // +0+
		j[3] = x0+y0+zp; // 00+
		j[4] = x0+yp+z0; // 0+0
		j[5] = xp+yp+z0; // ++0
		j[6] = xp+yp+zp; // +++
		j[7] = x0+yp+zp; // 0++
		uchar flags_cell = 0u; // check with cheap flags if the isosurface goes through the current marching-cubes cell (~15% performance boost)
		for(uint i=0u; i<8u; i++) flags_cell |= flags[j[i]];
		if(!(flags_cell&(TYPE_S|TYPE_E|TYPE_I))) continue; // cell is entirely inside/outside of the isosurface
		float v[8];
		for(uint i=0u; i<8u; i++) v[i] = phi[j[i]];
		float3 triangles[15]; // maximum of 5 triangles with 3 vertices each
		const uint tn = marching_cubes(v, 0.5f, triangles); // run marching cubes algorithm
		if(tn==0u) continue; // if returned tn value is non-zero, iterate through triangles
		const float3 offset = (float3)((float)xyz.x+0.5f-0.5f*(float)Nx, (float)xyz.y+0.5f-0.5f*(float)Ny, (float)xyz.z+0.5f-0.5f*(float)Nz);
		for(uint i=0u; i<tn; i++) {
			const float3 p0 = triangles[3u*i   ]+offset;
			const float3 p1 = triangles[3u*i+1u]+offset;
			const float3 p2 = triangles[3u*i+2u]+offset;
			const float intersect = intersect_triangle_bidirectional(r, p0, p1, p2); // for each triangle, check ray-triangle intersection
			if(intersect>0.0f) { // intersection found (there can only be exactly 1 intersection)
				const uint xq =  ((uint)xyz.x   +2u)%Nx; // central difference stencil on each cube corner point
				const uint xm =  ((uint)xyz.x+Nx-1u)%Nx;
				const uint yq = (((uint)xyz.y   +2u)%Ny)*Nx;
				const uint ym = (((uint)xyz.y+Ny-1u)%Ny)*Nx;
				const uint zq = (((uint)xyz.z   +2u)%Nz)*Ny*Nx;
				const uint zm = (((uint)xyz.z+Nz-1u)%Nz)*Ny*Nx;
				float3 n[8];
				n[0] = (float3)(phi[xm+y0+z0]-v[1], phi[x0+ym+z0]-v[4], phi[x0+y0+zm]-v[3]); // central difference stencil on each cube corner point
				n[1] = (float3)(v[0]-phi[xq+y0+z0], phi[xp+ym+z0]-v[5], phi[xp+y0+zm]-v[2]); // compute normal vectors from gradient
				n[2] = (float3)(v[3]-phi[xq+y0+zp], phi[xp+ym+zp]-v[6], v[1]-phi[xp+y0+zq]); // normalize later during trilinear interpolation more efficiently
				n[3] = (float3)(phi[xm+y0+zp]-v[2], phi[x0+ym+zp]-v[7], v[0]-phi[x0+y0+zq]);
				n[4] = (float3)(phi[xm+yp+z0]-v[5], v[0]-phi[x0+yq+z0], phi[x0+yp+zm]-v[7]);
				n[5] = (float3)(v[4]-phi[xq+yp+z0], v[1]-phi[xp+yq+z0], phi[xp+yp+zm]-v[6]);
				n[6] = (float3)(v[7]-phi[xq+yp+zp], v[2]-phi[xp+yq+zp], v[5]-phi[xp+yp+zq]);
				n[7] = (float3)(phi[xm+yp+zp]-v[6], v[3]-phi[x0+yq+zp], v[4]-phi[x0+yp+zq]);
				const float3 p = r.origin+intersect*r.direction-offset; // intersection point minus offset
				const float x1=p.x-floor(p.x), y1=p.y-floor(p.y), z1=p.z-floor(p.z), x0=1.0f-x1, y0=1.0f-y1, z0=1.0f-z1; // calculate interpolation factors
				*normal = normalize(
					(x0*y0*z0*rsqrt(fma(n[0].x, n[0].x, fma(n[0].y, n[0].y, fma(n[0].z, n[0].z, 1E-9f)))))*n[0]+
					(x1*y0*z0*rsqrt(fma(n[1].x, n[1].x, fma(n[1].y, n[1].y, fma(n[1].z, n[1].z, 1E-9f)))))*n[1]+
					(x1*y0*z1*rsqrt(fma(n[2].x, n[2].x, fma(n[2].y, n[2].y, fma(n[2].z, n[2].z, 1E-9f)))))*n[2]+
					(x0*y0*z1*rsqrt(fma(n[3].x, n[3].x, fma(n[3].y, n[3].y, fma(n[3].z, n[3].z, 1E-9f)))))*n[3]+
					(x0*y1*z0*rsqrt(fma(n[4].x, n[4].x, fma(n[4].y, n[4].y, fma(n[4].z, n[4].z, 1E-9f)))))*n[4]+
					(x1*y1*z0*rsqrt(fma(n[5].x, n[5].x, fma(n[5].y, n[5].y, fma(n[5].z, n[5].z, 1E-9f)))))*n[5]+
					(x1*y1*z1*rsqrt(fma(n[6].x, n[6].x, fma(n[6].y, n[6].y, fma(n[6].z, n[6].z, 1E-9f)))))*n[6]+
					(x0*y1*z1*rsqrt(fma(n[7].x, n[7].x, fma(n[7].y, n[7].y, fma(n[7].z, n[7].z, 1E-9f)))))*n[7]
				); // perform normalization and trilinear interpolation
				return intersect; // intersection found, exit loop, process transmission ray
			}
		}
	}
	return -1.0f; // no intersection found
}
)+R(bool raytrace_phi_mirror(const ray ray_in, ray* ray_reflect, const global float* phi, const global uchar* flags, const global int* skybox, const uint Nx, const uint Ny, const uint Nz) { // only reflection
	float3 normal;
	float d = ray_grid_traverse(ray_in, phi, flags, &normal, Nx, Ny, Nz); // move ray through lattice, at each cell call marching_cubes
	if(d==-1.0f) return false; // no intersection found
	ray_reflect->origin = ray_in.origin+(d-0.0003163f)*ray_in.direction; // start intersection points a bit in front triangle to avoid self-reflection
	ray_reflect->direction = reflect(ray_in.direction, normal);
	return true;
}
)+R(bool raytrace_phi(const ray ray_in, ray* ray_reflect, ray* ray_transmit, float* reflectivity, float* transmissivity, const global float* phi, const global uchar* flags, const global int* skybox, const uint Nx, const uint Ny, const uint Nz) {
	float3 normal;
	float d = ray_grid_traverse(ray_in, phi, flags, &normal, Nx, Ny, Nz); // move ray through lattice, at each cell call marching_cubes
	if(d==-1.0f) return false; // no intersection found
	const float ray_in_normal = dot(ray_in.direction, normal);
	const bool is_inside = ray_in_normal>0.0f; // camera is in fluid
	ray_reflect->origin = ray_in.origin+(d-0.0003163f)*ray_in.direction; // start intersection points a bit in front triangle to avoid self-reflection
	ray_reflect->direction = reflect(ray_in.direction, normal); // compute reflection ray
	ray ray_internal; // compute internal ray and transmission ray
	ray_internal.origin = ray_in.origin+(d+0.0003163f)*ray_in.direction; // start intersection points a bit behind triangle to avoid self-transmission
	ray_internal.direction = refract(ray_in.direction, normal, def_n);
	const float wr = clamp(sq(cb(2.0f*acospi(fabs(ray_in_normal)))), 0.0f, 1.0f); // increase reflectivity if ray intersects surface at shallow angle
	if(is_inside) { // swap ray_reflect and ray_internal
		const float3 ray_internal_origin = ray_internal.origin;
		ray_internal.origin = ray_reflect->origin;
		ray_internal.direction = ray_reflect->direction;
		ray_reflect->origin = ray_internal_origin; // re-use internal ray origin
		ray_reflect->direction = refract(ray_in.direction, -normal, 1.0f/def_n); // compute refraction again: refract out of fluid
		if(sq(1.0f/def_n)-1.0f+sq(ray_in_normal)>=0.0f) { // refraction through Snell's window
			ray_transmit->origin = ray_reflect->origin; // reflection ray and transmission ray are the same
			ray_transmit->direction = ray_reflect->direction;
			*reflectivity = 0.0f;
			*transmissivity = exp(def_attenuation*d); // Beer-Lambert law
			return true;
		}
	}
	float d_internal = d;
	d = ray_grid_traverse(ray_internal, phi, flags, &normal, Nx, Ny, Nz); // 2nd ray-grid traversal call: refraction (camera outside) or total internal reflection (camera inside)
	ray_transmit->origin = d!=-1.0f ? ray_internal.origin+(d+0.0003163f)*ray_internal.direction : ray_internal.origin; // start intersection points a bit behind triangle to avoid self-transmission
	ray_transmit->direction = d!=-1.0f ? refract(ray_internal.direction, -normal, 1.0f/def_n) : ray_internal.direction; // internal ray intersects isosurface : internal ray does not intersect again
	*reflectivity = is_inside ? 0.0f : wr; // is_inside means camera is inside fluid, so this is a total internal reflection down here
	*transmissivity = d!=-1.0f ? exp(def_attenuation*((float)is_inside*d_internal+d)) : (float)(def_attenuation==0.0f); // Beer-Lambert law
	return true;
}
)+R(bool is_above_plane(const float3 point, const float3 plane_p, const float3 plane_n) {
	return dot(point-plane_p, plane_n)>=0.0f;
}
)+R(bool is_below_plane(const float3 point, const float3 plane_p, const float3 plane_n) {
	return dot(point-plane_p, plane_n)<=0.0f;
}
)+R(bool is_in_camera_frustrum(const float3 p, const float* camera_cache) { // returns true if point is located in camera frustrum
	const float vr = (float)((as_int(camera_cache[14])>>31)&0x1);
	const ray r00 = get_camray(0,                       0                       , camera_cache); // get 4 edge vectors of frustrum
	const ray r01 = get_camray((int)def_screen_width-1, 0                       , camera_cache);
	const ray r10 = get_camray(0,                       (int)def_screen_height-1, camera_cache);
	const ray r11 = get_camray((int)def_screen_width-1, (int)def_screen_height-1, camera_cache);
	const float3 plane_n_top    = cross(r00.direction, r01.direction); // get 4 frustrum planes
	const float3 plane_n_bottom = cross(r11.direction, r10.direction);
	const float3 plane_n_left   = cross(r10.direction, r00.direction);
	const float3 plane_n_right  = cross(r01.direction, r11.direction);
	const float3 plane_p_top    = r00.origin-2.0f*plane_n_top; // move frustrum planes outward by 2 units
	const float3 plane_p_bottom = r11.origin-2.0f*plane_n_bottom;
	const float3 plane_p_left   = r00.origin-(2.0f+16.0f*vr)*plane_n_left; // move frustrum planes outward by 2 units, for stereoscopic rendering a bit more
	const float3 plane_p_right  = r11.origin-(2.0f+16.0f*vr)*plane_n_right;
	return is_above_plane(p, plane_p_top, plane_n_top)&&is_above_plane(p, plane_p_bottom, plane_n_bottom)&&is_above_plane(p, plane_p_left, plane_n_left)&&is_above_plane(p, plane_p_right, plane_n_right);
}
)+"#endif"+R( // GRAPHICS



// ################################################## LBM code ##################################################

)+R(uint3 coordinates(const uint n) { // disassemble 1D index to 3D coordinates (n -> x,y,z)
	const uint t = n%(def_Nx*def_Ny);
	return (uint3)(t%def_Nx, t/def_Nx, n/(def_Nx*def_Ny)); // n = x+(y+z*Ny)*Nx
}
)+R(uint index(const uint3 xyz) { // assemble 1D index from 3D coordinates (x,y,z -> n)
	return xyz.x+(xyz.y+xyz.z*def_Ny)*def_Nx; // n = x+(y+z*Ny)*Nx
}
)+R(float3 position(const uint3 xyz) { // 3D coordinates to 3D position
	return (float3)((float)xyz.x+0.5f-0.5f*(float)def_Nx, (float)xyz.y+0.5f-0.5f*(float)def_Ny, (float)xyz.z+0.5f-0.5f*(float)def_Nz);
}
)+R(float3 mirror_position(const float3 p) { // mirror position into periodic boundaries
	float3 r;
	r.x = sign(p.x)*(fmod(fabs(p.x)+0.5f*(float)def_Nx, (float)def_Nx)-0.5f*(float)def_Nx);
	r.y = sign(p.y)*(fmod(fabs(p.y)+0.5f*(float)def_Ny, (float)def_Ny)-0.5f*(float)def_Ny);
	r.z = sign(p.z)*(fmod(fabs(p.z)+0.5f*(float)def_Nz, (float)def_Nz)-0.5f*(float)def_Nz);
	return r;
}
)+R(float3 mirror_distance(const float3 d) { // mirror distance vector into periodic boundaries
	return mirror_position(d);
}
)+R(bool is_halo(const uint n) {
	const uint3 xyz = coordinates(n);
	return ((def_Dx>1u)&(xyz.x==0u||xyz.x>=def_Nx-1u))||((def_Dy>1u)&(xyz.y==0u||xyz.y>=def_Ny-1u))||((def_Dz>1u)&(xyz.z==0u||xyz.z>=def_Nz-1u));
}
)+R(bool is_halo_q(const uint n) {
	const uint3 xyz = coordinates(n);
)+"#ifndef SURFACE"+R(
	return ((def_Dx>1u)&(xyz.x<=1u||xyz.x>=def_Nx-3u))||((def_Dy>1u)&(xyz.y<=1u||xyz.y>=def_Ny-3u))||((def_Dz>1u)&(xyz.z<=1u||xyz.z>=def_Nz-3u));
)+"#else"+R( // SURFACE
	return ((def_Dx>1u)&(xyz.x==0u||xyz.x>=def_Nx-2u))||((def_Dy>0u)&(xyz.y==1u||xyz.y>=def_Ny-2u))||((def_Dz>1u)&(xyz.z==0u||xyz.z>=def_Nz-2u)); // halo data is kept up-to-date with SURFACE extension, so allow using halo data for rendering
)+"#endif"+R( // SURFACE
} // is_halo_q()

)+R(float half_to_float_custom(const ushort x) { // custom 16-bit floating-point format, 1-4-11, exp-15, +-1.99951168, +-6.10351562E-5, +-2.98023224E-8, 3.612 digits
	const uint e = (x&0x7800)>>11; // exponent
	const uint m = (x&0x07FF)<<12; // mantissa
	const uint v = as_uint((float)m)>>23; // evil log2 bit hack to count leading zeros in denormalized format
	return as_float((x&0x8000)<<16 | (e!=0)*((e+112)<<23|m) | ((e==0)&(m!=0))*((v-37)<<23|((m<<(150-v))&0x007FF000))); // sign : normalized : denormalized
}
)+R(ushort float_to_half_custom(const float x) { // custom 16-bit floating-point format, 1-4-11, exp-15, +-1.99951168, +-6.10351562E-5, +-2.98023224E-8, 3.612 digits
	const uint b = as_uint(x)+0x00000800; // round-to-nearest-even: add last bit after truncated mantissa
	const uint e = (b&0x7F800000)>>23; // exponent
	const uint m = b&0x007FFFFF; // mantissa; in line below: 0x007FF800 = 0x00800000-0x00000800 = decimal indicator flag - initial rounding
	return (b&0x80000000)>>16 | (e>112)*((((e-112)<<11)&0x7800)|m>>12) | ((e<113)&(e>100))*((((0x007FF800+m)>>(124-e))+1)>>1); // sign : normalized : denormalized (assume [-2,2])
}

)+R(ulong index_f(const uint n, const uint i) { // 64-bit indexing (maximum 2^32 lattice points (1624^3 lattice resolution, 225GB)
	return (ulong)i*def_N+(ulong)n; // SoA (229% faster on GPU)
}
)+R(float c(const uint i) { // avoid constant keyword by encapsulating data in function which gets inlined by compiler
	const float c[3u*def_velocity_set] = {
)+"#if defined(D2Q9)"+R(
		0, 1,-1, 0, 0, 1,-1, 1,-1, // x
		0, 0, 0, 1,-1, 1,-1,-1, 1, // y
		0, 0, 0, 0, 0, 0, 0, 0, 0  // z
)+"#elif defined(D3Q15)"+R(
		0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, // x
		0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1, // y
		0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1, 1,-1  // z
)+"#elif defined(D3Q19)"+R(
		0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, // x
		0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, // y
		0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1  // z
)+"#elif defined(D3Q27)"+R(
		0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, // x
		0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1, // y
		0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1  // z
)+"#endif"+R( // D3Q27
	};
	return c[i];
}
)+R(float w(const uint i) { // avoid constant keyword by encapsulating data in function which gets inlined by compiler
	const float w[def_velocity_set] = { def_w0, // velocity set weights
)+"#if defined(D2Q9)"+R(
		def_ws, def_ws, def_ws, def_ws, def_we, def_we, def_we, def_we
)+"#elif defined(D3Q15)"+R(
		def_ws, def_ws, def_ws, def_ws, def_ws, def_ws,
		def_wc, def_wc, def_wc, def_wc, def_wc, def_wc, def_wc, def_wc
)+"#elif defined(D3Q19)"+R(
		def_ws, def_ws, def_ws, def_ws, def_ws, def_ws,
		def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we
)+"#elif defined(D3Q27)"+R(
		def_ws, def_ws, def_ws, def_ws, def_ws, def_ws,
		def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we,
		def_wc, def_wc, def_wc, def_wc, def_wc, def_wc, def_wc, def_wc
)+"#endif"+R( // D3Q27
	};
	return w[i];
}
)+R(void calculate_indices(const uint n, uint* x0, uint* xp, uint* xm, uint* y0, uint* yp, uint* ym, uint* z0, uint* zp, uint* zm) {
	const uint3 xyz = coordinates(n);
	*x0 =   xyz.x; // pre-calculate indices (periodic boundary conditions)
	*xp =  (xyz.x       +1u)%def_Nx;
	*xm =  (xyz.x+def_Nx-1u)%def_Nx;
	*y0 =   xyz.y                   *def_Nx;
	*yp = ((xyz.y       +1u)%def_Ny)*def_Nx;
	*ym = ((xyz.y+def_Ny-1u)%def_Ny)*def_Nx;
	*z0 =   xyz.z                   *def_Ny*def_Nx;
	*zp = ((xyz.z       +1u)%def_Nz)*def_Ny*def_Nx;
	*zm = ((xyz.z+def_Nz-1u)%def_Nz)*def_Ny*def_Nx;
} // calculate_indices()
)+R(void neighbors(const uint n, uint* j) { // calculate neighbor indices
	uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
	calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
	j[0] = n;
)+"#if defined(D2Q9)"+R(
	j[ 1] = xp+y0; j[ 2] = xm+y0; // +00 -00
	j[ 3] = x0+yp; j[ 4] = x0+ym; // 0+0 0-0
	j[ 5] = xp+yp; j[ 6] = xm+ym; // ++0 --0
	j[ 7] = xp+ym; j[ 8] = xm+yp; // +-0 -+0
)+"#elif defined(D3Q15)"+R(
	j[ 1] = xp+y0+z0; j[ 2] = xm+y0+z0; // +00 -00
	j[ 3] = x0+yp+z0; j[ 4] = x0+ym+z0; // 0+0 0-0
	j[ 5] = x0+y0+zp; j[ 6] = x0+y0+zm; // 00+ 00-
	j[ 7] = xp+yp+zp; j[ 8] = xm+ym+zm; // +++ ---
	j[ 9] = xp+yp+zm; j[10] = xm+ym+zp; // ++- --+
	j[11] = xp+ym+zp; j[12] = xm+yp+zm; // +-+ -+-
	j[13] = xm+yp+zp; j[14] = xp+ym+zm; // -++ +--
)+"#elif defined(D3Q19)"+R(
	j[ 1] = xp+y0+z0; j[ 2] = xm+y0+z0; // +00 -00
	j[ 3] = x0+yp+z0; j[ 4] = x0+ym+z0; // 0+0 0-0
	j[ 5] = x0+y0+zp; j[ 6] = x0+y0+zm; // 00+ 00-
	j[ 7] = xp+yp+z0; j[ 8] = xm+ym+z0; // ++0 --0
	j[ 9] = xp+y0+zp; j[10] = xm+y0+zm; // +0+ -0-
	j[11] = x0+yp+zp; j[12] = x0+ym+zm; // 0++ 0--
	j[13] = xp+ym+z0; j[14] = xm+yp+z0; // +-0 -+0
	j[15] = xp+y0+zm; j[16] = xm+y0+zp; // +0- -0+
	j[17] = x0+yp+zm; j[18] = x0+ym+zp; // 0+- 0-+
)+"#elif defined(D3Q27)"+R(
	j[ 1] = xp+y0+z0; j[ 2] = xm+y0+z0; // +00 -00
	j[ 3] = x0+yp+z0; j[ 4] = x0+ym+z0; // 0+0 0-0
	j[ 5] = x0+y0+zp; j[ 6] = x0+y0+zm; // 00+ 00-
	j[ 7] = xp+yp+z0; j[ 8] = xm+ym+z0; // ++0 --0
	j[ 9] = xp+y0+zp; j[10] = xm+y0+zm; // +0+ -0-
	j[11] = x0+yp+zp; j[12] = x0+ym+zm; // 0++ 0--
	j[13] = xp+ym+z0; j[14] = xm+yp+z0; // +-0 -+0
	j[15] = xp+y0+zm; j[16] = xm+y0+zp; // +0- -0+
	j[17] = x0+yp+zm; j[18] = x0+ym+zp; // 0+- 0-+
	j[19] = xp+yp+zp; j[20] = xm+ym+zm; // +++ ---
	j[21] = xp+yp+zm; j[22] = xm+ym+zp; // ++- --+
	j[23] = xp+ym+zp; j[24] = xm+yp+zm; // +-+ -+-
	j[25] = xm+yp+zp; j[26] = xp+ym+zm; // -++ +--
)+"#endif"+R( // D3Q27
} // neighbors()

)+R(float3 load_u(const uint n, const global float* u) {
	return (float3)(u[n], u[def_N+(ulong)n], u[2ul*def_N+(ulong)n]);
}
)+R(float3 closest_u(const float3 p, const global float* u) { // return velocity of closest lattice point to point p
	const uint x = (uint)(p.x+1.5f*(float)def_Nx)%def_Nx;
	const uint y = (uint)(p.y+1.5f*(float)def_Ny)%def_Ny;
	const uint z = (uint)(p.z+1.5f*(float)def_Nz)%def_Nz;
	const uint n = x+(y+z*def_Ny)*def_Nx;
	return load_u(n, u);
} // closest_u()
)+R(float3 interpolate_u(const float3 p, const global float* u) { // trilinear interpolation of velocity at point p
	const float xa=p.x-0.5f+1.5f*def_Nx, ya=p.y-0.5f+1.5f*def_Ny, za=p.z-0.5f+1.5f*def_Nz; // subtract lattice offsets
	const uint xb=(uint)xa, yb=(uint)ya, zb=(uint)za; // integer casting to find bottom left corner
	const float x1=xa-(float)xb, y1=ya-(float)yb, z1=za-(float)zb, x0=1.0f-x1, y0=1.0f-y1, z0=1.0f-z1; // calculate interpolation factors
	float3 un[8]; // velocities of unit cube corner points
	for(uint c=0u; c<8u; c++) { // count over eight corner points
		const uint i=(c&0x04u)>>2, j=(c&0x02u)>>1, k=c&0x01u; // disassemble c into corner indices ijk
		const uint x=(xb+i)%def_Nx, y=(yb+j)%def_Ny, z=(zb+k)%def_Nz; // calculate corner lattice positions
		const uint n = x+(y+z*def_Ny)*def_Nx; // calculate lattice linear index
		un[c] = load_u(n, u); // load velocity from lattice point
	}
	return (x0*y0*z0)*un[0]+(x0*y0*z1)*un[1]+(x0*y1*z0)*un[2]+(x0*y1*z1)*un[3]+(x1*y0*z0)*un[4]+(x1*y0*z1)*un[5]+(x1*y1*z0)*un[6]+(x1*y1*z1)*un[7]; // perform trilinear interpolation
} // interpolate_u()
)+R(float calculate_Q_cached(const float3* uj) { // Q-criterion
	const float duxdx=uj[0].x-uj[1].x, duydx=uj[0].y-uj[1].y, duzdx=uj[0].z-uj[1].z; // du/dx = (u2-u0)/2
	const float duxdy=uj[2].x-uj[3].x, duydy=uj[2].y-uj[3].y, duzdy=uj[2].z-uj[3].z;
	const float duxdz=uj[4].x-uj[5].x, duydz=uj[4].y-uj[5].y, duzdz=uj[4].z-uj[5].z;
	const float omega_xy=duxdy-duydx, omega_xz=duxdz-duzdx, omega_yz=duydz-duzdy; // antisymmetric tensor, omega_xx = omega_yy = omega_zz = 0
	const float s_xx2=duxdx, s_yy2=duydy, s_zz2=duzdz; // s_xx2 = s_xx/2, s_yy2 = s_yy/2, s_zz2 = s_zz/2
	const float s_xy=duxdy+duydx, s_xz=duxdz+duzdx, s_yz=duydz+duzdy; // symmetric tensor
	const float omega2 = sq(omega_xy)+sq(omega_xz)+sq(omega_yz); // ||omega||_2^2
	const float s2 = 2.0f*(sq(s_xx2)+sq(s_yy2)+sq(s_zz2))+sq(s_xy)+sq(s_xz)+sq(s_yz); // ||s||_2^2
	return 0.25f*(omega2-s2); // Q = 1/2*(||omega||_2^2-||s||_2^2), addidional factor 1/2 from cental finite differences of velocity
} // calculate_Q_cached()
)+R(float calculate_Q(const uint n, const global float* u) { // Q-criterion
	uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
	calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
	uint j[6];
	j[0] = xp+y0+z0; j[1] = xm+y0+z0; // +00 -00
	j[2] = x0+yp+z0; j[3] = x0+ym+z0; // 0+0 0-0
	j[4] = x0+y0+zp; j[5] = x0+y0+zm; // 00+ 00-
	float3 uj[6];
	for(uint i=0u; i<6u; i++) uj[i] = load_u(j[i], u);
	return calculate_Q_cached(uj);
} // calculate_Q()

)+R(void calculate_f_eq(const float rho, float ux, float uy, float uz, float* feq) { // calculate f_equilibrium from density and velocity field (perturbation method / DDF-shifting)
	const float c3=-3.0f*(sq(ux)+sq(uy)+sq(uz)), rhom1=rho-1.0f; // c3 = -2*sq(u)/(2*sq(c)), rhom1 is arithmetic optimization to minimize digit extinction
	ux *= 3.0f;
	uy *= 3.0f;
	uz *= 3.0f;
	feq[ 0] = def_w0*fma(rho, 0.5f*c3, rhom1); // 000 (identical for all velocity sets)
)+"#if defined(D2Q9)"+R(
	const float u0=ux+uy, u1=ux-uy; // these pre-calculations make manual unrolling require less FLOPs
	const float rhos=def_ws*rho, rhoe=def_we*rho, rhom1s=def_ws*rhom1, rhom1e=def_we*rhom1;
	feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
	feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
	feq[ 5] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), u0), rhom1e); feq[ 6] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), -u0), rhom1e); // ++0 --0
	feq[ 7] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), u1), rhom1e); feq[ 8] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), -u1), rhom1e); // +-0 -+0
)+"#elif defined(D3Q15)"+R(
	const float u0=ux+uy+uz, u1=ux+uy-uz, u2=ux-uy+uz, u3=-ux+uy+uz;
	const float rhos=def_ws*rho, rhoc=def_wc*rho, rhom1s=def_ws*rhom1, rhom1c=def_wc*rhom1;
	feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
	feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
	feq[ 5] = fma(rhos, fma(0.5f, fma(uz, uz, c3), uz), rhom1s); feq[ 6] = fma(rhos, fma(0.5f, fma(uz, uz, c3), -uz), rhom1s); // 00+ 00-
	feq[ 7] = fma(rhoc, fma(0.5f, fma(u0, u0, c3), u0), rhom1c); feq[ 8] = fma(rhoc, fma(0.5f, fma(u0, u0, c3), -u0), rhom1c); // +++ ---
	feq[ 9] = fma(rhoc, fma(0.5f, fma(u1, u1, c3), u1), rhom1c); feq[10] = fma(rhoc, fma(0.5f, fma(u1, u1, c3), -u1), rhom1c); // ++- --+
	feq[11] = fma(rhoc, fma(0.5f, fma(u2, u2, c3), u2), rhom1c); feq[12] = fma(rhoc, fma(0.5f, fma(u2, u2, c3), -u2), rhom1c); // +-+ -+-
	feq[13] = fma(rhoc, fma(0.5f, fma(u3, u3, c3), u3), rhom1c); feq[14] = fma(rhoc, fma(0.5f, fma(u3, u3, c3), -u3), rhom1c); // -++ +--
)+"#elif defined(D3Q19)"+R(
	const float u0=ux+uy, u1=ux+uz, u2=uy+uz, u3=ux-uy, u4=ux-uz, u5=uy-uz;
	const float rhos=def_ws*rho, rhoe=def_we*rho, rhom1s=def_ws*rhom1, rhom1e=def_we*rhom1;
	feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
	feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
	feq[ 5] = fma(rhos, fma(0.5f, fma(uz, uz, c3), uz), rhom1s); feq[ 6] = fma(rhos, fma(0.5f, fma(uz, uz, c3), -uz), rhom1s); // 00+ 00-
	feq[ 7] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), u0), rhom1e); feq[ 8] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), -u0), rhom1e); // ++0 --0
	feq[ 9] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), u1), rhom1e); feq[10] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), -u1), rhom1e); // +0+ -0-
	feq[11] = fma(rhoe, fma(0.5f, fma(u2, u2, c3), u2), rhom1e); feq[12] = fma(rhoe, fma(0.5f, fma(u2, u2, c3), -u2), rhom1e); // 0++ 0--
	feq[13] = fma(rhoe, fma(0.5f, fma(u3, u3, c3), u3), rhom1e); feq[14] = fma(rhoe, fma(0.5f, fma(u3, u3, c3), -u3), rhom1e); // +-0 -+0
	feq[15] = fma(rhoe, fma(0.5f, fma(u4, u4, c3), u4), rhom1e); feq[16] = fma(rhoe, fma(0.5f, fma(u4, u4, c3), -u4), rhom1e); // +0- -0+
	feq[17] = fma(rhoe, fma(0.5f, fma(u5, u5, c3), u5), rhom1e); feq[18] = fma(rhoe, fma(0.5f, fma(u5, u5, c3), -u5), rhom1e); // 0+- 0-+
)+"#elif defined(D3Q27)"+R(
	const float u0=ux+uy, u1=ux+uz, u2=uy+uz, u3=ux-uy, u4=ux-uz, u5=uy-uz, u6=ux+uy+uz, u7=ux+uy-uz, u8=ux-uy+uz, u9=-ux+uy+uz;
	const float rhos=def_ws*rho, rhoe=def_we*rho, rhoc=def_wc*rho, rhom1s=def_ws*rhom1, rhom1e=def_we*rhom1, rhom1c=def_wc*rhom1;
	feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
	feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
	feq[ 5] = fma(rhos, fma(0.5f, fma(uz, uz, c3), uz), rhom1s); feq[ 6] = fma(rhos, fma(0.5f, fma(uz, uz, c3), -uz), rhom1s); // 00+ 00-
	feq[ 7] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), u0), rhom1e); feq[ 8] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), -u0), rhom1e); // ++0 --0
	feq[ 9] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), u1), rhom1e); feq[10] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), -u1), rhom1e); // +0+ -0-
	feq[11] = fma(rhoe, fma(0.5f, fma(u2, u2, c3), u2), rhom1e); feq[12] = fma(rhoe, fma(0.5f, fma(u2, u2, c3), -u2), rhom1e); // 0++ 0--
	feq[13] = fma(rhoe, fma(0.5f, fma(u3, u3, c3), u3), rhom1e); feq[14] = fma(rhoe, fma(0.5f, fma(u3, u3, c3), -u3), rhom1e); // +-0 -+0
	feq[15] = fma(rhoe, fma(0.5f, fma(u4, u4, c3), u4), rhom1e); feq[16] = fma(rhoe, fma(0.5f, fma(u4, u4, c3), -u4), rhom1e); // +0- -0+
	feq[17] = fma(rhoe, fma(0.5f, fma(u5, u5, c3), u5), rhom1e); feq[18] = fma(rhoe, fma(0.5f, fma(u5, u5, c3), -u5), rhom1e); // 0+- 0-+
	feq[19] = fma(rhoc, fma(0.5f, fma(u6, u6, c3), u6), rhom1c); feq[20] = fma(rhoc, fma(0.5f, fma(u6, u6, c3), -u6), rhom1c); // +++ ---
	feq[21] = fma(rhoc, fma(0.5f, fma(u7, u7, c3), u7), rhom1c); feq[22] = fma(rhoc, fma(0.5f, fma(u7, u7, c3), -u7), rhom1c); // ++- --+
	feq[23] = fma(rhoc, fma(0.5f, fma(u8, u8, c3), u8), rhom1c); feq[24] = fma(rhoc, fma(0.5f, fma(u8, u8, c3), -u8), rhom1c); // +-+ -+-
	feq[25] = fma(rhoc, fma(0.5f, fma(u9, u9, c3), u9), rhom1c); feq[26] = fma(rhoc, fma(0.5f, fma(u9, u9, c3), -u9), rhom1c); // -++ +--
)+"#endif"+R( // D3Q27
} // calculate_f_eq()

)+R(void calculate_rho_u(const float* f, float* rhon, float* uxn, float* uyn, float* uzn) { // calculate density and velocity fields from fi
	float rho=f[0], ux, uy, uz;
	for(uint i=1u; i<def_velocity_set; i++) rho += f[i]; // calculate density from fi
	rho += 1.0f; // add 1.0f last to avoid digit extinction effects when summing up fi (perturbation method / DDF-shifting)
)+"#if defined(D2Q9)"+R(
	ux = f[1]-f[2]+f[5]-f[6]+f[7]-f[8]; // calculate velocity from fi (alternating + and - for best accuracy)
	uy = f[3]-f[4]+f[5]-f[6]+f[8]-f[7];
	uz = 0.0f;
)+"#elif defined(D3Q15)"+R(
	ux = f[ 1]-f[ 2]+f[ 7]-f[ 8]+f[ 9]-f[10]+f[11]-f[12]+f[14]-f[13]; // calculate velocity from fi (alternating + and - for best accuracy)
	uy = f[ 3]-f[ 4]+f[ 7]-f[ 8]+f[ 9]-f[10]+f[12]-f[11]+f[13]-f[14];
	uz = f[ 5]-f[ 6]+f[ 7]-f[ 8]+f[10]-f[ 9]+f[11]-f[12]+f[13]-f[14];
)+"#elif defined(D3Q19)"+R(
	ux = f[ 1]-f[ 2]+f[ 7]-f[ 8]+f[ 9]-f[10]+f[13]-f[14]+f[15]-f[16]; // calculate velocity from fi (alternating + and - for best accuracy)
	uy = f[ 3]-f[ 4]+f[ 7]-f[ 8]+f[11]-f[12]+f[14]-f[13]+f[17]-f[18];
	uz = f[ 5]-f[ 6]+f[ 9]-f[10]+f[11]-f[12]+f[16]-f[15]+f[18]-f[17];
)+"#elif defined(D3Q27)"+R(
	ux = f[ 1]-f[ 2]+f[ 7]-f[ 8]+f[ 9]-f[10]+f[13]-f[14]+f[15]-f[16]+f[19]-f[20]+f[21]-f[22]+f[23]-f[24]+f[26]-f[25]; // calculate velocity from fi (alternating + and - for best accuracy)
	uy = f[ 3]-f[ 4]+f[ 7]-f[ 8]+f[11]-f[12]+f[14]-f[13]+f[17]-f[18]+f[19]-f[20]+f[21]-f[22]+f[24]-f[23]+f[25]-f[26];
	uz = f[ 5]-f[ 6]+f[ 9]-f[10]+f[11]-f[12]+f[16]-f[15]+f[18]-f[17]+f[19]-f[20]+f[22]-f[21]+f[23]-f[24]+f[25]-f[26];
)+"#endif"+R( // D3Q27
	*rhon = rho;
	*uxn = ux/rho;
	*uyn = uy/rho;
	*uzn = uz/rho;
} // calculate_rho_u()

)+"#ifdef VOLUME_FORCE"+R(
)+R(void calculate_forcing_terms(const float ux, const float uy, const float uz, const float fx, const float fy, const float fz, float* Fin) { // calculate volume force terms Fin from velocity field (Guo forcing, Krueger p.233f)
)+"#ifdef D2Q9"+R(
	const float uF = -0.33333334f*fma(ux, fx, uy*fy); // 2D
)+"#else"+R( // D2Q9
	const float uF = -0.33333334f*fma(ux, fx, fma(uy, fy, uz*fz)); // 3D
)+"#endif"+R( // D2Q9
	Fin[0] = 9.0f*def_w0*uF ; // 000 (identical for all velocity sets)
	for(uint i=1u; i<def_velocity_set; i++) { // loop is entirely unrolled by compiler, no unnecessary FLOPs are happening
		Fin[i] = 9.0f*w(i)*fma(c(i)*fx+c(def_velocity_set+i)*fy+c(2u*def_velocity_set+i)*fz, c(i)*ux+c(def_velocity_set+i)*uy+c(2u*def_velocity_set+i)*uz+0.33333334f, uF);
	}
} // calculate_forcing_terms()
)+"#endif"+R( // VOLUME_FORCE

)+"#ifdef MOVING_BOUNDARIES"+R(
)+R(void apply_moving_boundaries(float* fhn, const uint* j, const global float* u, const global uchar* flags) { // apply Dirichlet velocity boundaries if necessary (Krueger p.180, rho_solid=1)
	uint ji; // reads velocities of only neighboring boundary nodes, which do not change during simulation
	for(uint i=1u; i<def_velocity_set; i+=2u) { // loop is entirely unrolled by compiler, no unnecessary memory access is happening
		const float w6 = -6.0f*w(i); // w(i) = w(i+1) if i is odd
		ji = j[i+1u]; fhn[i   ] = (flags[ji]&TYPE_BO)==TYPE_S ? fma(w6, c(i+1u)*u[ji]+c(def_velocity_set+i+1u)*u[def_N+(ulong)ji]+c(2u*def_velocity_set+i+1u)*u[2ul*def_N+(ulong)ji], fhn[i   ]) : fhn[i   ]; // boundary : regular
		ji = j[i   ]; fhn[i+1u] = (flags[ji]&TYPE_BO)==TYPE_S ? fma(w6, c(i   )*u[ji]+c(def_velocity_set+i   )*u[def_N+(ulong)ji]+c(2u*def_velocity_set+i   )*u[2ul*def_N+(ulong)ji], fhn[i+1u]) : fhn[i+1u];
	}
} // apply_moving_boundaries()
)+"#endif"+R( // MOVING_BOUNDARIES

)+"#ifdef SURFACE"+R(
)+R(void average_neighbors_non_gas(const uint n, const global float* rho, const global float* u, const global uchar* flags, float* rhon, float* uxn, float* uyn, float* uzn) { // calculate average density and velocity of neighbors of node n
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	float rhot=0.0f, uxt=0.0f, uyt=0.0f, uzt=0.0f, counter=0.0f; // average over all fluid/interface neighbors
	for(uint i=1u; i<def_velocity_set; i++) {
		const uchar flagsji_sus = flags[j[i]]&(TYPE_SU|TYPE_S); // extract SURFACE flags
		if(flagsji_sus==TYPE_F||flagsji_sus==TYPE_I||flagsji_sus==TYPE_IF) { // fluid or interface or (interface->fluid) neighbor
			counter += 1.0f;
			rhot += rho[               j[i]];
			uxt  += u[                 j[i]];
			uyt  += u[    def_N+(ulong)j[i]];
			uzt  += u[2ul*def_N+(ulong)j[i]];
		}
	}
	*rhon = counter>0.0f ? rhot/counter : 1.0f;
	*uxn  = counter>0.0f ? uxt /counter : 0.0f;
	*uyn  = counter>0.0f ? uyt /counter : 0.0f;
	*uzn  = counter>0.0f ? uzt /counter : 0.0f;
}
)+R(void average_neighbors_fluid(const uint n, const global float* rho, const global float* u, const global uchar* flags, float* rhon, float* uxn, float* uyn, float* uzn) { // calculate average density and velocity of neighbors of node n
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	float rhot=0.0f, uxt=0.0f, uyt=0.0f, uzt=0.0f, counter=0.0f; // average over all fluid/interface neighbors
	for(uint i=1u; i<def_velocity_set; i++) {
		const uchar flagsji_su = flags[j[i]]&TYPE_SU;
		if(flagsji_su==TYPE_F) { // fluid neighbor
			counter += 1.0f;
			rhot += rho[               j[i]];
			uxt  += u[                 j[i]];
			uyt  += u[    def_N+(ulong)j[i]];
			uzt  += u[2ul*def_N+(ulong)j[i]];
		}
	}
	*rhon = counter>0.0f ? rhot/counter : 1.0f;
	*uxn  = counter>0.0f ? uxt /counter : 0.0f;
	*uyn  = counter>0.0f ? uyt /counter : 0.0f;
	*uzn  = counter>0.0f ? uzt /counter : 0.0f;
}
)+R(float calculate_phi(const float rhon, const float massn, const uchar flagsn) { // calculate fill level
	return flagsn&TYPE_F ? 1.0f : flagsn&TYPE_I ? rhon>0.0f ? clamp(massn/rhon, 0.0f, 1.0f) : 0.5f : 0.0f;
}
)+R(float3 calculate_normal_py(const float* phij) { // calculate surface normal vector (Parker-youngs approximation, more accurate, works only for D3Q27 neighborhood)
	float3 n; // normal vector
)+"#ifdef D2Q9"+R(
	n.x = 2.0f*(phij[2]-phij[1])+phij[6]-phij[5]+phij[8]-phij[7];
	n.y = 2.0f*(phij[4]-phij[3])+phij[6]-phij[5]+phij[7]-phij[8];
	n.z = 0.0f;
)+"#else"+R( // D2Q9
	n.x = 4.0f*(phij[ 2]-phij[ 1])+2.0f*(phij[ 8]-phij[ 7]+phij[10]-phij[ 9]+phij[14]-phij[13]+phij[16]-phij[15])+phij[20]-phij[19]+phij[22]-phij[21]+phij[24]-phij[23]+phij[25]-phij[26];
	n.y = 4.0f*(phij[ 4]-phij[ 3])+2.0f*(phij[ 8]-phij[ 7]+phij[12]-phij[11]+phij[13]-phij[14]+phij[18]-phij[17])+phij[20]-phij[19]+phij[22]-phij[21]+phij[23]-phij[24]+phij[26]-phij[25];
	n.z = 4.0f*(phij[ 6]-phij[ 5])+2.0f*(phij[10]-phij[ 9]+phij[12]-phij[11]+phij[15]-phij[16]+phij[17]-phij[18])+phij[20]-phij[19]+phij[21]-phij[22]+phij[24]-phij[23]+phij[26]-phij[25];
)+"#endif"+R( // D2Q9
	return normalize(n);
}
)+R(float plic_cube_reduced(const float V, const float n1, const float n2, const float n3) { // optimized solution from SZ and Kawano, source: https://doi.org/10.3390/computation10020021
	const float n12=n1+n2, n3V=n3*V;
	if(n12<=2.0f*n3V) return n3V+0.5f*n12; // case (5)
	const float sqn1=sq(n1), n26=6.0f*n2, v1=sqn1/n26; // after case (5) check n2>0 is true
	if(v1<=n3V && n3V<v1+0.5f*(n2-n1)) return 0.5f*(n1+sqrt(sqn1+8.0f*n2*(n3V-v1))); // case (2)
	const float V6 = n1*n26*n3V;
	if(n3V<v1) return cbrt(V6); // case (1)
	const float v3 = n3<n12 ? (sq(n3)*(3.0f*n12-n3)+sqn1*(n1-3.0f*n3)+sq(n2)*(n2-3.0f*n3))/(n1*n26) : 0.5f*n12; // after case (2) check n1>0 is true
	const float sqn12=sqn1+sq(n2), V6cbn12=V6-cb(n1)-cb(n2);
	const bool case34 = n3V<v3; // true: case (3), false: case (4)
	const float a = case34 ? V6cbn12 : 0.5f*(V6cbn12-cb(n3));
	const float b = case34 ?   sqn12 : 0.5f*(sqn12+sq(n3));
	const float c = case34 ?     n12 : 0.5f;
	const float t = sqrt(sq(c)-b);
	return c-2.0f*t*sin(0.33333334f*asin((cb(c)-0.5f*a-1.5f*b*c)/cb(t)));
}
)+R(float plic_cube(const float V0, const float3 n) { // unit cube - plane intersection: volume V0 in [0,1], normal vector n -> plane offset d0
	const float ax=fabs(n.x), ay=fabs(n.y), az=fabs(n.z), V=0.5f-fabs(V0-0.5f), l=ax+ay+az; // eliminate symmetry cases, normalize n using L1 norm
	const float n1 = fmin(fmin(ax, ay), az)/l;
	const float n3 = fmax(fmax(ax, ay), az)/l;
	const float n2 = fdim(1.0f, n1+n3); // ensure n2>=0
	const float d = plic_cube_reduced(V, n1, n2, n3); // calculate PLIC with reduced symmetry
	return l*copysign(0.5f-d, V0-0.5f); // rescale result and apply symmetry for V0>0.5
}
)+R(void get_remaining_neighbor_phij(const uint n, const float* phit, const global float* phi, float* phij) { // get remaining phij for D3Q27 neighborhood
)+"#ifndef D3Q27"+R(
	uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
	calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
)+"#endif"+R( // D3Q27
)+"#if defined(D3Q15)"+R(
	uint j[12]; // calculate neighbor indices
	j[ 0] = xp+yp+z0; j[ 1] = xm+ym+z0; // ++0 --0
	j[ 2] = xp+y0+zp; j[ 3] = xm+y0+zm; // +0+ -0-
	j[ 4] = x0+yp+zp; j[ 5] = x0+ym+zm; // 0++ 0--
	j[ 6] = xp+ym+z0; j[ 7] = xm+yp+z0; // +-0 -+0
	j[ 8] = xp+y0+zm; j[ 9] = xm+y0+zp; // +0- -0+
	j[10] = x0+yp+zm; j[11] = x0+ym+zp; // 0+- 0-+
	for(uint i=0u; i<7u; i++) phij[i] = phit[i];
	for(uint i=7u; i<19u; i++) phij[i] = phi[j[i-7u]];
	for(uint i=19u; i<27u; i++) phij[i] = phit[i-12u];
)+"#elif defined(D3Q19)"+R(
	uint j[8]; // calculate remaining neighbor indices
	j[0] = xp+yp+zp; j[1] = xm+ym+zm; // +++ ---
	j[2] = xp+yp+zm; j[3] = xm+ym+zp; // ++- --+
	j[4] = xp+ym+zp; j[5] = xm+yp+zm; // +-+ -+-
	j[6] = xm+yp+zp; j[7] = xp+ym+zm; // -++ +--
	for(uint i=0u; i<19u; i++) phij[i] = phit[i];
	for(uint i=19u; i<27u; i++) phij[i] = phi[j[i-19u]];
)+"#elif defined(D3Q27)"+R(
	for(uint i=0u; i<def_velocity_set; i++) phij[i] = phit[i];
)+"#endif"+R( // D3Q27
}
)+R(float c_D3Q27(const uint i) { // avoid constant keyword by encapsulating data in function which gets inlined by compiler
	const float c[3*27] = {
		0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, // x
		0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1, // y
		0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1  // z
	};
	return c[i];
}
)+R(float curvature_calculation(const uint n, const float* phit, const global float* phi) { // calculate surface curvature, always use D3Q27 stencil here, source: https://doi.org/10.3390/computation10020021
)+"#ifndef D2Q9"+R(
	float phij[27];
	get_remaining_neighbor_phij(n, phit, phi, phij); // complete neighborhood from whatever velocity set is selected to D3Q27
	const float3 bz = calculate_normal_py(phij); // new coordinate system: bz is normal to surface, bx and by are tangent to surface
	const float3 rn = (float3)(0.56270900f, 0.32704452f, 0.75921047f); // random normalized vector that is just by random chance not collinear with bz
	const float3 by = normalize(cross(bz, rn)); // normalize() is necessary here because bz and rn are not perpendicular
	const float3 bx = cross(by, bz);
	uint number = 0; // number of neighboring interface points
	float3 p[24]; // number of neighboring interface points is less or equal than than 26 minus 1 gas and minus 1 fluid point = 24
	const float center_offset = plic_cube(phij[0], bz); // calculate z-offset PLIC of center point only once
	for(uint i=1u; i<27u; i++) { // iterate over neighbors, no loop unrolling here (50% better perfoemance without loop unrolling)
		if(phij[i]>0.0f&&phij[i]<1.0f) { // limit neighbors to interface nodes
			const float3 ei = (float3)(c_D3Q27(i), c_D3Q27(27u+i), c_D3Q27(2u*27u+i)); // assume neighbor normal vector is the same as center normal vector
			const float offset = plic_cube(phij[i], bz)-center_offset;
			p[number++] = (float3)(dot(ei, bx), dot(ei, by), dot(ei, bz)+offset); // do coordinate system transformation into (x, y, f(x,y)) and apply PLIC pffsets
		}
	}
	float M[25], x[5]={0.0f,0.0f,0.0f,0.0f,0.0f}, b[5]={0.0f,0.0f,0.0f,0.0f,0.0f};
	for(uint i=0u; i<25u; i++) M[i] = 0.0f;
	for(uint i=0u; i<number; i++) { // f(x,y)=A*x2+B*y2+C*x*y+H*x+I*y, x=(A,B,C,H,I), Q=(x2,y2,x*y,x,y), M*x=b, M=Q*Q^T, b=Q*z
		const float x=p[i].x, y=p[i].y, z=p[i].z, x2=x*x, y2=y*y, x3=x2*x, y3=y2*y;
		/**/M[ 0]+=x2*x2; M[ 1]+=x2*y2; M[ 2]+=x3*y ; M[ 3]+=x3   ; M[ 4]+=x2*y ; b[0]+=x2   *z;
		/*M[ 5]+=x2*y2;*/ M[ 6]+=y2*y2; M[ 7]+=x *y3; M[ 8]+=x *y2; M[ 9]+=   y3; b[1]+=   y2*z;
		/*M[10]+=x3*y ; M[11]+=x *y3;*/ M[12]+=x2*y2; M[13]+=x2*y ; M[14]+=x *y2; b[2]+=x *y *z;
		/*M[15]+=x3   ; M[16]+=x *y2; M[17]+=x2*y ;*/ M[18]+=x2   ; M[19]+=x *y ; b[3]+=x    *z;
		/*M[20]+=x2*y ; M[21]+=   y3; M[22]+=x *y2; M[23]+=x *y ;*/ M[24]+=   y2; b[4]+=   y *z;
	}
	for(uint i=1u; i<5u; i++) { // use symmetry of matrix to save arithmetic operations
		for(uint j=0u; j<i; j++) M[i*5+j] = M[j*5+i];
	}
	if(number>=5u) lu_solve(M, x, b, 5, 5);
	else lu_solve(M, x, b, 5, min(5u, number)); // cannot do loop unrolling here -> slower -> extra if-else to avoid slowdown
	const float A=x[0], B=x[1], C=x[2], H=x[3], I=x[4];
	const float K = (A*(I*I+1.0f)+B*(H*H+1.0f)-C*H*I)*cb(rsqrt(H*H+I*I+1.0f)); // mean curvature of Monge patch (x, y, f(x, y))
)+"#else"+R( // D2Q9
	const float3 by = calculate_normal_py(phit); // new coordinate system: bz is normal to surface, bx and by are tangent to surface
	const float3 bx = cross(by, (float3)(0.0f, 0.0f, 1.0f)); // normalize() is necessary here because bz and rn are not perpendicular
	uint number = 0u; // number of neighboring interface points
	float2 p[6]; // number of neighboring interface points is less or equal than than 8 minus 1 gas and minus 1 fluid point = 6
	const float center_offset = plic_cube(phit[0], by); // calculate z-offset PLIC of center point only once
	for(uint i=1u; i<9u; i++) { // iterate over neighbors, no loop unrolling here (50% better perfoemance without loop unrolling)
		if(phit[i]>0.0f&&phit[i]<1.0f) { // limit neighbors to interface nodes
			const float3 ei = (float3)(c_D2Q9(i), c_D2Q9(9+i), 0.0f); // assume neighbor normal vector is the same as center normal vector
			const float offset = plic_cube(phit[i], by)-center_offset;
			p[number++] = (float2)(dot(ei, bx), dot(ei, by)+offset); // do coordinate system transformation into (x, f(x)) and apply PLIC pffsets
		}
	}
	float M[4]={0.0f,0.0f,0.0f,0.0f}, x[2]={0.0f,0.0f}, b[2]={0.0f,0.0f};
	for(uint i=0u; i<number; i++) { // f(x,y)=A*x2+H*x, x=(A,H), Q=(x2,x), M*x=b, M=Q*Q^T, b=Q*z
		const float x=p[i].x, y=p[i].y, x2=x*x, x3=x2*x;
		/**/M[0]+=x2*x2; M[1]+=x3; b[0]+=x2*y;
		/*M[2]+=x3   ;*/ M[3]+=x2; b[1]+=x *y;
	}
	M[2] = M[1]; // use symmetry of matrix to save arithmetic operations
	if(number>=2u) lu_solve(M, x, b, 2, 2);
	else lu_solve(M, x, b, 2, min(2u, number)); // cannot do loop unrolling here -> slower -> extra if-else to avoid slowdown
	const float A=x[0], H=x[1];
	const float K = 2.0f*A*cb(rsqrt(H*H+1.0f)); // mean curvature of Monge patch (x, f(x)), note that curvature definition in 2D is different than 3D (additional factor 2)
)+"#endif"+R( // D2Q9
	return clamp(K, -1.0f, 1.0f); // prevent extreme pressures in the case of almost degenerate matrices
}
)+"#endif"+R( // SURFACE

)+"#ifdef TEMPERATURE"+R(
)+R(void neighbors_temperature(const uint n, uint* j7) { // calculate neighbor indices
	uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
	calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
	j7[0] = n;
	j7[1] = xp+y0+z0; j7[2] = xm+y0+z0; // +00 -00
	j7[3] = x0+yp+z0; j7[4] = x0+ym+z0; // 0+0 0-0
	j7[5] = x0+y0+zp; j7[6] = x0+y0+zm; // 00+ 00-
}
)+R(void calculate_g_eq(const float T, const float ux, const float uy, const float uz, float* geq) { // calculate g_equilibrium from density and velocity field (perturbation method / DDF-shifting)
	const float wsT4=0.5f*T, wsTm1=0.125f*(T-1.0f); // 0.125f*T*4.0f (straight directions in D3Q7), wsTm1 is arithmetic optimization to minimize digit extinction, lattice speed of sound is 1/2 for D3Q7 and not 1/sqrt(3)
	geq[0] = fma(0.25f, T, -0.25f); // 000
	geq[1] = fma(wsT4, ux, wsTm1); geq[2] = fma(wsT4, -ux, wsTm1); // +00 -00, source: http://dx.doi.org/10.1016/j.ijheatmasstransfer.2009.11.014
	geq[3] = fma(wsT4, uy, wsTm1); geq[4] = fma(wsT4, -uy, wsTm1); // 0+0 0-0
	geq[5] = fma(wsT4, uz, wsTm1); geq[6] = fma(wsT4, -uz, wsTm1); // 00+ 00-
}
)+R(void load_g(const uint n, float* ghn, const global fpxx* gi, const uint* j7, const ulong t) {
	ghn[0] = load(gi, index_f(n, 0u)); // Esoteric-Pull
	for(uint i=1u; i<7u; i+=2u) {
		ghn[i   ] = load(gi, index_f(n    , t%2ul ? i    : i+1u));
		ghn[i+1u] = load(gi, index_f(j7[i], t%2ul ? i+1u : i   ));
	}
}
)+R(void store_g(const uint n, const float* ghn, global fpxx* gi, const uint* j7, const ulong t) {
	store(gi, index_f(n, 0u), ghn[0]); // Esoteric-Pull
	for(uint i=1u; i<7u; i+=2u) {
		store(gi, index_f(j7[i], t%2ul ? i+1u : i   ), ghn[i   ]);
		store(gi, index_f(n    , t%2ul ? i    : i+1u), ghn[i+1u]);
	}
}
)+"#endif"+R( // TEMPERATURE

)+R(void load_f(const uint n, float* fhn, const global fpxx* fi, const uint* j, const ulong t) {
	fhn[0] = load(fi, index_f(n, 0u)); // Esoteric-Pull
	for(uint i=1u; i<def_velocity_set; i+=2u) {
		fhn[i   ] = load(fi, index_f(n   , t%2ul ? i    : i+1u));
		fhn[i+1u] = load(fi, index_f(j[i], t%2ul ? i+1u : i   ));
	}
}
)+R(void store_f(const uint n, const float* fhn, global fpxx* fi, const uint* j, const ulong t) {
	store(fi, index_f(n, 0u), fhn[0]); // Esoteric-Pull
	for(uint i=1u; i<def_velocity_set; i+=2u) {
		store(fi, index_f(j[i], t%2ul ? i+1u : i   ), fhn[i   ]);
		store(fi, index_f(n   , t%2ul ? i    : i+1u), fhn[i+1u]);
	}
}

)+"#ifdef SURFACE"+R(
)+R(void load_f_outgoing(const uint n, float* fon, const global fpxx* fi, const uint* j, const ulong t) { // load outgoing DDFs, even: 1:1 like stream-out odd, odd: 1:1 like stream-out even
	for(uint i=1u; i<def_velocity_set; i+=2u) { // Esoteric-Pull
		fon[i   ] = load(fi, index_f(j[i], t%2ul ? i    : i+1u));
		fon[i+1u] = load(fi, index_f(n   , t%2ul ? i+1u : i   ));
	}
}
)+R(void store_f_reconstructed(const uint n, const float* fhn, global fpxx* fi, const uint* j, const ulong t, const uchar* flagsj_su) { // store reconstructed gas DDFs, even: 1:1 like stream-in even, odd: 1:1 like stream-in odd
	for(uint i=1u; i<def_velocity_set; i+=2u) { // Esoteric-Pull
		if(flagsj_su[i+1u]==TYPE_G) store(fi, index_f(n   , t%2ul ? i    : i+1u), fhn[i   ]); // only store reconstructed gas DDFs to locations from which
		if(flagsj_su[i   ]==TYPE_G) store(fi, index_f(j[i], t%2ul ? i+1u : i   ), fhn[i+1u]); // they are going to be streamed in during next stream_collide()
	}
}
)+"#endif"+R( // SURFACE



)+R(kernel void initialize)+"("+R(global fpxx* fi, const global float* rho, global float* u, global uchar* flags // ) { // initialize LBM
)+"#ifdef SURFACE"+R(
	, global float* mass, global float* massex, global float* phi // argument order is important
)+"#endif"+R( // SURFACE
)+"#ifdef TEMPERATURE"+R(
	, global fpxx* gi, const global float* T // argument order is important
)+"#endif"+R( // TEMPERATURE
)+") {"+R( // initialize()
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute initialize() on halo
	uchar flagsn = flags[n];
	const uchar flagsn_bo = flagsn&TYPE_BO; // extract boundary flags
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	uchar flagsj[def_velocity_set]; // cache neighbor flags for multiple readings
	for(uint i=1u; i<def_velocity_set; i++) flagsj[i] = flags[j[i]];
	if(flagsn_bo==TYPE_S) { // node is solid
		bool TYPE_ONLY_S = true; // has only solid neighbors
		for(uint i=1u; i<def_velocity_set; i++) TYPE_ONLY_S = TYPE_ONLY_S&&(flagsj[i]&TYPE_BO)==TYPE_S;
		if(TYPE_ONLY_S) {
			u[                 n] = 0.0f; // reset velocity for solid lattice points with only boundary neighbors
			u[    def_N+(ulong)n] = 0.0f;
			u[2ul*def_N+(ulong)n] = 0.0f;
		}
)+"#ifndef MOVING_BOUNDARIES"+R(
		if(flagsn_bo==TYPE_S) {
			u[                 n] = 0.0f; // reset velocity for all solid lattice points
			u[    def_N+(ulong)n] = 0.0f;
			u[2ul*def_N+(ulong)n] = 0.0f;
		}
)+"#else"+R( // MOVING_BOUNDARIES
	} else if(flagsn_bo!=TYPE_E) { // local lattice point is not solid and not equilibrium boundary
		bool next_to_moving_boundary = false;
		for(uint i=1u; i<def_velocity_set; i++) {
			next_to_moving_boundary = next_to_moving_boundary||((flagsj[i]&TYPE_BO)==TYPE_S&&(u[j[i]]!=0.0f||u[def_N+(ulong)j[i]]!=0.0f||u[2ul*def_N+(ulong)j[i]]!=0.0f));
		}
		flags[n] = flagsn = next_to_moving_boundary ? flagsn|TYPE_MS : flagsn&~TYPE_MS; // mark/unmark nodes next to TYPE_S nodes with velocity!=0 with TYPE_MS
)+"#endif"+R( // MOVING_BOUNDARIES
	}
	float feq[def_velocity_set]; // f_equilibrium
	calculate_f_eq(rho[n], u[n], u[def_N+(ulong)n], u[2ul*def_N+(ulong)n], feq);
)+"#ifdef SURFACE"+R( // automatically generate the interface layer between fluid and gas
	{ // separate block to avoid variable name conflicts
		float phin = phi[n];
		if(!(flagsn&(TYPE_S|TYPE_E|TYPE_T|TYPE_F|TYPE_I))) flagsn = (flagsn&~TYPE_SU)|TYPE_G; // change all non-fluid and non-interface flags to gas
		if((flagsn&TYPE_SU)==TYPE_G) { // node with updated flags is gas
			bool change = false; // check if node has to be changed to interface
			for(uint i=1u; i<def_velocity_set; i++) change = change||(flagsj[i]&TYPE_SU)==TYPE_F; // if neighbor flag fluid is set, the node must be interface
			if(change) { // create interface automatically if phi has not explicitely defined for the interface layer
				flagsn = (flagsn&~TYPE_SU)|TYPE_I; // node must be interface
				phin = 0.5f;
				float rhon, uxn, uyn, uzn; // initialize interface nodes with average density/velocity of fluid neighbors
				average_neighbors_fluid(n, rho, u, flags, &rhon, &uxn, &uyn, &uzn); // get average rho/u from all fluid neighbors
				calculate_f_eq(rhon, uxn, uyn, uzn, feq); // calculate equilibrium DDFs
			}
		}
		if((flagsn&TYPE_SU)==TYPE_G) { // node with updated flags is still gas
			u[                 n] = 0.0f; // reset velocity for gas nodes
			u[    def_N+(ulong)n] = 0.0f;
			u[2ul*def_N+(ulong)n] = 0.0f;
			phin = 0.0f;
		} else if((flagsn&TYPE_SU)==TYPE_I && (phin<0.0f||phin>1.0f)) {
			phin = 0.5f; // node should be interface, but phi was invalid
		} else if((flagsn&TYPE_SU)==TYPE_F) {
			phin = 1.0f;
		}
		phi[n] = phin;
		mass[n] = phin*rho[n];
		massex[n] = 0.0f; // reset excess mass
		flags[n] = flagsn;
	}
)+"#endif"+R( // SURFACE
)+"#ifdef TEMPERATURE"+R(
	{ // separate block to avoid variable name conflicts
		float geq[7];
		calculate_g_eq(T[n], u[n], u[def_N+(ulong)n], u[2ul*def_N+(ulong)n], geq);
		uint j7[7]; // neighbors of D3Q7 subset
		neighbors_temperature(n, j7);
		store_g(n, geq, gi, j7, 1ul);
	}
)+"#endif"+R( // TEMPERATURE
	store_f(n, feq, fi, j, 1ul); // write to fi
} // initialize()

)+"#ifdef MOVING_BOUNDARIES"+R(
)+R(kernel void update_moving_boundaries(const global float* u, global uchar* flags) { // mark/unmark nodes next to TYPE_S nodes with velocity!=0 with TYPE_MS
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute update_moving_boundaries() on halo
	const uchar flagsn = flags[n];
	const uchar flagsn_bo = flagsn&TYPE_BO; // extract boundary flags
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	uchar flagsj[def_velocity_set]; // cache neighbor flags for multiple readings
	for(uint i=1u; i<def_velocity_set; i++) flagsj[i] = flags[j[i]];
	if(flagsn_bo!=TYPE_S&&flagsn_bo!=TYPE_E&&!(flagsn&TYPE_T)) { // local lattice point is not solid and not equilibrium boundary and not temperature boundary
		bool next_to_moving_boundary = false;
		for(uint i=1u; i<def_velocity_set; i++) {
			next_to_moving_boundary = next_to_moving_boundary||((u[j[i]]!=0.0f||u[def_N+(ulong)j[i]]!=0.0f||u[2ul*def_N+(ulong)j[i]]!=0.0f)&&(flagsj[i]&TYPE_BO)==TYPE_S);
		}
		flags[n] = next_to_moving_boundary ? flagsn|TYPE_MS : flagsn&~TYPE_MS; // mark/unmark nodes next to TYPE_S nodes with velocity!=0 with TYPE_MS
	}
} // update_moving_boundaries()
)+"#endif"+R( // MOVING_BOUNDARIES



)+R(kernel void stream_collide)+"("+R(global fpxx* fi, global float* rho, global float* u, global uchar* flags, const ulong t, const float fx, const float fy, const float fz // ) { // main LBM kernel
)+"#ifdef FORCE_FIELD"+R(
	, const global float* F // argument order is important
)+"#endif"+R( // FORCE_FIELD
)+"#ifdef SURFACE"+R(
	, const global float* mass // argument order is important
)+"#endif"+R( // SURFACE
)+"#ifdef TEMPERATURE"+R(
	, global fpxx* gi, global float* T // argument order is important
)+"#endif"+R( // TEMPERATURE
)+") {"+R( // stream_collide()
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute stream_collide() on halo
	const uchar flagsn = flags[n]; // cache flags[n] for multiple readings
	const uchar flagsn_bo=flagsn&TYPE_BO, flagsn_su=flagsn&TYPE_SU; // extract boundary and surface flags
	if(flagsn_bo==TYPE_S||flagsn_su==TYPE_G) return; // if node is solid boundary or gas, just return

	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices

	float fhn[def_velocity_set]; // local DDFs
	load_f(n, fhn, fi, j, t); // perform streaming (part 2)

)+"#ifdef MOVING_BOUNDARIES"+R(
	if(flagsn_bo==TYPE_MS) apply_moving_boundaries(fhn, j, u, flags); // apply Dirichlet velocity boundaries if necessary (reads velocities of only neighboring boundary nodes, which do not change during simulation)
)+"#endif"+R( // MOVING_BOUNDARIES

	float rhon, uxn, uyn, uzn; // calculate local density and velocity for collision
)+"#ifndef EQUILIBRIUM_BOUNDARIES"+R(
	calculate_rho_u(fhn, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fi
)+"#else"+R( // EQUILIBRIUM_BOUNDARIES
	if(flagsn_bo==TYPE_E) {
		rhon = rho[               n]; // apply preset velocity/density
		uxn  = u[                 n];
		uyn  = u[    def_N+(ulong)n];
		uzn  = u[2ul*def_N+(ulong)n];
	} else {
		calculate_rho_u(fhn, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fi
	}
)+"#endif"+R( // EQUILIBRIUM_BOUNDARIES
	float fxn=fx, fyn=fy, fzn=fz; // force starts as constant volume force, can be modified before call of calculate_forcing_terms(...)
	float Fin[def_velocity_set]; // forcing terms

)+"#ifdef FORCE_FIELD"+R(
	{ // separate block to avoid variable name conflicts
		fxn += F[                 n]; // apply force field
		fyn += F[    def_N+(ulong)n];
		fzn += F[2ul*def_N+(ulong)n];
	}
)+"#endif"+R( // FORCE_FIELD

)+"#ifdef SURFACE"+R(
	if(flagsn_su==TYPE_I) { // node was interface, eventually initiate flag change
		bool TYPE_NO_F=true, TYPE_NO_G=true; // temporary flags for no fluid or gas neighbors
		for(uint i=1u; i<def_velocity_set; i++) {
			const uchar flagsji_su = flags[j[i]]&TYPE_SU; // extract SURFACE flags
			TYPE_NO_F = TYPE_NO_F&&flagsji_su!=TYPE_F;
			TYPE_NO_G = TYPE_NO_G&&flagsji_su!=TYPE_G;
		}
		const float massn = mass[n]; // load mass
		     if(massn>rhon || TYPE_NO_G) flags[n] = (flagsn&~TYPE_SU)|TYPE_IF; // set flag interface->fluid
		else if(massn<0.0f || TYPE_NO_F) flags[n] = (flagsn&~TYPE_SU)|TYPE_IG; // set flag interface->gas
	}
)+"#endif"+R( // SURFACE

)+"#ifdef TEMPERATURE"+R(
	{ // separate block to avoid variable name conflicts
		uint j7[7]; // neighbors of D3Q7 subset
		neighbors_temperature(n, j7);
		float ghn[7]; // read from gA and stream to gh (D3Q7 subset, periodic boundary conditions)
		load_g(n, ghn, gi, j7, t); // perform streaming (part 2)
		float Tn;
		if(flagsn&TYPE_T) {
			Tn = T[n]; // apply preset temperature
		} else {
			Tn = 0.0f;
			for(uint i=0u; i<7u; i++) Tn += ghn[i]; // calculate temperature from g
			Tn += 1.0f; // add 1.0f last to avoid digit extinction effects when summing up gi (perturbation method / DDF-shifting)
		}
		float geq[7]; // cache f_equilibrium[n]
		calculate_g_eq(Tn, uxn, uyn, uzn, geq); // calculate equilibrium DDFs
		if(flagsn&TYPE_T) {
			for(uint i=0u; i<7u; i++) ghn[i] = geq[i]; // just write geq to ghn (no collision)
		} else {
)+"#ifdef UPDATE_FIELDS"+R(
			T[n] = Tn; // update temperature field
)+"#endif"+R( // UPDATE_FIELDS
			for(uint i=0u; i<7u; i++) ghn[i] = fma(1.0f-def_w_T, ghn[i], def_w_T*geq[i]); // perform collision
		}
		store_g(n, ghn, gi, j7, t); // perform streaming (part 1)
		fxn -= fx*def_beta*(Tn-def_T_avg);
		fyn -= fy*def_beta*(Tn-def_T_avg);
		fzn -= fz*def_beta*(Tn-def_T_avg);
	}
)+"#endif"+R( // TEMPERATURE

	{ // separate block to avoid variable name conflicts
)+"#ifdef VOLUME_FORCE"+R( // apply force and collision operator, write to fi in video memory
		const float rho2 = 0.5f/rhon; // apply external volume force (Guo forcing, Krueger p.233f)
		uxn = clamp(fma(fxn, rho2, uxn), -def_c, def_c); // limit velocity (for stability purposes)
		uyn = clamp(fma(fyn, rho2, uyn), -def_c, def_c); // force term: F*dt/(2*rho)
		uzn = clamp(fma(fzn, rho2, uzn), -def_c, def_c);
		calculate_forcing_terms(uxn, uyn, uzn, fxn, fyn, fzn, Fin); // calculate volume force terms Fin from velocity field (Guo forcing, Krueger p.233f)
)+"#else"+R( // VOLUME_FORCE
		uxn = clamp(uxn, -def_c, def_c); // limit velocity (for stability purposes)
		uyn = clamp(uyn, -def_c, def_c); // force term: F*dt/(2*rho)
		uzn = clamp(uzn, -def_c, def_c);
		for(uint i=0u; i<def_velocity_set; i++) Fin[i] = 0.0f;
)+"#endif"+R( // VOLUME_FORCE
	}

)+"#ifndef EQUILIBRIUM_BOUNDARIES"+R(
)+"#ifdef UPDATE_FIELDS"+R(
	rho[               n] = rhon; // update density field
	u[                 n] = uxn; // update velocity field
	u[    def_N+(ulong)n] = uyn;
	u[2ul*def_N+(ulong)n] = uzn;
)+"#endif"+R( // UPDATE_FIELDS
)+"#else"+R( // EQUILIBRIUM_BOUNDARIES
)+"#ifdef UPDATE_FIELDS"+R(
	if(flagsn_bo!=TYPE_E) { // only update fields for non-TYPE_E nodes
		rho[               n] = rhon; // update density field
		u[                 n] = uxn; // update velocity field
		u[    def_N+(ulong)n] = uyn;
		u[2ul*def_N+(ulong)n] = uzn;
	}
)+"#endif"+R( // UPDATE_FIELDS
)+"#endif"+R( // EQUILIBRIUM_BOUNDARIES

	float feq[def_velocity_set]; // equilibrium DDFs
	calculate_f_eq(rhon, uxn, uyn, uzn, feq); // calculate equilibrium DDFs
	float w = def_w; // LBM relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)

)+"#ifdef SUBGRID"+R(
	{ // Smagorinsky-Lilly subgrid turbulence model, source: https://arxiv.org/pdf/comp-gas/9401004.pdf, in the eq. below (26), it is "tau_0" not "nu_0", and "sqrt(2)/rho" (they call "rho" "n") is missing
		const float tau0 = 1.0f/w; // source 2: https://youtu.be/V8ydRrdCzl0
		float Hxx=0.0f, Hyy=0.0f, Hzz=0.0f, Hxy=0.0f, Hxz=0.0f, Hyz=0.0f; // non-equilibrium stress tensor
		for(uint i=1u; i<def_velocity_set; i++) {
			const float fneqi = fhn[i]-feq[i];
			const float cxi=c(i), cyi=c(def_velocity_set+i), czi=c(2u*def_velocity_set+i);
			Hxx += cxi*cxi*fneqi; //Hyx += cyi*cxi*fneqi; Hzx += czi*cxi*fneqi; // symmetric tensor
			Hxy += cxi*cyi*fneqi; Hyy += cyi*cyi*fneqi; //Hzy += czi*cyi*fneqi;
			Hxz += cxi*czi*fneqi; Hyz += cyi*czi*fneqi; Hzz += czi*czi*fneqi;
		}
		const float Q = sq(Hxx)+sq(Hyy)+sq(Hzz)+2.0f*(sq(Hxy)+sq(Hxz)+sq(Hyz)); // Q = H*H, turbulent eddy viscosity nut = (C*Delta)^2*|S|, intensity of local strain rate tensor |S|=sqrt(2*S*S)
		w = 2.0f/(tau0+sqrt(sq(tau0)+0.76421222f*sqrt(Q)/rhon)); // 0.76421222 = 18*sqrt(2)*(C*Delta)^2, C = 1/pi*(2/(3*CK))^(3/4) = Smagorinsky-Lilly constant, CK = 3/2 = Kolmogorov constant, Delta = 1 = lattice constant
	} // modity LBM relaxation rate by increasing effective viscosity in regions of high strain rate (add turbulent eddy viscosity), nu_eff = nu_0+nu_t
)+"#endif"+R( // SUBGRID

)+"#if defined(SRT)"+R(
)+"#ifdef VOLUME_FORCE"+R(
	const float c_tau = fma(w, -0.5f, 1.0f);
	for(uint i=0u; i<def_velocity_set; i++) Fin[i] *= c_tau;
)+"#endif"+R( // VOLUME_FORCE
)+"#ifndef EQUILIBRIUM_BOUNDARIES"+R(
	for(uint i=0u; i<def_velocity_set; i++) fhn[i] = fma(1.0f-w, fhn[i], fma(w, feq[i], Fin[i])); // perform collision (SRT)
)+"#else"+R( // EQUILIBRIUM_BOUNDARIES
	for(uint i=0u; i<def_velocity_set; i++) fhn[i] = flagsn_bo==TYPE_E ? feq[i] : fma(1.0f-w, fhn[i], fma(w, feq[i], Fin[i])); // perform collision (SRT)
)+"#endif"+R( // EQUILIBRIUM_BOUNDARIES
)+"#elif defined(TRT)"+R(
	const float wp = w; // TRT: inverse of "+" relaxation time
	const float wm = 1.0f/(0.1875f/(1.0f/w-0.5f)+0.5f); // TRT: inverse of "-" relaxation time wm = 1.0f/(0.1875f/(3.0f*nu)+0.5f), nu = (1.0f/w-0.5f)/3.0f;
)+"#ifdef VOLUME_FORCE"+R(
	const float c_taup=fma(wp, -0.25f, 0.5f), c_taum=fma(wm, -0.25f, 0.5f); // source: https://arxiv.org/pdf/1901.08766.pdf
	float Fib[def_velocity_set]; // F_bar
	Fib[0] = Fin[0];
	for(uint i=1u; i<def_velocity_set; i+=2u) {
		Fib[i   ] = Fin[i+1u];
		Fib[i+1u] = Fin[i   ];
	}
	for(uint i=0u; i<def_velocity_set; i++) Fin[i] = fma(c_taup, Fin[i]+Fib[i], c_taum*(Fin[i]-Fib[i]));
)+"#endif"+R( // VOLUME_FORCE
	float fhb[def_velocity_set]; // fhn in inverse directions
	float feb[def_velocity_set]; // feq in inverse directions
	fhb[0] = fhn[0];
	feb[0] = feq[0];
	for(uint i=1u; i<def_velocity_set; i+=2u) {
		fhb[i   ] = fhn[i+1u];
		fhb[i+1u] = fhn[i   ];
		feb[i   ] = feq[i+1u];
		feb[i+1u] = feq[i   ];
	}
)+"#ifndef EQUILIBRIUM_BOUNDARIES"+R(
	for(uint i=0u; i<def_velocity_set; i++) fhn[i] = fma(0.5f*wp, feq[i]-fhn[i]+feb[i]-fhb[i], fma(0.5f*wm, feq[i]-feb[i]-fhn[i]+fhb[i], fhn[i]+Fin[i])); // perform collision (TRT)
)+"#else"+R( // EQUILIBRIUM_BOUNDARIES
	for(uint i=0u; i<def_velocity_set; i++) fhn[i] = flagsn_bo==TYPE_E ? feq[i] : fma(0.5f*wp, feq[i]-fhn[i]+feb[i]-fhb[i], fma(0.5f*wm, feq[i]-feb[i]-fhn[i]+fhb[i], fhn[i]+Fin[i])); // perform collision (TRT)
)+"#endif"+R( // EQUILIBRIUM_BOUNDARIES
)+"#endif"+R( // TRT

	store_f(n, fhn, fi, j, t); // perform streaming (part 1)
} // stream_collide()

)+"#ifdef SURFACE"+R(
)+R(kernel void surface_0(global fpxx* fi, const global float* rho, const global float* u, const global uchar* flags, global float* mass, const global float* massex, const global float* phi, const ulong t, const float fx, const float fy, const float fz) { // capture outgoing DDFs before streaming
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute surface_0() on halo
	const uchar flagsn = flags[n]; // cache flags[n] for multiple readings
	const uchar flagsn_bo=flagsn&TYPE_BO, flagsn_su=flagsn&TYPE_SU; // extract boundary and surface flags
	if(flagsn_bo==TYPE_S||flagsn_su==TYPE_G) return; // node processed here is fluid or interface

	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	float fhn[def_velocity_set]; // incoming DDFs
	load_f(n, fhn, fi, j, t); // load incoming DDFs
	float fon[def_velocity_set]; // outgoing DDFs
	fon[0] = fhn[0]; // fon[0] is already loaded in fhn[0]
	load_f_outgoing(n, fon, fi, j, t); // load outgoing DDFs

	float massn = mass[n];
	for(uint i=1u; i<def_velocity_set; i++) {
		massn += massex[j[i]]; // distribute excess mass from last step which is stored in neighbors
	}
	if(flagsn_su==TYPE_F) { // node is fluid
		for(uint i=1u; i<def_velocity_set; i++) massn += fhn[i]-fon[i]; // neighbor is fluid or interface node
	} else if(flagsn_su==TYPE_I) { // node is interface
		float phij[def_velocity_set]; // cache fill level of neighbor lattice points
		for(uint i=1u; i<def_velocity_set; i++) phij[i] = phi[j[i]]; // cache fill level of neighbor lattice points
		float rhon, uxn, uyn, uzn, rho_laplace=0.0f; // no surface tension if rho_laplace is not overwritten later
)+"#ifndef EQUILIBRIUM_BOUNDARIES"+R(
		calculate_rho_u(fon, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fon (not fhn)
)+"#else"+R( // EQUILIBRIUM_BOUNDARIES
		if(flagsn_bo==TYPE_E) {
			rhon = rho[               n]; // apply preset velocity/density
			uxn  = u[                 n];
			uyn  = u[    def_N+(ulong)n];
			uzn  = u[2ul*def_N+(ulong)n];
		} else {
			calculate_rho_u(fon, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fon (not fhn)
		}
)+"#endif"+R( // EQUILIBRIUM_BOUNDARIES
		uxn = clamp(uxn, -def_c, def_c); // limit velocity (for stability purposes)
		uyn = clamp(uyn, -def_c, def_c);
		uzn = clamp(uzn, -def_c, def_c);
		phij[0] = calculate_phi(rhon, massn, flagsn); // don't load phi[n] from memory, instead recalculate it with mass corrected by excess mass
		rho_laplace = def_6_sigma==0.0f ? 0.0f : def_6_sigma*curvature_calculation(n, phij, phi); // surface tension least squares fit (PLIC, most accurate)
		float feg[def_velocity_set]; // reconstruct f from neighbor gas lattice points
		const float rho2tmp = 0.5f/rhon; // apply external volume force (Guo forcing, Krueger p.233f)
		const float uxntmp = clamp(fma(fx, rho2tmp, uxn), -def_c, def_c); // limit velocity (for stability purposes)
		const float uyntmp = clamp(fma(fy, rho2tmp, uyn), -def_c, def_c); // force term: F*dt/(2*rho)
		const float uzntmp = clamp(fma(fz, rho2tmp, uzn), -def_c, def_c);
		calculate_f_eq(1.0f-rho_laplace, uxntmp, uyntmp, uzntmp, feg); // calculate gas equilibrium DDFs with constant ambient pressure
		uchar flagsj_su[def_velocity_set]; // cache neighbor flags for multiple readings
		for(uint i=1u; i<def_velocity_set; i++) flagsj_su[i] = flags[j[i]]&TYPE_SU;
		for(uint i=1u; i<def_velocity_set; i+=2u) { // calculate mass exchange between current node and fluid/interface nodes
			massn += flagsj_su[i   ]&(TYPE_F|TYPE_I) ? flagsj_su[i   ]==TYPE_F ? fhn[i+1]-fon[i   ] : 0.5f*(phij[i   ]+phij[0])*(fhn[i+1 ]-fon[i   ]) : 0.0f; // neighbor is fluid or interface node
			massn += flagsj_su[i+1u]&(TYPE_F|TYPE_I) ? flagsj_su[i+1u]==TYPE_F ? fhn[i  ]-fon[i+1u] : 0.5f*(phij[i+1u]+phij[0])*(fhn[i   ]-fon[i+1u]) : 0.0f; // fluid : interface : gas
		}
		for(uint i=1u; i<def_velocity_set; i+=2u) { // calculate reconstructed gas DDFs
			fhn[i   ] = feg[i+1u]-fon[i+1u]+feg[i   ];
			fhn[i+1u] = feg[i   ]-fon[i   ]+feg[i+1u];
		}
		store_f_reconstructed(n, fhn, fi, j, t, flagsj_su); // store reconstructed gas DDFs that are streamed in during the following stream_collide()
	}
	mass[n] = massn;
}
)+R(kernel void surface_1(global uchar* flags) { // prevent neighbors from interface->fluid nodes to become/be gas nodes
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N) return; // execute surface_1() also on halo
	const uchar flagsn_sus = flags[n]&(TYPE_SU|TYPE_S); // extract SURFACE flags
	if(flagsn_sus==TYPE_IF) { // flag interface->fluid is set
		uint j[def_velocity_set]; // neighbor indices
		neighbors(n, j); // calculate neighbor indices
		for(uint i=1u; i<def_velocity_set; i++) {
			const uchar flagsji = flags[j[i]];
			const uchar flagsji_su = flagsji&(TYPE_SU|TYPE_S); // extract SURFACE flags
			const uchar flagsji_r = flagsji&~TYPE_SU; // extract all non-SURFACE flags
			if(flagsji_su==TYPE_IG) flags[j[i]] = flagsji_r|TYPE_I; // prevent interface neighbor nodes from becoming gas
			else if(flagsji_su==TYPE_G) flags[j[i]] = flagsji_r|TYPE_GI; // neighbor node was gas and must change to interface
		}
	}
} // possible types at the end of surface_1(): TYPE_F / TYPE_I / TYPE_G / TYPE_IF / TYPE_IG / TYPE_GI
)+R(kernel void surface_2(global fpxx* fi, const global float* rho, const global float* u, global uchar* flags, const ulong t) {  // apply flag changes and calculate excess mass
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N) return; // execute surface_2() also on halo
	const uchar flagsn_sus = flags[n]&(TYPE_SU|TYPE_S); // extract SURFACE flags
	if(flagsn_sus==TYPE_GI) { // initialize the fi of gas nodes that should become interface
		float rhon, uxn, uyn, uzn; // average over all fluid/interface neighbors
		average_neighbors_non_gas(n, rho, u, flags, &rhon, &uxn, &uyn, &uzn); // get average rho/u from all fluid/interface neighbors
		float feq[def_velocity_set];
		calculate_f_eq(rhon, uxn, uyn, uzn, feq); // calculate equilibrium DDFs
		uint j[def_velocity_set];
		neighbors(n, j);
		store_f(n, feq, fi, j, t); // write feq to fi in video memory
	} else if(flagsn_sus==TYPE_IG) { // flag interface->gas is set
		uint j[def_velocity_set]; // neighbor indices
		neighbors(n, j); // calculate neighbor indices
		for(uint i=1u; i<def_velocity_set; i++) {
			const uchar flagsji = flags[j[i]];
			const uchar flagsji_su = flagsji&(TYPE_SU|TYPE_S); // extract SURFACE flags
			const uchar flagsji_r = flagsji&~TYPE_SU; // extract all non-SURFACE flags
			if(flagsji_su==TYPE_F||flagsji_su==TYPE_IF) {
				flags[j[i]] = flagsji_r|TYPE_I; // prevent fluid or interface neighbors that turn to fluid from being/becoming fluid
			}
		}
	}
} // possible types at the end of surface_2(): TYPE_F / TYPE_I / TYPE_G / TYPE_IF / TYPE_IG / TYPE_GI
)+R(kernel void surface_3(const global float* rho, global uchar* flags, global float* mass, global float* massex, global float* phi) { // apply flag changes and calculate excess mass
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute surface_3() on halo
	const uchar flagsn_sus = flags[n]&(TYPE_SU|TYPE_S); // extract SURFACE flags
	if(flagsn_sus&TYPE_S) return;
	const float rhon = rho[n]; // density of node n
	float massn = mass[n]; // mass of node n
	float massexn = 0.0f; // excess mass of node n
	float phin = 0.0f;
	if(flagsn_sus==TYPE_F) { // regular fluid node
		massexn = massn-rhon; // dump mass-rho difference into excess mass
		massn = rhon; // fluid node mass has to equal rho
		phin = 1.0f;
	} else if(flagsn_sus==TYPE_I) { // regular interface node
		massexn = massn>rhon ? massn-rhon : massn<0.0f ? massn : 0.0f; // allow interface nodes with mass>rho or mass<0
		massn = clamp(massn, 0.0f, rhon);
		phin = calculate_phi(rhon, massn, TYPE_I); // calculate fill level for next step (only necessary for interface nodes)
	} else if(flagsn_sus==TYPE_G) { // regular gas node
		massexn = massn; // dump remaining mass into excess mass
		massn = 0.0f;
		phin = 0.0f;
	} else if(flagsn_sus==TYPE_IF) { // flag interface->fluid is set
		flags[n] = (flags[n]&~TYPE_SU)|TYPE_F; // node becomes fluid
		massexn = massn-rhon; // dump mass-rho difference into excess mass
		massn = rhon; // fluid node mass has to equal rho
		phin = 1.0f; // set phi[n] to 1.0f for fluid nodes
	} else if(flagsn_sus==TYPE_IG) { // flag interface->gas is set
		flags[n] = (flags[n]&~TYPE_SU)|TYPE_G; // node becomes gas
		massexn = massn; // dump remaining mass into excess mass
		massn = 0.0f; // gas mass has to be zero
		phin = 0.0f; // set phi[n] to 0.0f for gas nodes
	} else if(flagsn_sus==TYPE_GI) { // flag gas->interface is set
		flags[n] = (flags[n]&~TYPE_SU)|TYPE_I; // node becomes interface
		massexn = massn>rhon ? massn-rhon : massn<0.0f ? massn : 0.0f; // allow interface nodes with mass>rho or mass<0
		massn = clamp(massn, 0.0f, rhon);
		phin = calculate_phi(rhon, massn, TYPE_I); // calculate fill level for next step (only necessary for interface nodes)
	}
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	uint counter = 0u; // count (fluid|interface) neighbors
	for(uint i=1u; i<def_velocity_set; i++) { // simple model: distribute excess mass equally to all interface and fluid neighbors
		const uchar flagsji_su = flags[j[i]]&(TYPE_SU|TYPE_S); // extract SURFACE flags
		counter += (uint)(flagsji_su==TYPE_F||flagsji_su==TYPE_I||flagsji_su==TYPE_IF||flagsji_su==TYPE_GI); // avoid branching
	}
	massn += counter>0u ? 0.0f : massexn; // if excess mass can't be distributed to neighboring interface or fluid nodes, add it to local mass (ensure mass conservation)
	massexn = counter>0u ? massexn/(float)counter : 0.0f; // divide excess mass up for all interface or fluid neighbors
	mass[n] = massn; // update mass
	massex[n] = massexn; // update excess mass
	phi[n] = phin; // update phi
} // possible types at the end of surface_3(): TYPE_F / TYPE_I / TYPE_G
)+"#endif"+R( // SURFACE

)+R(kernel void update_fields)+"("+R(const global fpxx* fi, global float* rho, global float* u, const global uchar* flags, const ulong t, const float fx, const float fy, const float fz // ) { // calculate fields from DDFs
)+"#ifdef FORCE_FIELD"+R(
	, const global float* F // argument order is important
)+"#endif"+R( // FORCE_FIELD
)+"#ifdef TEMPERATURE"+R(
	, const global fpxx* gi, global float* T // argument order is important
)+"#endif"+R( // TEMPERATURE
)+") {"+R( // update_fields()
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute update_fields() on halo
	const uchar flagsn = flags[n];
	const uchar flagsn_bo=flagsn&TYPE_BO, flagsn_su=flagsn&TYPE_SU; // extract boundary and surface flags
	if(flagsn_bo==TYPE_S||flagsn_su==TYPE_G) return; // don't update fields for boundary or gas lattice points

	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	float fhn[def_velocity_set]; // local DDFs
	load_f(n, fhn, fi, j, t); // perform streaming (part 2)

)+"#ifdef MOVING_BOUNDARIES"+R(
	if(flagsn_bo==TYPE_MS) apply_moving_boundaries(fhn, j, u, flags); // apply Dirichlet velocity boundaries if necessary (reads velocities of only neighboring boundary nodes, which do not change during simulation)
)+"#endif"+R( // MOVING_BOUNDARIES

	float rhon, uxn, uyn, uzn; // calculate local density and velocity for collision
	calculate_rho_u(fhn, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fi
	float fxn=fx, fyn=fy, fzn=fz; // force starts as constant volume force, can be modified before call of calculate_forcing_terms(...)

)+"#ifdef FORCE_FIELD"+R(
	{ // separate block to avoid variable name conflicts
		fxn += F[                 n]; // apply force field
		fyn += F[    def_N+(ulong)n];
		fzn += F[2ul*def_N+(ulong)n];
	}
)+"#endif"+R( // FORCE_FIELD

)+"#ifdef TEMPERATURE"+R(
	{ // separate block to avoid variable name conflicts
		uint j7[7]; // neighbors of D3Q7 subset
		neighbors_temperature(n, j7);
		float ghn[7]; // read from gA and stream to gh (D3Q7 subset, periodic boundary conditions)
		load_g(n, ghn, gi, j7, t); // perform streaming (part 2)
		float Tn;
		if(flagsn&TYPE_T) {
			Tn = T[n]; // apply preset temperature
		} else {
			Tn = 0.0f;
			for(uint i=0u; i<7u; i++) Tn += ghn[i]; // calculate temperature from g
			Tn += 1.0f; // add 1.0f last to avoid digit extinction effects when summing up gi (perturbation method / DDF-shifting)
			T[n] = Tn; // update temperature field
		}
		fxn -= fx*def_beta*(Tn-def_T_avg);
		fyn -= fy*def_beta*(Tn-def_T_avg);
		fzn -= fz*def_beta*(Tn-def_T_avg);
	}
)+"#endif"+R( // TEMPERATURE

	{ // separate block to avoid variable name conflicts
)+"#ifdef VOLUME_FORCE"+R( // apply force and collision operator, write to fi in video memory
		const float rho2 = 0.5f/rhon; // apply external volume force (Guo forcing, Krueger p.233f)
		uxn = clamp(fma(fxn, rho2, uxn), -def_c, def_c); // limit velocity (for stability purposes)
		uyn = clamp(fma(fyn, rho2, uyn), -def_c, def_c); // force term: F*dt/(2*rho)
		uzn = clamp(fma(fzn, rho2, uzn), -def_c, def_c);
)+"#else"+R( // VOLUME_FORCE
		uxn = clamp(uxn, -def_c, def_c); // limit velocity (for stability purposes)
		uyn = clamp(uyn, -def_c, def_c); // force term: F*dt/(2*rho)
		uzn = clamp(uzn, -def_c, def_c);
)+"#endif"+R( // VOLUME_FORCE
	}

)+"#ifndef EQUILIBRIUM_BOUNDARIES"+R(
	rho[               n] = rhon; // update density field
	u[                 n] = uxn; // update velocity field
	u[    def_N+(ulong)n] = uyn;
	u[2ul*def_N+(ulong)n] = uzn;
)+"#else"+R( // EQUILIBRIUM_BOUNDARIES
	if(flagsn_bo!=TYPE_E) { // only update fields for non-TYPE_E nodes
		rho[               n] = rhon; // update density field
		u[                 n] = uxn; // update velocity field
		u[    def_N+(ulong)n] = uyn;
		u[2ul*def_N+(ulong)n] = uzn;
	}
)+"#endif"+R( // EQUILIBRIUM_BOUNDARIES
} // update_fields()

)+"#ifdef FORCE_FIELD"+R(
)+R(kernel void calculate_force_on_boundaries(const global fpxx* fi, const global uchar* flags, const ulong t, global float* F) { // calculate force from the fluid on solid boundaries from fi directly
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute calculate_force_on_boundaries() on halo
	if((flags[n]&TYPE_BO)!=TYPE_S) return; // only continue for solid boundary nodes
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	float fhn[def_velocity_set]; // local DDFs
	load_f(n, fhn, fi, j, t); // perform streaming (part 2)
	float Fb=1.0f, fx=0.0f, fy=0.0f, fz=0.0f;
	calculate_rho_u(fhn, &Fb, &fx, &fy, &fz); // abuse calculate_rho_u() method for calculating force
	F[                 n] = 2.0f*fx*Fb; // 2 times because fi are reflected on solid boundary nodes (bounced-back)
	F[    def_N+(ulong)n] = 2.0f*fy*Fb;
	F[2ul*def_N+(ulong)n] = 2.0f*fz*Fb;
} // calculate_force_on_boundaries()
)+R(kernel void reset_force_field(global float* F) { // reset force field
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N) return; // execute reset_force_field() also on halo
	F[                 n] = 0.0f;
	F[    def_N+(ulong)n] = 0.0f;
	F[2ul*def_N+(ulong)n] = 0.0f;
} // reset_force_field()
)+R(void spread_force(volatile global float* F, const float3 p, const float3 Fn) {
	const float xa=p.x-0.5f+1.5f*def_Nx, ya=p.y-0.5f+1.5f*def_Ny, za=p.z-0.5f+1.5f*def_Nz; // subtract lattice offsets
	const uint xb=(uint)xa, yb=(uint)ya, zb=(uint)za; // integer casting to find bottom left corner
	const float x1=xa-(float)xb, y1=ya-(float)yb, z1=za-(float)zb; // calculate interpolation factors
	for(uint c=0u; c<8u; c++) { // count over eight corner points
		const uint i=(c&0x04u)>>2, j=(c&0x02u)>>1, k=c&0x01u; // disassemble c into corner indices ijk
		const uint x=(xb+i)%def_Nx, y=(yb+j)%def_Ny, z=(zb+k)%def_Nz; // calculate corner lattice positions
		const uint n = x+(y+z*def_Ny)*def_Nx; // calculate lattice linear index
		const float d = (1.0f-fabs(x1-(float)i))*(1.0f-fabs(y1-(float)j))*(1.0f-fabs(z1-(float)k)); // force spreading
		atomic_add_f(&F[                 n], Fn.x*d); // F[                 n] += Fn.x*d;
		atomic_add_f(&F[    def_N+(ulong)n], Fn.y*d); // F[    def_N+(ulong)n] += Fn.y*d;
		atomic_add_f(&F[2ul*def_N+(ulong)n], Fn.z*d); // F[2ul*def_N+(ulong)n] += Fn.z*d;
	}
} // spread_force()
)+"#endif"+R( // FORCE_FIELD

)+"#ifdef PARTICLES"+R(
)+R(kernel void integrate_particles)+"("+R(global float* particles, const global float* u // ) {
)+"#ifdef FORCE_FIELD"+R(
	, volatile global float* F, const float fx, const float fy, const float fz
)+"#endif"+R( // FORCE_FIELD
)+") {"+R( // integrate_particles()
	const uint n = get_global_id(0); // index of membrane points
	if(n>=(uint)def_particles_N) return;
	const float3 p0 = (float3)(particles[n], particles[def_particles_N+(ulong)n], particles[2ul*def_particles_N+(ulong)n]); // cache particle position
)+"#ifdef FORCE_FIELD"+R(
	if(def_particles_rho!=1.0f) {
		const float drho = def_particles_rho-1.0f; // density difference leads to particle buoyancy
		float3 Fn = (float3)(fx*drho, fy*drho, fz*drho); // F = F_p+F_f = (m_p-m_f)*g = (rho_p-rho_f)*g*V
		spread_force(F, p0, Fn); // do force spreading
	}
)+"#endif"+R( // FORCE_FIELD
	const float3 un = interpolate_u(mirror_position(p0), u); // trilinear interpolation of velocity at point p
	const float3 p = mirror_position(p0+un); // advect particles
	particles[                           n] = p.x;
	particles[    def_particles_N+(ulong)n] = p.y;
	particles[2ul*def_particles_N+(ulong)n] = p.z;
} // integrate_particles()
)+"#endif"+R( // PARTICLES



)+R(uint get_area(const uint direction) {
	const uint A[3] = { def_Ax, def_Ay, def_Az };
	return A[direction];
}
)+R(uint index_extract_p(const uint a, const uint direction) {
	const uint3 coordinates[3] = { (uint3)(def_Nx-2u, a%def_Ny, a/def_Ny), (uint3)(a/def_Nz, def_Ny-2u, a%def_Nz), (uint3)(a%def_Nx, a/def_Nx, def_Nz-2u) };
	return index(coordinates[direction]);
}
)+R(uint index_extract_m(const uint a, const uint direction) {
	const uint3 coordinates[3] = { (uint3)(       1u, a%def_Ny, a/def_Ny), (uint3)(a/def_Nz,        1u, a%def_Nz), (uint3)(a%def_Nx, a/def_Nx,        1u) };
	return index(coordinates[direction]);
}
)+R(uint index_insert_p(const uint a, const uint direction) {
	const uint3 coordinates[3] = { (uint3)(def_Nx-1u, a%def_Ny, a/def_Ny), (uint3)(a/def_Nz, def_Ny-1u, a%def_Nz), (uint3)(a%def_Nx, a/def_Nx, def_Nz-1u) };
	return index(coordinates[direction]);
}
)+R(uint index_insert_m(const uint a, const uint direction) {
	const uint3 coordinates[3] = { (uint3)(       0u, a%def_Ny, a/def_Ny), (uint3)(a/def_Nz,        0u, a%def_Nz), (uint3)(a%def_Nx, a/def_Nx,        0u) };
	return index(coordinates[direction]);
}

)+R(uint index_transfer(const uint side_i) {
	const uchar index_transfer_data[2u*def_dimensions*def_transfers] = {
)+"#if defined(D2Q9)"+R(
		1,  5,  7, // xp
		2,  6,  8, // xm
		3,  5,  8, // yp
		4,  6,  7  // ym
)+"#elif defined(D3Q15)"+R(
		1,  7,  9, 11, 14, // xp
		2,  8, 10, 12, 13, // xm
		3,  7,  9, 12, 13, // yp
		4,  8, 10, 11, 14, // ym
		5,  7, 10, 11, 13, // zp
		6,  8,  9, 12, 14  // zm
)+"#elif defined(D3Q19)"+R(
		1,  7,  9, 13, 15, // xp
		2,  8, 10, 14, 16, // xm
		3,  7, 11, 14, 17, // yp
		4,  8, 12, 13, 18, // ym
		5,  9, 11, 16, 18, // zp
		6, 10, 12, 15, 17  // zm
)+"#elif defined(D3Q27)"+R(
		1,  7,  9, 13, 15, 19, 21, 23, 26, // xp
		2,  8, 10, 14, 16, 20, 22, 24, 25, // xm
		3,  7, 11, 14, 17, 19, 21, 24, 25, // yp
		4,  8, 12, 13, 18, 20, 22, 23, 26, // ym
		5,  9, 11, 16, 18, 19, 22, 23, 25, // zp
		6, 10, 12, 15, 17, 20, 21, 24, 26  // zm
)+"#endif"+R( // D3Q27
	};
	return (uint)index_transfer_data[side_i];
}
)+R(void extract_fi(const uint a, const uint A, const uint n, const uint side, const ulong t, global fpxx_copy* transfer_buffer, const global fpxx_copy* fi) {
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	for(uint b=0u; b<def_transfers; b++) {
		const uint i = index_transfer(side*def_transfers+b);
		const ulong index = index_f(i%2u ? j[i] : n, t%2ul ? (i%2u ? i+1u : i-1u) : i); // Esoteric-Pull: standard store, or streaming part 1/2
		transfer_buffer[b*A+a] = fi[index]; // fpxx_copy allows direct copying without decompression+compression
	}
}
)+R(void insert_fi(const uint a, const uint A, const uint n, const uint side, const ulong t, const global fpxx_copy* transfer_buffer, global fpxx_copy* fi) {
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	for(uint b=0u; b<def_transfers; b++) {
		const uint i = index_transfer(side*def_transfers+b);
		const ulong index = index_f(i%2u ? n : j[i-1u], t%2ul ? i : (i%2u ? i+1u : i-1u)); // Esoteric-Pull: standard load, or streaming part 2/2
		fi[index] = transfer_buffer[b*A+a]; // fpxx_copy allows direct copying without decompression+compression
	}
}
)+R(kernel void transfer_extract_fi(const uint direction, const ulong t, global fpxx_copy* transfer_buffer_p, global fpxx_copy* transfer_buffer_m, const global fpxx_copy* fi) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	extract_fi(a, A, index_extract_p(a, direction), 2u*direction+0u, t, transfer_buffer_p, fi);
	extract_fi(a, A, index_extract_m(a, direction), 2u*direction+1u, t, transfer_buffer_m, fi);
}
)+R(kernel void transfer__insert_fi(const uint direction, const ulong t, const global fpxx_copy* transfer_buffer_p, const global fpxx_copy* transfer_buffer_m, global fpxx_copy* fi) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	insert_fi(a, A, index_insert_p(a, direction), 2u*direction+0u, t, transfer_buffer_p, fi);
	insert_fi(a, A, index_insert_m(a, direction), 2u*direction+1u, t, transfer_buffer_m, fi);
}

)+R(void extract_rho_u_flags(const uint a, const uint A, const uint n, global char* transfer_buffer, const global float* rho, const global float* u, const global uchar* flags) {
	((global float*)transfer_buffer)[      a] = rho[               n];
	((global float*)transfer_buffer)[    A+a] = u[                 n];
	((global float*)transfer_buffer)[ 2u*A+a] = u[    def_N+(ulong)n];
	((global float*)transfer_buffer)[ 3u*A+a] = u[2ul*def_N+(ulong)n];
	((global uchar*)transfer_buffer)[16u*A+a] = flags[             n];
}
)+R(void insert_rho_u_flags(const uint a, const uint A, const uint n, const global char* transfer_buffer, global float* rho, global float* u, global uchar* flags) {
	rho[               n] = ((const global float*)transfer_buffer)[      a];
	u[                 n] = ((const global float*)transfer_buffer)[    A+a];
	u[    def_N+(ulong)n] = ((const global float*)transfer_buffer)[ 2u*A+a];
	u[2ul*def_N+(ulong)n] = ((const global float*)transfer_buffer)[ 3u*A+a];
	flags[             n] = ((const global uchar*)transfer_buffer)[16u*A+a];
}
)+R(kernel void transfer_extract_rho_u_flags(const uint direction, const ulong t, global char* transfer_buffer_p, global char* transfer_buffer_m, const global float* rho, const global float* u, const global uchar* flags) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	extract_rho_u_flags(a, A, index_extract_p(a, direction), transfer_buffer_p, rho, u, flags);
	extract_rho_u_flags(a, A, index_extract_m(a, direction), transfer_buffer_m, rho, u, flags);
}
)+R(kernel void transfer__insert_rho_u_flags(const uint direction, const ulong t, const global char* transfer_buffer_p, const global char* transfer_buffer_m, global float* rho, global float* u, global uchar* flags) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	insert_rho_u_flags(a, A, index_insert_p(a, direction), transfer_buffer_p, rho, u, flags);
	insert_rho_u_flags(a, A, index_insert_m(a, direction), transfer_buffer_m, rho, u, flags);
}

)+"#ifdef SURFACE"+R(
)+R(kernel void transfer_extract_flags(const uint direction, const ulong t, global uchar* transfer_buffer_p, global uchar* transfer_buffer_m, const global uchar* flags) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	transfer_buffer_p[a] = flags[index_extract_p(a, direction)];
	transfer_buffer_m[a] = flags[index_extract_m(a, direction)];
}
)+R(kernel void transfer__insert_flags(const uint direction, const ulong t, const global uchar* transfer_buffer_p, const global uchar* transfer_buffer_m, global uchar* flags) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	flags[index_insert_p(a, direction)] = transfer_buffer_p[a];
	flags[index_insert_m(a, direction)] = transfer_buffer_m[a];
}

)+R(void extract_phi_massex_flags(const uint a, const uint A, const uint n, global char* transfer_buffer, const global float* phi, const global float* massex, const global uchar* flags) {
	((global float*)transfer_buffer)[     a] = phi   [n];
	((global float*)transfer_buffer)[   A+a] = massex[n];
	((global uchar*)transfer_buffer)[8u*A+a] = flags [n];
}
)+R(void insert_phi_massex_flags(const uint a, const uint A, const uint n, const global char* transfer_buffer, global float* phi, global float* massex, global uchar* flags) {
	phi   [n] = ((global float*)transfer_buffer)[     a];
	massex[n] = ((global float*)transfer_buffer)[   A+a];
	flags [n] = ((global uchar*)transfer_buffer)[8u*A+a];
}
)+R(kernel void transfer_extract_phi_massex_flags(const uint direction, const ulong t, global char* transfer_buffer_p, global char* transfer_buffer_m, const global float* phi, const global float* massex, const global uchar* flags) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	extract_phi_massex_flags(a, A, index_extract_p(a, direction), transfer_buffer_p, phi, massex, flags);
	extract_phi_massex_flags(a, A, index_extract_m(a, direction), transfer_buffer_m, phi, massex, flags);
}
)+R(kernel void transfer__insert_phi_massex_flags(const uint direction, const ulong t, const global char* transfer_buffer_p, const global char* transfer_buffer_m, global float* phi, global float* massex, global uchar* flags) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	insert_phi_massex_flags(a, A, index_insert_p(a, direction), transfer_buffer_p, phi, massex, flags);
	insert_phi_massex_flags(a, A, index_insert_m(a, direction), transfer_buffer_m, phi, massex, flags);
}
)+"#endif"+R( // SURFACE

)+"#ifdef TEMPERATURE"+R(
)+R(void extract_gi(const uint a, const uint n, const uint side, const ulong t, global fpxx_copy* transfer_buffer, const global fpxx_copy* gi) {
	uint j7[7u]; // neighbor indices
	neighbors_temperature(n, j7); // calculate neighbor indices
	const uint i = side+1u;
	const ulong index = index_f(i%2u ? j7[i] : n, t%2ul ? (i%2u ? i+1u : i-1u) : i); // Esoteric-Pull: standard store, or streaming part 1/2
	transfer_buffer[a] = gi[index]; // fpxx_copy allows direct copying without decompression+compression
}
)+R(void insert_gi(const uint a, const uint n, const uint side, const ulong t, const global fpxx_copy* transfer_buffer, global fpxx_copy* gi) {
	uint j7[7u]; // neighbor indices
	neighbors_temperature(n, j7); // calculate neighbor indices
	const uint i = side+1u;
	const ulong index = index_f(i%2u ? n : j7[i-1u], t%2ul ? i : (i%2u ? i+1u : i-1u)); // Esoteric-Pull: standard load, or streaming part 2/2
	gi[index] = transfer_buffer[a]; // fpxx_copy allows direct copying without decompression+compression
}
)+R(kernel void transfer_extract_gi(const uint direction, const ulong t, global fpxx_copy* transfer_buffer_p, global fpxx_copy* transfer_buffer_m, const global fpxx_copy* gi) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	extract_gi(a, index_extract_p(a, direction), 2u*direction+0u, t, transfer_buffer_p, gi);
	extract_gi(a, index_extract_m(a, direction), 2u*direction+1u, t, transfer_buffer_m, gi);
}
)+R(kernel void transfer__insert_gi(const uint direction, const ulong t, const global fpxx_copy* transfer_buffer_p, const global fpxx_copy* transfer_buffer_m, global fpxx_copy* gi) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	insert_gi(a, index_insert_p(a, direction), 2u*direction+0u, t, transfer_buffer_p, gi);
	insert_gi(a, index_insert_m(a, direction), 2u*direction+1u, t, transfer_buffer_m, gi);
}
)+"#endif"+R( // TEMPERATURE



)+R(kernel void voxelize_mesh(const uint direction, global fpxx* fi, global float* u, global uchar* flags, const ulong t, const uchar flag, const global float* p0, const global float* p1, const global float* p2, const global float* bbu) { // voxelize triangle mesh
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	const uint triangle_number = as_uint(bbu[0]);
	const float x0=bbu[ 1], y0=bbu[ 2], z0=bbu[ 3], x1=bbu[ 4], y1=bbu[ 5], z1=bbu[ 6];
	const float cx=bbu[ 7], cy=bbu[ 8], cz=bbu[ 9], ux=bbu[10], uy=bbu[11], uz=bbu[12], rx=bbu[13], ry=bbu[14], rz=bbu[15];

	const uint3 xyz = direction==0u ? (uint3)((uint)max(0, (int)x0-def_Ox), a%def_Ny, a/def_Ny) : direction==1u ? (uint3)(a/def_Nz, (uint)max(0, (int)y0-def_Oy), a%def_Nz) : (uint3)(a%def_Nx, a/def_Nx, (uint)max(0, (int)z0-def_Oz));
	const float3 p = position(xyz);
	const float3 offset = (float3)(0.5f*(float)((def_Nx-2u*(def_Dx>1u))*def_Dx)-0.5f, 0.5f*(float)((def_Ny-2u*(def_Dy>1u))*def_Dy)-0.5f, 0.5f*(float)((def_Nz-2u*(def_Dz>1u))*def_Dz)-0.5f)+(float3)(def_domain_offset_x, def_domain_offset_y, def_domain_offset_z);
	const float3 r_origin = p+offset;
	const float3 r_direction = (float3)((float)(direction==0u), (float)(direction==1u), (float)(direction==2u));
	uint intersections=0u, intersections_check=0u;
	ushort distances[64]; // allow up to 64 mesh intersections
	const bool condition = direction==0u ? r_origin.y<y0||r_origin.z<z0||r_origin.y>=y1||r_origin.z>=z1 : direction==1u ? r_origin.x<x0||r_origin.z<z0||r_origin.x>=x1||r_origin.z>=z1 : r_origin.x<x0||r_origin.y<y0||r_origin.x>=x1||r_origin.y>=y1;

	volatile local uint workgroup_condition; // use local memory optimization (~25% faster)
	workgroup_condition = 1u;
	barrier(CLK_LOCAL_MEM_FENCE);
	atomic_and(&workgroup_condition, (uint)condition);
	barrier(CLK_LOCAL_MEM_FENCE);
	const bool workgroup_all = (bool)workgroup_condition;
	if(workgroup_all) return; // return only if the entire workgroup is outside of the bounding-box of the mesh
	const uint lid = get_local_id(0);
	local float3 cache_p0[def_workgroup_size];
	local float3 cache_p1[def_workgroup_size];
	local float3 cache_p2[def_workgroup_size];
	for(uint i=0u; i<triangle_number; i+=def_workgroup_size) {
		const uint tx=3u*(i+lid), ty=tx+1u, tz=ty+1u;
		cache_p0[lid] = (float3)(p0[tx], p0[ty], p0[tz]);
		cache_p1[lid] = (float3)(p1[tx], p1[ty], p1[tz]);
		cache_p2[lid] = (float3)(p2[tx], p2[ty], p2[tz]);
		barrier(CLK_LOCAL_MEM_FENCE);
		for(int j=0; j<def_workgroup_size&&i+j<triangle_number; j++) {
			const float3 p0i=cache_p0[j], p1i=cache_p1[j], p2i=cache_p2[j];
			const float3 u=p1i-p0i, v=p2i-p0i, w=r_origin-p0i, h=cross(r_direction, v), q=cross(w, u); // bidirectional ray-triangle intersection (Moeller-Trumbore algorithm)
			const float f=1.0f/dot(u, h), s=f*dot(w, h), t=f*dot(r_direction, q), d=f*dot(v, q);
			if(s>=0.0f&&s<1.0f&&t>=0.0f&&s+t<1.0f) { // ray-triangle intersection ahead or behind
				if(d>0.0f) { // ray-triangle intersection ahead
					if(intersections<64u&&d<65536.0f) distances[intersections] = (ushort)d; // store distance to intersection in array as ushort
					intersections++;
				} else { // ray-triangle intersection behind
					intersections_check++; // cast a second ray to check if starting point is really inside (error correction)
				}
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if(condition) return; // extra workgroup threads outside of the bounding-box are not needed anymore, so return /**/

	/*if(condition) return; // don't use local memory (this also runs on old OpenCL 1.0 GPUs)
	for(uint i=0u; i<triangle_number; i++) {
		const float3 p0i = (float3)(p0[3u*i], p0[3u*i+1u], p0[3u*i+2u]);
		const float3 p1i = (float3)(p1[3u*i], p1[3u*i+1u], p1[3u*i+2u]);
		const float3 p2i = (float3)(p2[3u*i], p2[3u*i+1u], p2[3u*i+2u]);
		const float3 u=p1i-p0i, v=p2i-p0i, w=r_origin-p0i, h=cross(r_direction, v), q=cross(w, u); // bidirectional ray-triangle intersection (Moeller-Trumbore algorithm)
		const float f=1.0f/dot(u, h), s=f*dot(w, h), t=f*dot(r_direction, q), d=f*dot(v, q);
		if(s>=0.0f&&s<1.0f&&t>=0.0f&&s+t<1.0f) { // ray-triangle intersection ahead or behind
			if(d>0.0f) { // ray-triangle intersection ahead
				if(intersections<64u&&d<65536.0f) distances[intersections] = (ushort)d; // store distance to intersection in array as ushort
				intersections++;
			} else { // ray-triangle intersection behind
				intersections_check++; // cast a second ray to check if starting point is really inside (error correction)
			}
		}
	}/**/

	for(int i=1; i<(int)intersections; i++) { // insertion-sort distances
		ushort t = distances[i];
		int j = i-1;
		while(distances[j]>t&&j>=0) {
			distances[j+1] = distances[j];
			j--;
		}
		distances[j+1] = t;
	}

	bool inside = (intersections%2u)&&(intersections_check%2u);
	const bool set_u = sq(ux)+sq(uy)+sq(uz)+sq(rx)+sq(ry)+sq(rz)>0.0f;
	uint intersection = intersections%2u!=intersections_check%2u; // iterate through column, start with 0 regularly, start with 1 if forward and backward intersection count evenness differs (error correction)
	const uint h0 = direction==0u ? xyz.x : direction==1u ? xyz.y : xyz.z;
	const uint hmax = direction==0u ? (uint)clamp((int)x1-def_Ox, 0, (int)def_Nx-1) : direction==1u ? (uint)clamp((int)y1-def_Oy, 0, (int)def_Ny-1) : (uint)clamp((int)z1-def_Oz, 0, (int)def_Nz-1);
	const uint hmesh = h0+(uint)distances[intersections-1u];
	for(uint h=h0; h<hmax; h++) {
		while(intersection<intersections&&h>h0+(uint)distances[intersection]) {
			inside = !inside; // passed mesh intersection, so switch inside/outside state
			intersection++;
		}
		inside &= (intersection<intersections&&h<hmesh); // point must be outside if there are no more ray-mesh intersections ahead (error correction)
		const ulong n = index((uint3)(direction==0u?h:xyz.x, direction==1u?h:xyz.y, direction==2u?h:xyz.z));
		uchar flagsn = flags[n];
		if(inside) {
			flagsn = (flagsn&~TYPE_BO)|flag;
			if(set_u) {
				const float3 p = position(coordinates(n))+offset;
				const float3 un = (float3)(ux, uy, uz)+cross((float3)(cx, cy, cz)-p, (float3)(rx, ry, rz));
				u[                 n] = un.x;
				u[    def_N+(ulong)n] = un.y;
				u[2ul*def_N+(ulong)n] = un.z;
			}
		} else {
			if(set_u) {
				const float3 un = (float3)(u[n], u[def_N+(ulong)n], u[2ul*def_N+(ulong)n]); // for velocity voxelization, only clear moving boundaries
				if((flagsn&TYPE_BO)==TYPE_S) { // reconstruct DDFs when boundary point is converted to fluid
					uint j[def_velocity_set]; // neighbor indices
					neighbors(n, j); // calculate neighbor indices
					float feq[def_velocity_set]; // f_equilibrium
					calculate_f_eq(1.0f, un.x, un.y, un.z, feq);
					store_f(n, feq, fi, j, t); // write to fi
				}
				if(sq(un.x)+sq(un.y)+sq(un.z)>0.0f) {
					flagsn = (flagsn&TYPE_BO)==TYPE_MS ? flagsn&~TYPE_MS : flagsn&~flag;
				}
			} else {
				flagsn = (flagsn&TYPE_BO)==TYPE_MS ? flagsn&~TYPE_MS : flagsn&~flag;
			}
		}
		flags[n] = flagsn;
	}
} // voxelize_mesh()

)+R(kernel void unvoxelize_mesh(global uchar* flags, const uchar flag, float x0, float y0, float z0, float x1, float y1, float z1) { // remove voxelized triangle mesh
	const uint n = get_global_id(0);
	const float3 p = position(coordinates(n))+(float3)(0.5f*(float)((def_Nx-2u*(def_Dx>1u))*def_Dx)-0.5f, 0.5f*(float)((def_Ny-2u*(def_Dy>1u))*def_Dy)-0.5f, 0.5f*(float)((def_Nz-2u*(def_Dz>1u))*def_Dz)-0.5f)+(float3)(def_domain_offset_x, def_domain_offset_y, def_domain_offset_z);
	if(p.x>=x0-1.0f&&p.y>=y0-1.0f&&p.z>=z0-1.0f&&p.x<=x1+1.0f&&p.y<=y1+1.0f&&p.z<=z1+1.0f) flags[n] &= ~flag;
} // unvoxelize_mesh()



// ################################################## graphics code ##################################################

)+"#ifdef GRAPHICS"+R(

)+"#ifndef FORCE_FIELD"+R( // render flags as grid
)+R(kernel void graphics_flags(const global uchar* flags, const global float* camera, global int* bitmap, global int* zbuffer) {
)+"#else"+R( // FORCE_FIELD
)+R(kernel void graphics_flags(const global uchar* flags, const global float* camera, global int* bitmap, global int* zbuffer, const global float* F) {
)+"#endif"+R( // FORCE_FIELD
	const uint n = get_global_id(0);
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute graphics_flags() on halo
	const uchar flagsn = flags[n]; // cache flags
	const uchar flagsn_bo = flagsn&TYPE_BO; // extract boundary flags
	if(flagsn==0u||flagsn==TYPE_G) return; // don't draw regular fluid nodes
	//if(flagsn&TYPE_SU) return; // don't draw surface
	float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
	for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
	uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
	calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
	const uint3 xyz = coordinates(n);
	const float3 p = position(xyz);
	const int c =  // coloring scheme
		flagsn_bo==TYPE_S ? COLOR_S : // solid boundary
		((flagsn&TYPE_T)&&flagsn_bo==TYPE_E) ? color_mix(COLOR_T, COLOR_E, 0.5f) : // both temperature boundary and equilibrium boundary
		((flagsn&TYPE_T)&&flagsn_bo==TYPE_MS) ? color_mix(COLOR_T, COLOR_M, 0.5f) : // both temperature boundary and moving boundary
		flagsn&TYPE_T ? COLOR_T : // temperature boundary
		flagsn_bo==TYPE_E ? COLOR_E : // equilibrium boundary
		flagsn_bo==TYPE_MS ? COLOR_M : // moving boundary
		flagsn&TYPE_F ? COLOR_F : // fluid
		flagsn&TYPE_I ? COLOR_I : // interface
		flagsn&TYPE_X ? COLOR_X : // reserved type X
		flagsn&TYPE_Y ? COLOR_Y : // reserved type Y
		COLOR_0; // regular or gas node
	//draw_point(p, c, camera_cache, bitmap, zbuffer); // draw one pixel for every boundary node
	uint t;
	t = xp+y0+z0; const bool not_xp = xyz.x<def_Nx-1u && flagsn==flags[t] && !is_halo(t); // +00
	t = xm+y0+z0; const bool not_xm = xyz.x>       0u && flagsn==flags[t] && !is_halo(t); // -00
	t = x0+yp+z0; const bool not_yp = xyz.y<def_Ny-1u && flagsn==flags[t] && !is_halo(t); // 0+0
	t = x0+ym+z0; const bool not_ym = xyz.y>       0u && flagsn==flags[t] && !is_halo(t); // 0-0
	t = x0+y0+zp; const bool not_zp = xyz.z<def_Nz-1u && flagsn==flags[t] && !is_halo(t); // 00+
	t = x0+y0+zm; const bool not_zm = xyz.z>       0u && flagsn==flags[t] && !is_halo(t); // 00-
	const float3 p0 = (float3)(p.x-0.5f, p.y-0.5f, p.z-0.5f); // ---
	const float3 p1 = (float3)(p.x+0.5f, p.y+0.5f, p.z+0.5f); // +++
	const float3 p2 = (float3)(p.x-0.5f, p.y-0.5f, p.z+0.5f); // --+
	const float3 p3 = (float3)(p.x+0.5f, p.y+0.5f, p.z-0.5f); // ++-
	const float3 p4 = (float3)(p.x-0.5f, p.y+0.5f, p.z-0.5f); // -+-
	const float3 p5 = (float3)(p.x+0.5f, p.y-0.5f, p.z+0.5f); // +-+
	const float3 p6 = (float3)(p.x+0.5f, p.y-0.5f, p.z-0.5f); // +--
	const float3 p7 = (float3)(p.x-0.5f, p.y+0.5f, p.z+0.5f); // -++
	if(!(not_xm||not_ym)) draw_line(p0, p2, c, camera_cache, bitmap, zbuffer); // to draw the entire surface, replace || by &&
	if(!(not_xm||not_zm)) draw_line(p0, p4, c, camera_cache, bitmap, zbuffer);
	if(!(not_ym||not_zm)) draw_line(p0, p6, c, camera_cache, bitmap, zbuffer);
	if(!(not_xp||not_yp)) draw_line(p1, p3, c, camera_cache, bitmap, zbuffer);
	if(!(not_xp||not_zp)) draw_line(p1, p5, c, camera_cache, bitmap, zbuffer);
	if(!(not_yp||not_zp)) draw_line(p1, p7, c, camera_cache, bitmap, zbuffer);
	if(!(not_ym||not_zp)) draw_line(p2, p5, c, camera_cache, bitmap, zbuffer);
	if(!(not_xm||not_zp)) draw_line(p2, p7, c, camera_cache, bitmap, zbuffer);
	if(!(not_yp||not_zm)) draw_line(p3, p4, c, camera_cache, bitmap, zbuffer);
	if(!(not_xp||not_zm)) draw_line(p3, p6, c, camera_cache, bitmap, zbuffer);
	if(!(not_xm||not_yp)) draw_line(p4, p7, c, camera_cache, bitmap, zbuffer);
	if(!(not_xp||not_ym)) draw_line(p5, p6, c, camera_cache, bitmap, zbuffer);
)+"#ifdef FORCE_FIELD"+R(
	if(flagsn_bo==TYPE_S) {
		const float3 Fn = def_scale_F*(float3)(F[n], F[def_N+(ulong)n], F[2ul*def_N+(ulong)n]);
		const float Fnl = length(Fn);
		if(Fnl>0.0f) {
			const int c = iron_color(255.0f*Fnl); // color boundaries depending on the force on them
			draw_line(p, p+5.0f*Fn, c, camera_cache, bitmap, zbuffer); // draw colored force vectors
		}
	}
)+"#endif"+R( // FORCE_FIELD
}/**/

/*)+"#ifndef FORCE_FIELD"+R( // render solid boundaries with marching-cubes
)+R(kernel void graphics_flags(const global uchar* flags, const global float* camera, global int* bitmap, global int* zbuffer) {
)+"#else"+R( // FORCE_FIELD
)+R(kernel void graphics_flags(const global uchar* flags, const global float* camera, global int* bitmap, global int* zbuffer, const global float* F) {
)+"#endif"+R( // FORCE_FIELD
	const uint n = get_global_id(0);
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute graphics_flags() on halo
	const uint3 xyz = coordinates(n);
	if(xyz.x>=def_Nx-1u||xyz.y>=def_Ny-1u||xyz.z>=def_Nz-1u) return;
	//if(xyz.x==0u||xyz.y==0u||xyz.z==0u||xyz.x>=def_Nx-2u||xyz.y>=def_Ny-2u||xyz.z>=def_Nz-2u) return;
	uint j[8];
	const uint x0 =  xyz.x; // cube stencil
	const uint xp =  xyz.x+1u;
	const uint y0 =  xyz.y    *def_Nx;
	const uint yp = (xyz.y+1u)*def_Nx;
	const uint z0 =  xyz.z    *def_Ny*def_Nx;
	const uint zp = (xyz.z+1u)*def_Ny*def_Nx;
	j[0] = n       ; // 000
	j[1] = xp+y0+z0; // +00
	j[2] = xp+y0+zp; // +0+
	j[3] = x0+y0+zp; // 00+
	j[4] = x0+yp+z0; // 0+0
	j[5] = xp+yp+z0; // ++0
	j[6] = xp+yp+zp; // +++
	j[7] = x0+yp+zp; // 0++
	float v[8];
	for(uint i=0u; i<8u; i++) v[i] = (float)((flags[j[i]]&TYPE_BO)==TYPE_S);
	float3 triangles[15]; // maximum of 5 triangles with 3 vertices each
	const uint tn = marching_cubes(v, 0.5f, triangles); // run marching cubes algorithm
	if(tn==0u) return;
	float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
	for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
	const float3 offset = (float3)((float)xyz.x+0.5f-0.5f*(float)def_Nx, (float)xyz.y+0.5f-0.5f*(float)def_Ny, (float)xyz.z+0.5f-0.5f*(float)def_Nz);
	for(uint i=0u; i<tn; i++) {
		const float3 p0 = triangles[3u*i   ]+offset;
		const float3 p1 = triangles[3u*i+1u]+offset;
		const float3 p2 = triangles[3u*i+2u]+offset;
		const float3 p=(p0+p1+p2)/3.0f, normal=cross(p1-p0, p2-p0);
		const int c = lighting(191<<16|191<<8|191, p, normal, camera_cache);
		draw_triangle(p0, p1, p2, c, camera_cache, bitmap, zbuffer);
	}
)+"#ifdef FORCE_FIELD"+R(
	const uchar flagsn_bo = flags[n]&TYPE_BO;
	const float3 p = position(xyz);
	if(flagsn_bo==TYPE_S) {
		const float3 Fn = def_scale_F*(float3)(F[n], F[def_N+(ulong)n], F[2ul*def_N+(ulong)n]);
		const float Fnl = length(Fn);
		if(Fnl>0.0f) {
			const int c = iron_color(255.0f*Fnl); // color boundaries depending on the force on them
			draw_line(p, p+5.0f*Fn, c, camera_cache, bitmap, zbuffer); // draw colored force vectors
		}
	}
)+"#endif"+R( // FORCE_FIELD
}/**/

)+R(kernel void graphics_field(const global uchar* flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer) {
	const uint n = get_global_id(0);
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute graphics_field() on halo
)+"#ifndef MOVING_BOUNDARIES"+R(
	if(flags[n]&(TYPE_S|TYPE_E|TYPE_I|TYPE_G)) return;
)+"#else"+R( // EQUILIBRIUM_BOUNDARIES
	if(flags[n]&(TYPE_I|TYPE_G)) return;
)+"#endif"+R( // EQUILIBRIUM_BOUNDARIES
	float3 un = load_u(n, u); // cache velocity
	const float ul = length(un);
	if(def_scale_u*ul<0.1f) return; // don't draw lattice points where the velocity is lower than this threshold
	float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
	for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
	const float3 p = position(coordinates(n));
	const int c = iron_color(255.0f*def_scale_u*ul); // coloring by velocity
	draw_line(p-(0.5f/ul)*un, p+(0.5f/ul)*un, c, camera_cache, bitmap, zbuffer);
}

)+"#ifndef GRAPHICS_TEMPERATURE"+R(
)+R(kernel void graphics_streamline(const global uchar* flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer) {
)+"#else"+R( // GRAPHICS_TEMPERATURE
)+R(kernel void graphics_streamline(const global uchar* flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer, const global float* T) {
)+"#endif"+R( // GRAPHICS_TEMPERATURE
	const uint n = get_global_id(0);
)+"#ifndef D2Q9"+R(
	if(n>=(def_Nx/def_streamline_sparse)*(def_Ny/def_streamline_sparse)*(def_Nz/def_streamline_sparse)) return;
	const uint z = n/((def_Nx/def_streamline_sparse)*(def_Ny/def_streamline_sparse)); // disassemble 1D index to 3D coordinates
	const uint t = n%((def_Nx/def_streamline_sparse)*(def_Ny/def_streamline_sparse));
	const uint y = t/(def_Nx/def_streamline_sparse);
	const uint x = t%(def_Nx/def_streamline_sparse);
	float3 p = (float)def_streamline_sparse*((float3)((float)x+0.5f, (float)y+0.5f, (float)z+0.5f))-0.5f*((float3)((float)def_Nx, (float)def_Ny, (float)def_Nz));
)+"#else"+R( // D2Q9
	if(n>=(def_Nx/def_streamline_sparse)*(def_Ny/def_streamline_sparse)) return;
	const uint y = n/(def_Nx/def_streamline_sparse); // disassemble 1D index to 3D coordinates
	const uint x = n%(def_Nx/def_streamline_sparse);
	float3 p = ((float3)((float)def_streamline_sparse*((float)x+0.5f), (float)def_streamline_sparse*((float)y+0.5f), 0.5f))-0.5f*((float3)((float)def_Nx, (float)def_Ny, (float)def_Nz));
)+"#endif"+R( // D2Q9
	float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
	for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
	const float hLx=0.5f*(float)(def_Nx-2u*(def_Dx>1u)), hLy=0.5f*(float)(def_Ny-2u*(def_Dy>1u)), hLz=0.5f*(float)(def_Nz-2u*(def_Dz>1u));
	//draw_circle(p, 0.5f*def_streamline_sparse, 0xFFFFFF, camera_cache, bitmap, zbuffer);
	for(float dt=-1.0f; dt<=1.0f; dt+=2.0f) { // integrate forward and backward in time
		float3 p0, p1=p;
		for(uint l=0u; l<def_streamline_length/2u; l++) {
			const uint x = (uint)(p1.x+1.5f*(float)def_Nx)%def_Nx;
			const uint y = (uint)(p1.y+1.5f*(float)def_Ny)%def_Ny;
			const uint z = (uint)(p1.z+1.5f*(float)def_Nz)%def_Nz;
			const uint n = x+(y+z*def_Ny)*def_Nx;
			if(flags[n]&(TYPE_S|TYPE_E|TYPE_I|TYPE_G)) return;
			const float3 un = load_u(n, u); // interpolate_u(p1, u)
			const float ul = length(un);
			p0 = p1;
			p1 += (dt/ul)*un; // integrate forward in time
			if(def_scale_u*ul<0.1f||p1.x<-hLx||p1.x>hLx||p1.y<-hLy||p1.y>hLy||p1.z<-hLz||p1.z>hLz) break;
)+"#ifndef GRAPHICS_TEMPERATURE"+R(
			const int c = iron_color(255.0f*def_scale_u*ul);
)+"#else"+R( // GRAPHICS_TEMPERATURE
			const int c = iron_color(167.0f+255.0f*(T[n]-def_T_avg));
)+"#endif"+R( // GRAPHICS_TEMPERATURE
			draw_line(p0, p1, c, camera_cache, bitmap, zbuffer);
		}
	}
}

)+R(kernel void graphics_q_field(const global uchar* flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer) {
	const uint n = get_global_id(0);
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute graphics_q_field() on halo
	if(flags[n]&(TYPE_S|TYPE_E|TYPE_I|TYPE_G)) return;
	float3 un = load_u(n, u); // cache velocity
	const float ul = length(un);
	const float Q = calculate_Q(n, u);
	if(Q<def_scale_Q_min||ul==0.0f) return; // don't draw lattice points where the velocity is very low
	float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
	for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
	const float3 p = position(coordinates(n));
	const int c = rainbow_color(255.0f*def_scale_u*ul); // coloring by velocity
	draw_line(p-(0.5f/ul)*un, p+(0.5f/ul)*un, c, camera_cache, bitmap, zbuffer);
}

)+R(kernel void graphics_q(const global uchar* flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer) {
	const uint n = get_global_id(0);
	if(is_halo_q(n)) return; // don't execute graphics_q_field() on marching-cubes halo
	const uint3 xyz = coordinates(n);
	if(xyz.x>=def_Nx-1u||xyz.y>=def_Ny-1u||xyz.z>=def_Nz-1u) return;
	const uint x0 =  xyz.x; // cube stencil
	const uint xp =  xyz.x+1u;
	const uint y0 =  xyz.y    *def_Nx;
	const uint yp = (xyz.y+1u)*def_Nx;
	const uint z0 =  xyz.z    *def_Ny*def_Nx;
	const uint zp = (xyz.z+1u)*def_Ny*def_Nx;
	const uint xq =  (xyz.x       +2u)%def_Nx; // central difference stencil on each cube corner point
	const uint xm =  (xyz.x+def_Nx-1u)%def_Nx;
	const uint yq = ((xyz.y       +2u)%def_Ny)*def_Nx;
	const uint ym = ((xyz.y+def_Ny-1u)%def_Ny)*def_Nx;
	const uint zq = ((xyz.z       +2u)%def_Nz)*def_Ny*def_Nx;
	const uint zm = ((xyz.z+def_Nz-1u)%def_Nz)*def_Ny*def_Nx;
	uint j[32];
	j[ 0] = n       ; // 000 // cube stencil
	j[ 1] = xp+y0+z0; // +00
	j[ 2] = xp+y0+zp; // +0+
	j[ 3] = x0+y0+zp; // 00+
	j[ 4] = x0+yp+z0; // 0+0
	j[ 5] = xp+yp+z0; // ++0
	j[ 6] = xp+yp+zp; // +++
	j[ 7] = x0+yp+zp; // 0++
	j[ 8] = xm+y0+z0; // -00 // central difference stencil on each cube corner point
	j[ 9] = x0+ym+z0; // 0-0
	j[10] = x0+y0+zm; // 00-
	j[11] = xq+y0+z0; // #00
	j[12] = xp+ym+z0; // +-0
	j[13] = xp+y0+zm; // +0-
	j[14] = xq+y0+zp; // #0+
	j[15] = xp+ym+zp; // +-+
	j[16] = xp+y0+zq; // +0#
	j[17] = xm+y0+zp; // -0+
	j[18] = x0+ym+zp; // 0-+
	j[19] = x0+y0+zq; // 00#
	j[20] = xm+yp+z0; // -+0
	j[21] = x0+yq+z0; // 0#0
	j[22] = x0+yp+zm; // 0+-
	j[23] = xq+yp+z0; // #+0
	j[24] = xp+yq+z0; // +#0
	j[25] = xp+yp+zm; // ++-
	j[26] = xq+yp+zp; // #++
	j[27] = xp+yq+zp; // +#+
	j[28] = xp+yp+zq; // ++#
	j[29] = xm+yp+zp; // -++
	j[30] = x0+yq+zp; // 0#+
	j[31] = x0+yp+zq; // 0+#
	uchar flags_cell = 0u;
	for(uint i=0u; i<32u; i++) flags_cell |= flags[j[i]];
	if(flags_cell&(TYPE_S|TYPE_E|TYPE_I|TYPE_G)) return;
	float3 uj[8], u0[6], u1[6], u2[6], u3[6], u4[6], u5[6], u6[6], u7[6]; // don't load any velocity twice from global memory
	for(uint i=0u; i<8u; i++) uj[i] = load_u(j[i], u);
	u0[0]=      uj[ 1]    ; u0[1]=load_u(j[ 8], u); u0[2]=      uj[ 4]    ; u0[3]=load_u(j[ 9], u); u0[4]=      uj[ 3]    ; u0[5]=load_u(j[10], u);
	u1[0]=load_u(j[11], u); u1[1]=      uj[ 0]    ; u1[2]=      uj[ 5]    ; u1[3]=load_u(j[12], u); u1[4]=      uj[ 2]    ; u1[5]=load_u(j[13], u);
	u2[0]=load_u(j[14], u); u2[1]=      uj[ 3]    ; u2[2]=      uj[ 6]    ; u2[3]=load_u(j[15], u); u2[4]=load_u(j[16], u); u2[5]=      uj[ 1]    ;
	u3[0]=      uj[ 2]    ; u3[1]=load_u(j[17], u); u3[2]=      uj[ 7]    ; u3[3]=load_u(j[18], u); u3[4]=load_u(j[19], u); u3[5]=      uj[ 0]    ;
	u4[0]=      uj[ 5]    ; u4[1]=load_u(j[20], u); u4[2]=load_u(j[21], u); u4[3]=      uj[ 0]    ; u4[4]=      uj[ 7]    ; u4[5]=load_u(j[22], u);
	u5[0]=load_u(j[23], u); u5[1]=      uj[ 4]    ; u5[2]=load_u(j[24], u); u5[3]=      uj[ 1]    ; u5[4]=      uj[ 6]    ; u5[5]=load_u(j[25], u);
	u6[0]=load_u(j[26], u); u6[1]=      uj[ 7]    ; u6[2]=load_u(j[27], u); u6[3]=      uj[ 2]    ; u6[4]=load_u(j[28], u); u6[5]=      uj[ 5]    ;
	u7[0]=      uj[ 6]    ; u7[1]=load_u(j[29], u); u7[2]=load_u(j[30], u); u7[3]=      uj[ 3]    ; u7[4]=load_u(j[31], u); u7[5]=      uj[ 4]    ;
	float v[8];
	v[0] = calculate_Q_cached(u0);
	v[1] = calculate_Q_cached(u1);
	v[2] = calculate_Q_cached(u2);
	v[3] = calculate_Q_cached(u3);
	v[4] = calculate_Q_cached(u4);
	v[5] = calculate_Q_cached(u5);
	v[6] = calculate_Q_cached(u6);
	v[7] = calculate_Q_cached(u7);
	float3 triangles[15]; // maximum of 5 triangles with 3 vertices each
	const uint tn = marching_cubes(v, def_scale_Q_min, triangles); // run marching cubes algorithm
	if(tn==0u) return;
	float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
	for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
	const float3 offset = (float3)((float)xyz.x+0.5f-0.5f*(float)def_Nx, (float)xyz.y+0.5f-0.5f*(float)def_Ny, (float)xyz.z+0.5f-0.5f*(float)def_Nz);
	for(uint i=0u; i<tn; i++) {
		const float3 p0 = triangles[3u*i   ]; // triangle coordinates in [0,1] (local cell)
		const float3 p1 = triangles[3u*i+1u];
		const float3 p2 = triangles[3u*i+2u];
		const float3 normal = cross(p1-p0, p2-p0);
		int c0, c1, c2;
		{
			const float x1=p0.x, y1=p0.y, z1=p0.z, x0=1.0f-x1, y0=1.0f-y1, z0=1.0f-z1; // calculate interpolation factors
			const float3 ui = (x0*y0*z0)*uj[0]+(x1*y0*z0)*uj[1]+(x1*y0*z1)*uj[2]+(x0*y0*z1)*uj[3]+(x0*y1*z0)*uj[4]+(x1*y1*z0)*uj[5]+(x1*y1*z1)*uj[6]+(x0*y1*z1)*uj[7]; // perform trilinear interpolation
			c0 = lighting(rainbow_color(255.0f*def_scale_u*length(ui)), p0+offset, normal, camera_cache); // rainbow_color(255.0f*def_scale_u*length(ui));
		} {
			const float x1=p1.x, y1=p1.y, z1=p1.z, x0=1.0f-x1, y0=1.0f-y1, z0=1.0f-z1; // calculate interpolation factors
			const float3 ui = (x0*y0*z0)*uj[0]+(x1*y0*z0)*uj[1]+(x1*y0*z1)*uj[2]+(x0*y0*z1)*uj[3]+(x0*y1*z0)*uj[4]+(x1*y1*z0)*uj[5]+(x1*y1*z1)*uj[6]+(x0*y1*z1)*uj[7]; // perform trilinear interpolation
			c1 = lighting(rainbow_color(255.0f*def_scale_u*length(ui)), p1+offset, normal, camera_cache); // rainbow_color(255.0f*def_scale_u*length(ui));
		} {
			const float x1=p2.x, y1=p2.y, z1=p2.z, x0=1.0f-x1, y0=1.0f-y1, z0=1.0f-z1; // calculate interpolation factors
			const float3 ui = (x0*y0*z0)*uj[0]+(x1*y0*z0)*uj[1]+(x1*y0*z1)*uj[2]+(x0*y0*z1)*uj[3]+(x0*y1*z0)*uj[4]+(x1*y1*z0)*uj[5]+(x1*y1*z1)*uj[6]+(x0*y1*z1)*uj[7]; // perform trilinear interpolation
			c2 = lighting(rainbow_color(255.0f*def_scale_u*length(ui)), p2+offset, normal, camera_cache); // rainbow_color(255.0f*def_scale_u*length(ui));
		}
		draw_triangle_interpolated(p0+offset, p1+offset, p2+offset, c0, c1, c2, camera_cache, bitmap, zbuffer); // draw triangle with interpolated colors
	}
}

)+"#ifdef SURFACE"+R(
)+R(kernel void graphics_rasterize_phi(const global float* phi, const global float* camera, global int* bitmap, global int* zbuffer) { // marching cubes
	const uint n = get_global_id(0);
	const uint3 xyz = coordinates(n);
	if(xyz.x>=def_Nx-1u||xyz.y>=def_Ny-1u||xyz.z>=def_Nz-1u) return;
	uint j[8];
	const uint x0 =  xyz.x; // cube stencil
	const uint xp =  xyz.x+1u;
	const uint y0 =  xyz.y    *def_Nx;
	const uint yp = (xyz.y+1u)*def_Nx;
	const uint z0 =  xyz.z    *def_Ny*def_Nx;
	const uint zp = (xyz.z+1u)*def_Ny*def_Nx;
	j[0] = n       ; // 000
	j[1] = xp+y0+z0; // +00
	j[2] = xp+y0+zp; // +0+
	j[3] = x0+y0+zp; // 00+
	j[4] = x0+yp+z0; // 0+0
	j[5] = xp+yp+z0; // ++0
	j[6] = xp+yp+zp; // +++
	j[7] = x0+yp+zp; // 0++
	float v[8];
	for(uint i=0u; i<8u; i++) v[i] = phi[j[i]];
	float3 triangles[15]; // maximum of 5 triangles with 3 vertices each
	const uint tn = marching_cubes(v, 0.5f, triangles); // run marching cubes algorithm
	if(tn==0u) return;
	float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
	for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
	const float3 offset = (float3)((float)xyz.x+0.5f-0.5f*(float)def_Nx, (float)xyz.y+0.5f-0.5f*(float)def_Ny, (float)xyz.z+0.5f-0.5f*(float)def_Nz);
	for(uint i=0u; i<tn; i++) {
		const float3 p0 = triangles[3u*i   ]+offset;
		const float3 p1 = triangles[3u*i+1u]+offset;
		const float3 p2 = triangles[3u*i+2u]+offset;
		const float3 p=(p0+p1+p2)/3.0f, normal=cross(p1-p0, p2-p0);
		const int c = lighting(55<<16|155<<8|255, p, normal, camera_cache);
		draw_triangle(p0, p1, p2, c, camera_cache, bitmap, zbuffer);
		//draw_line(p0, p1, c, camera_cache, bitmap, zbuffer); // wireframe rendering
		//draw_line(p0, p2, c, camera_cache, bitmap, zbuffer);
		//draw_line(p1, p2, c, camera_cache, bitmap, zbuffer);
	}
}

)+R(int raytrace_phi_next_ray(const ray reflection, const ray transmission, const float reflectivity, const float transmissivity, const global float* phi, const global uchar* flags, const global int* skybox) {
	int color_reflect=0, color_transmit=0;
	ray reflection_next, transmission_next;
	float reflection_reflectivity, reflection_transmissivity, transmission_reflectivity, transmission_transmissivity;
	if(raytrace_phi(reflection, &reflection_next, &transmission_next, &reflection_reflectivity, &reflection_transmissivity, phi, flags, skybox, def_Nx, def_Ny, def_Nz)) {
		color_reflect = last_ray_reflectivity(reflection_next, transmission_next, reflection_reflectivity, reflection_transmissivity, skybox);
	} else {
		color_reflect = skybox_color(reflection, skybox);
	}
	if(raytrace_phi(transmission, &reflection_next, &transmission_next, &transmission_reflectivity, &transmission_transmissivity, phi, flags, skybox, def_Nx, def_Ny, def_Nz)) {
		color_transmit = last_ray_reflectivity(reflection_next, transmission_next, transmission_reflectivity, transmission_transmissivity, skybox);
	} else {
		color_transmit = skybox_color(transmission, skybox);
	}
	return color_mix(color_reflect, color_mix(color_transmit, def_absorption_color, transmissivity), reflectivity);
}
)+R(int raytrace_phi_next_ray_mirror(const ray reflection, const global float* phi, const global uchar* flags, const global int* skybox) {
	int color_reflect = 0;
	ray reflection_next;
	if(raytrace_phi_mirror(reflection, &reflection_next, phi, flags, skybox, def_Nx, def_Ny, def_Nz)) {
		color_reflect = skybox_color(reflection_next, skybox);
	} else {
		color_reflect = skybox_color(reflection, skybox);
	}
	return color_reflect;
}

)+R(kernel void graphics_raytrace_phi(const global float* phi, const global uchar* flags, const global int* skybox, const global float* camera, global int* bitmap) { // marching cubes
	const uint gid = get_global_id(0); // workgroup size alignment is critical
	const uint lid = get_local_id(0); // make workgropus not horizontal stripes of pixels, but 8x8 rectangular (close to square) tiles
	const uint lsi = get_local_size(0); // (50% performance boost due to more coalesced memory access)
	const uint tile_width=8u, tile_height=lsi/tile_width, tiles_x=def_screen_width/tile_width;
	const int lx=lid%tile_width, ly=lid/tile_width;
	const int tx=(gid/lsi)%tiles_x, ty=(gid/lsi)/tiles_x;
	const int x=tx*tile_width+lx, y=ty*tile_height+ly;
	const uint n = x+y*def_screen_width;
	float camera_cache[15]; // cache parameters in case the kernel draws more than one shape
	for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
	ray camray = get_camray(x, y, camera_cache);
	const float distance = intersect_cuboid(camray, (float3)(0.0f, 0.0f, 0.0f), (float)def_Nx, (float)def_Ny, (float)def_Nz);
	camray.origin = camray.origin+fmax(distance, 0.0f)*camray.direction;
	ray reflection, transmission; // reflection and transmission
	float reflectivity, transmissivity;
	int pixelcolor = 0;
	if(raytrace_phi(camray, &reflection, &transmission, &reflectivity, &transmissivity, phi, flags, skybox, def_Nx, def_Ny, def_Nz)) {
		pixelcolor = last_ray_reflectivity(reflection, transmission, reflectivity, transmissivity, skybox); // 1 ray pass
		//pixelcolor = raytrace_phi_next_ray(reflection, transmission, reflectivity, transmissivity, phi, flags, skybox); // 2 ray passes
	} else {
		pixelcolor = skybox_color(camray, skybox);
	}
	//if(raytrace_phi_mirror(camray, &reflection, phi, flags, skybox, def_Nx, def_Ny, def_Nz)) { // reflection only
	//	//pixelcolor = skybox_color(reflection, skybox); // 1 ray pass
	//	pixelcolor = raytrace_phi_next_ray_mirror(reflection, phi, flags, skybox); // 2 ray passes
	//} else {
	//	pixelcolor = skybox_color(camray, skybox);
	//}
	bitmap[n] = pixelcolor; // no zbuffer required
}
)+"#endif"+R( // SURFACE

)+"#ifdef PARTICLES"+R(
)+R(kernel void graphics_particles(const global float* particles, const global float* camera, global int* bitmap, global int* zbuffer) {
	const uint n = get_global_id(0);
	if(n>=(uint)def_particles_N) return;
	float camera_cache[15]; // cache parameters in case the kernel draws more than one shape
	for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
	const int c = COLOR_P; // coloring scheme
	const float3 p = (float3)(particles[n], particles[def_particles_N+(ulong)n], particles[2ul*def_particles_N+(ulong)n]);
	//draw_circle(p, 0.5f, c, camera_cache, bitmap, zbuffer);
	draw_point(p, c, camera_cache, bitmap, zbuffer);
}
)+"#endif"+R( // PARTICLES
)+"#endif"+R( // GRAPHICS



);} // ############################################################### end of OpenCL C code #####################################################################