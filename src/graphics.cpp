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

const uint light_sources_N = 100u; // maximal number of light sources
float3 light_sources[light_sources_N]; // coordinates of light sources
uint light_sources_n = 0u; // number of light sources

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
bool is_off_screen(const int x, const int y) { // check if point is off screen
	return x<0||x>=(int)camera.width  ||y<0||y>=(int)camera.height; // entire screen
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
void set_light(const uint i, const float3& position) {
	if(i<light_sources_N) {
		light_sources[i] = position;
		light_sources_n = max(light_sources_n, i+1u);
	}
}
int shading(const int color, const float3& p, const float3& normal, const bool translucent=false) {
	const float nl2 = sq(normal.x)+sq(normal.y)+sq(normal.z); // only one sqrt instead of two
	float br = 0.0f;
	for(uint i=0u; i<light_sources_n; i++) {
		const float3 d = light_sources[i]-p; // direction of light source
		const float dl2 = sq(d.x)+sq(d.y)+sq(d.z);
		const float bri = dot(d, normal)/sqrt(nl2*dl2);
		br = fmax(br, translucent ? fabs(bri) : bri);
	}
	return color_mul(color, fmax(1.25f*br, 0.3f));
}
ulong get_font_pixels(const int character) {
	ulong pixels[224] = { // font data (my own 6x11 monospace font)
		0x0000000000000000ull, 0x2082082080082000ull, 0x5145000000000000ull, 0x514F94514F945000ull, 0x21CAA870AA9C2000ull, 0x4AA50421052A9000ull, 0x2145085628A27400ull, 0x2082000000000000ull,
		0x0842082082040800ull, 0x4081041041084000ull, 0xA9C72A0000000000ull, 0x000208F882000000ull, 0x0000000000041080ull, 0x000000F800000000ull, 0x00000000000C3000ull, 0x0821042104208000ull,
		0x7228A28A28A27000ull, 0x0862928820820800ull, 0x722882108420F800ull, 0x7228823028A27000ull, 0x10C51493E1041000ull, 0xFA083C082082F000ull, 0x722820F228A27000ull, 0xF821042084104000ull,
		0x7228A27228A27000ull, 0x7228A27820A27000ull, 0x00030C00030C0000ull, 0x00030C00030C1080ull, 0x0021084081020000ull, 0x00003E03E0000000ull, 0x0204081084200000ull, 0x7228842080082000ull,
		0x7A1B75D75D6E81E0ull, 0x20851453E8A28800ull, 0xF228A2F228A2F000ull, 0x7228208208227000ull, 0xF228A28A28A2F000ull, 0xFA0820F20820F800ull, 0xFA0820F208208000ull, 0x722820B228A27000ull,
		0x8A28A2FA28A28800ull, 0xF88208208208F800ull, 0xF820820820A27000ull, 0x8A4928C289248800ull, 0x820820820820F800ull, 0x8B6DAAAA28A28800ull, 0x8B2CAAAAA9A68800ull, 0xFA28A28A28A2F800ull,
		0xF228A2F208208000ull, 0x7228A28A28A46800ull, 0xF228A2F289228800ull, 0x7228207020A27000ull, 0xF882082082082000ull, 0x8A28A28A28A27000ull, 0x8A28A25145082000ull, 0x8A28A28AAA945000ull,
		0x8A25142145228800ull, 0x8A25142082082000ull, 0xF82104210420F800ull, 0x3882082082083800ull, 0x8204102041020800ull, 0xE08208208208E000ull, 0x2148800000000000ull, 0x000000000000F800ull,
		0x2040000000000000ull, 0x0000007027A27000ull, 0x820820F228A2F000ull, 0x0000007A08207800ull, 0x0820827A28A27800ull, 0x000000722FA07000ull, 0x1882087082082000ull, 0x0000007A28A27827ull,
		0x820820F228A28800ull, 0x0002002082082000ull, 0x0002002082082284ull, 0x4104104946144800ull, 0x4104104104103000ull, 0x000000F2AAAAA800ull, 0x000000F228A28800ull, 0x0000007228A27000ull,
		0x0000007124927104ull, 0x0000007249247041ull, 0x0000007124104000ull, 0x0000003907827000ull, 0x2082087082081800ull, 0x0000008A28A27800ull, 0x0000008A25142000ull, 0x0000008AA7145000ull,
		0x0000008942148800ull, 0x0000008A25142084ull, 0x000000F84210F800ull, 0x0841042041040800ull, 0x2082082082082000ull, 0x8104102104108000ull, 0x000010A840000000ull, 0x0000000000000000ull,
		0x7228208208227084ull, 0x0005008A28A27800ull, 0x004200722FA07000ull, 0x0085007027A27000ull, 0x0005007027A27000ull, 0x0102007027A27000ull, 0x2142007027A27000ull, 0x0000007A08207884ull,
		0x008500722FA07000ull, 0x000500722FA07000ull, 0x010200722FA07000ull, 0x0005002082082000ull, 0x0085002082082000ull, 0x0102002082082000ull, 0xA8851453E8A28800ull, 0x21421453E8A28800ull,
		0x108FA083C820F800ull, 0x000000D89FE85C00ull, 0x5D4B2CFA49249C00ull, 0x0085007228A27000ull, 0x0005007228A27000ull, 0x0102007228A27000ull, 0x0085008A28A27800ull, 0x0102008A28A27800ull,
		0x0005008A25142084ull, 0x500FA28A28A2F800ull, 0x5008A28A28A27000ull, 0x002726AB27200000ull, 0x390410F10410F800ull, 0x082726AB27208000ull, 0x0008942148800000ull, 0x188208708208C000ull,
		0x0042007027A27000ull, 0x0042002082082000ull, 0x0042007228A27000ull, 0x0042008A28A27800ull, 0x42A1007124924800ull, 0x42A122CAAAA68800ull, 0x604624700F800000ull, 0x7228A2700F800000ull,
		0x2080082108A27000ull, 0x01E86DAEDAE17800ull, 0x000000F820000000ull, 0x4B25142164A4B800ull, 0x4B25142125AF8800ull, 0x2080082082082000ull, 0x0002529122400000ull, 0x0009122529000000ull,
		0x9004802409004800ull, 0xA95222568095A120ull, 0xA95A95A95A95A950ull, 0x2082082082082082ull, 0x2082082382082082ull, 0x108008514FA28800ull, 0x214008514FA28800ull, 0x408008514FA28800ull,
		0x01E867A699E17800ull, 0x514514D04D145145ull, 0x5145145145145145ull, 0x000000F04D145145ull, 0x514514D04F000000ull, 0x08472AA2A7108000ull, 0x8A251421C21C2000ull, 0x0000000382082082ull,
		0x20820820F0000000ull, 0x20820823F0000000ull, 0x00000003F2082082ull, 0x20820820F2082082ull, 0x00000003F0000000ull, 0x20820823F2082082ull, 0x42A1007027A27000ull, 0x42A108514FA28800ull,
		0x5145145D07C00000ull, 0x0000007D05D45145ull, 0x514514DC0FC00000ull, 0x000000FC0DD45145ull, 0x5145145D05D45145ull, 0x000000FC0FC00000ull, 0x514514DC0DD45145ull, 0x0227228A27220000ull,
		0xA10A047228A27000ull, 0x712492E924927000ull, 0x214FA083C820F800ull, 0x500FA083C820F800ull, 0x408FA083C820F800ull, 0x000000608208F800ull, 0x10803E208208F800ull, 0x21403E208208F800ull,
		0x01403E208208F800ull, 0x2082082380000000ull, 0x00000000F2082082ull, 0xFFFFFFFFFFFFFFF0ull, 0x00000003FFFFFFF0ull, 0x2082080002082080ull, 0x40803E208208F800ull, 0xFFFFFFFC00000000ull,
		0x10803E8A28A2F800ull, 0x3124944925104000ull, 0x21403E8A28A2F800ull, 0x40803E8A28A2F800ull, 0x42A1007228A27000ull, 0x42A13E8A28A2F800ull, 0x0000008A28A2FA08ull, 0x82082CCA28A2F208ull,
		0x000820F228A2F208ull, 0x1088A28A28A2F800ull, 0x2140228A28A2F800ull, 0x4088A28A28A2F800ull, 0x0042008A25142084ull, 0x12A8945082082000ull, 0xF800000000000000ull, 0x1080000000000000ull,
		0x0000007000000000ull, 0x00823E20803E0000ull, 0x00000000003E03E0ull, 0xC49C8AD0A3974800ull, 0x7BAEBA6820A27000ull, 0x31240C49230248C0ull, 0x000200F802000000ull, 0x0000000000002080ull,
		0x7147000000000000ull, 0x0140000000000000ull, 0x0000002000000000ull, 0x2182087000000000ull, 0x6046046000000000ull, 0x6042107000000000ull, 0x000FFFFFFFC00000ull, 0x0000000000000000ull
	};
	return pixels[clamp(character+256*(character<0)-32, 0, 223)];
}
void draw(const int x, const int y, const float z, const int color, const int stereo) {
	const int index=x+y*(int)camera.width, iz=(int)(z*(2147483647.0f/10000.0f));
	if(!is_off_screen(x, y, stereo)&&iz>camera.zbuffer[index]) {
		camera.zbuffer[index] = iz;
		camera.bitmap[index] = color; // only draw if point is on screen and first in zbuffer
	}
}
void convert_pixel(const float3& p, const int color, const int stereo) {
	int rx, ry; float rz;
	if(convert(rx, ry, rz, p, stereo) && !is_off_screen(rx, ry, stereo)) {
		draw(rx, ry, rz, color, stereo);
	}
}
void convert_circle(const float3& p, const float r, const int color, const int stereo) {
	int rx, ry; float rz;
	if(convert(rx, ry, rz, p, stereo)) {
		const float rs = camera.zoom*camera.dis/(camera.dis-rz*camera.zoom);
		const int radius = (int)(rs*r+0.5f);
		switch(stereo) {
			default: if(rx<                   -radius||rx>=(int)camera.width  +radius || ry<-radius||ry>=(int)camera.height+radius) return; break; // cancel drawing if circle is off screen
			case -1: if(rx<                   -radius||rx>=(int)camera.width/2+radius || ry<-radius||ry>=(int)camera.height+radius) return; break;
			case +1: if(rx<(int)camera.width/2-radius||rx>=(int)camera.width  +radius || ry<-radius||ry>=(int)camera.height+radius) return; break;
		}
		int d=-radius, x=radius, y=0; // Bresenham algorithm for circle
		while(x>=y) {
			draw(rx+x, ry+y, rz, color, stereo);
			draw(rx-x, ry+y, rz, color, stereo);
			draw(rx+x, ry-y, rz, color, stereo);
			draw(rx-x, ry-y, rz, color, stereo);
			draw(rx+y, ry+x, rz, color, stereo);
			draw(rx-y, ry+x, rz, color, stereo);
			draw(rx+y, ry-x, rz, color, stereo);
			draw(rx-y, ry-x, rz, color, stereo);
			d += 2*y+1;
			y++;
			if(d>0) d-=2*(--x);
		}
	}
}
void convert_line(const float3& p0, const float3& p1, const int color, const int stereo) {
	int r0x, r0y, r1x, r1y; float r0z, r1z;
	if(convert(r0x, r0y, r0z, p0, stereo) && convert(r1x, r1y, r1z, p1, stereo)
		&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo) && !intersects_screen(r0x, r0y, r1x, r1y, stereo))) {
		int x=r0x, y=r0y; // Bresenham algorithm
		const float z = 0.5f*(r0z+r1z); // approximate line z position for each pixel to be equal
		const int dx= abs(r1x-r0x), sx=2*(r0x<r1x)-1;
		const int dy=-abs(r1y-r0y), sy=2*(r0y<r1y)-1;
		int err = dx+dy;
		while(x!=r1x||y!=r1y) {
			draw(x, y, z, color, stereo);
			const int e2 = 2*err;
			if(e2>dy) { err+=dy; x+=sx; }
			if(e2<dx) { err+=dx; y+=sy; }
		}
	}
}
void convert_triangle(const float3& p0, const float3& p1, const float3& p2, const int color, const int stereo) {
	int r0x, r0y, r1x, r1y, r2x, r2y; float r0z, r1z, r2z;
	if(convert(r0x, r0y, r0z, p0, stereo) && convert(r1x, r1y, r1z, p1, stereo) && convert(r2x, r2y, r2z, p2, stereo)
		&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo) && is_off_screen(r2x, r2y, stereo)
		&& !intersects_screen(r0x, r0y, r1x, r1y, stereo) && !intersects_screen(r1x, r1y, r2x, r2y, stereo) && !intersects_screen(r2x, r2y, r0x, r0y, stereo))) {
		if(r0y==r1y&&r0y==r2y) return; // return for degenerate triangles
		if(r0y>r1y) { std::swap(r0x, r1x); std::swap(r0y, r1y); } // sort vertices ascending by y
		if(r0y>r2y) { std::swap(r0x, r2x); std::swap(r0y, r2y); }
		if(r1y>r2y) { std::swap(r1x, r2x); std::swap(r1y, r2y); }
		const float z = (r0z+r1z+r2z)/3.0f; // approximate triangle z position for each pixel to be equal
		for(int y=r0y; y<r1y; y++) { // Bresenham algorithm (lower triangle half)
			const int xA = r0x+(r2x-r0x)*(y-r0y)/(r2y-r0y);
			const int xB = r0x+(r1x-r0x)*(y-r0y)/(r1y-r0y);
			for(int x=min(xA, xB); x<max(xA, xB); x++) {
				draw(x, y, z, color, stereo);
			}
		}
		for(int y=r1y; y<r2y; y++) { // Bresenham algorithm (upper triangle half)
			const int xA = r0x+(r2x-r0x)*(y-r0y)/(r2y-r0y);
			const int xB = r1x+(r2x-r1x)*(y-r1y)/(r2y-r1y);
			for(int x=min(xA, xB); x<max(xA, xB); x++) {
				draw(x, y, z, color, stereo);
			}
		}
	}
}
void convert_triangle_interpolated(const float3& p0, const float3& p1, const float3& p2, int c0, int c1, int c2, const int stereo) {
	int r0x, r0y, r1x, r1y, r2x, r2y; float r0z, r1z, r2z;
	if(convert(r0x, r0y, r0z, p0, stereo) && convert(r1x, r1y, r1z, p1, stereo) && convert(r2x, r2y, r2z, p2, stereo)
		&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo) && is_off_screen(r2x, r2y, stereo)
		&& !intersects_screen(r0x, r0y, r1x, r1y, stereo) && !intersects_screen(r1x, r1y, r2x, r2y, stereo) && !intersects_screen(r2x, r2y, r0x, r0y, stereo))) {
		if(r0y==r1y&&r0y==r2y) return; // return for degenerate triangles
		if(r0y>r1y) { std::swap(r0x, r1x); std::swap(r0y, r1y); std::swap(c0, c1); } // sort vertices ascending by y
		if(r0y>r2y) { std::swap(r0x, r2x); std::swap(r0y, r2y); std::swap(c0, c2); }
		if(r1y>r2y) { std::swap(r1x, r2x); std::swap(r1y, r2y); std::swap(c1, c2); }
		const float z = (r0z+r1z+r2z)/3.0f; // approximate triangle z position for each pixel to be equal
		const float d = (float)((r1y-r2y)*(r0x-r2x)+(r2x-r1x)*(r0y-r2y));
		for(int y=r0y; y<r1y; y++) { // Bresenham algorithm (lower triangle half)
			const int xA = r0x+(r2x-r0x)*(y-r0y)/(r2y-r0y);
			const int xB = r0x+(r1x-r0x)*(y-r0y)/(r1y-r0y);
			for(int x=min(xA, xB); x<max(xA, xB); x++) {
				const float w0 = clamp((float)((r1y-r2y)*(x-r2x)+(r2x-r1x)*(y-r2y))/d, 0.0f, 1.0f); // barycentric coordinates
				const float w1 = clamp((float)((r2y-r0y)*(x-r2x)+(r0x-r2x)*(y-r2y))/d, 0.0f, 1.0f);
				const float w2 = clamp(1.0f-w0-w1, 0.0f, 1.0f);
				const int color = color_mix_3(c0, c1, c2, w0, w1, w2); // interpolate color
				draw(x, y, z, color, stereo);
			}
		}
		for(int y=r1y; y<r2y; y++) { // Bresenham algorithm (upper triangle half)
			const int xA = r0x+(r2x-r0x)*(y-r0y)/(r2y-r0y);
			const int xB = r1x+(r2x-r1x)*(y-r1y)/(r2y-r1y);
			for(int x=min(xA, xB); x<max(xA, xB); x++) {
				const float w0 = clamp((float)((r1y-r2y)*(x-r2x)+(r2x-r1x)*(y-r2y))/d, 0.0f, 1.0f); // barycentric coordinates
				const float w1 = clamp((float)((r2y-r0y)*(x-r2x)+(r0x-r2x)*(y-r2y))/d, 0.0f, 1.0f);
				const float w2 = clamp(1.0f-w0-w1, 0.0f, 1.0f);
				const int color = color_mix_3(c0, c1, c2, w0, w1, w2); // interpolate color
				draw(x, y, z, color, stereo);
			}
		}
	}
}
void convert_text(const float3& p, const string& s, const float r, const int color, const int stereo) {
	int rx, ry; float rz;
	if(convert(rx, ry, rz, p, stereo)) {
		const float rs = camera.zoom*camera.dis/(camera.dis-rz*camera.zoom);
		const int radius = (int)(rs*r+0.5f);
		const float tr = fmax(0.85f*radius, 2.0f);
		rx += 4+(int)tr;
		ry += 3+(int)tr;
		for(int i=0u; i<(int)length(s); i++) {
			const int character = (int)s[i];
			const ulong pixels = get_font_pixels(character);
			for(int k=0; k<64; k++) {
				if((pixels>>(63-k))&1) draw(rx+i*6+k%6+(character==113), ry+k/6, rz, color, stereo);
			}
		}
	}
}

void draw_pixel(const int x, const int y, const int color) {
	if(!is_off_screen(x, y)) camera.bitmap[x+y*(int)camera.width] = color; // only draw if point is on screen
}
void draw_circle(const int x, const int y, const int r, const int color) {
	int d=-r, dx=r, dy=0; // Bresenham algorithm for circle
	while(dx>=dy) {
		draw_pixel(x+dx, y+dy, color);
		draw_pixel(x-dx, y+dy, color);
		draw_pixel(x+dx, y-dy, color);
		draw_pixel(x-dx, y-dy, color);
		draw_pixel(x+dy, y+dx, color);
		draw_pixel(x-dy, y+dx, color);
		draw_pixel(x+dy, y-dx, color);
		draw_pixel(x-dy, y-dx, color);
		d += 2*dy+1;
		dy++;
		if(d>0) d-=2*(--dx);
	}
}
void draw_line(const int x0, const int y0, const int x1, const int y1, const int color) {
	int x=x0, y=y0; // Bresenham algorithm
	const int dx= abs(x1-x0), sx=2*(x0<x1)-1;
	const int dy=-abs(y1-y0), sy=2*(y0<y1)-1;
	int err = dx+dy;
	while(x!=x1||y!=y1) {
		draw_pixel(x, y, color);
		const int e2 = 2*err;
		if(e2>dy) { err+=dy; x+=sx; }
		if(e2<dx) { err+=dx; y+=sy; }
	}
}
void draw_triangle(const int x0, const int y0, const int x1, const int y1, const int x2, const int y2, const int color) {
	int r0x=x0, r0y=y0, r1x=x1, r1y=y1, r2x=x2, r2y=y2;
	if(r0y==r1y&&r0y==r2y) return; // return for degenerate triangles
	if(r0y>r1y) { std::swap(r0x, r1x); std::swap(r0y, r1y); } // sort vertices ascending by y
	if(r0y>r2y) { std::swap(r0x, r2x); std::swap(r0y, r2y); }
	if(r1y>r2y) { std::swap(r1x, r2x); std::swap(r1y, r2y); }
	for(int y=r0y; y<r1y; y++) { // Bresenham algorithm (lower triangle half)
		const int xA = r0x+(r2x-r0x)*(y-r0y)/(r2y-r0y);
		const int xB = r0x+(r1x-r0x)*(y-r0y)/(r1y-r0y);
		for(int x=min(xA, xB); x<max(xA, xB); x++) {
			draw_pixel(x, y, color);
		}
	}
	for(int y=r1y; y<r2y; y++) { // Bresenham algorithm (upper triangle half)
		const int xA = r0x+(r2x-r0x)*(y-r0y)/(r2y-r0y);
		const int xB = r1x+(r2x-r1x)*(y-r1y)/(r2y-r1y);
		for(int x=min(xA, xB); x<max(xA, xB); x++) {
			draw_pixel(x, y, color);
		}
	}
}
void draw_rectangle(const int x0, const int y0, const int x1, const int y1, const int color) {
	for(int dy=y0; dy<y1; dy++) { // Bresenham algorithm
		for(int dx=x0; dx<x1; dx++) {
			draw_pixel(dx, dy, color);
		}
	}
}
void draw_text(const int x, const int y, const string& s, const int color) {
	for(int i=0; i<(int)length(s); i++) {
		const int character = (int)s[i];
		const ulong pixels = get_font_pixels(character);
		for(int k=0; k<64; k++) {
			if((pixels>>(63-k))&1) draw_pixel(x+i*6+k%6+(character==113), y+k/6, color);
		}
	}
}
void draw_label(const int x, const int y, const string& s, const int color) {
	draw_text(x, y, s, color);
	if(camera.vr) {
		if(x-camera.width/2>0) {
			draw_text(x-(int)camera.width/2, y, s, color);
		}
		if(x+(int)camera.width/2<(int)camera.width) {
			draw_text(x+(int)camera.width/2, y, s, color);
		}
	}
}
void draw_line_label(const int x0, const int y0, const int x1, const int y1, const int color) {
	draw_line(x0, y0, x1, y1, color);
	if(camera.vr) {
		if((x0+x1)/2-camera.width/2>0) {
			draw_line(x0-(int)camera.width/2, y0, x1-(int)camera.width/2, y1, color);
		}
		if((x0+x1)/2+(int)camera.width/2<(int)camera.width) {
			draw_line(x0+(int)camera.width/2, y0, x1+(int)camera.width/2, y1, color);
		}
	}
}
void draw_bitmap(const int* bitmap) {
	std::copy(bitmap, bitmap+(int)camera.width*(int)camera.height, camera.bitmap);
}

void draw_pixel(const float3& p, const int color) {
	if(!camera.vr) {
		convert_pixel(p, color,  0);
	} else {
		convert_pixel(p, color, -1);
		convert_pixel(p, color, +1);
	}
}
void draw_circle(const float3& p, const float r, const int color) {
	if(!camera.vr) {
		convert_circle(p, r, color,  0);
	} else {
		convert_circle(p, r, color, -1);
		convert_circle(p, r, color, +1);
	}
}
void draw_line(const float3& p0, const float3& p1, const int color) {
	if(!camera.vr) {
		convert_line(p0, p1, color,  0);
	} else {
		convert_line(p0, p1, color, -1);
		convert_line(p0, p1, color, +1);
	}
}
void draw_triangle(const float3& p0, const float3& p1, const float3& p2, const int color, const bool translucent) { // points clockwise from above
	const int cl = shading(color, (p0+p1+p2)/3.0f, cross(p1-p0, p2-p0), translucent);
	if(!camera.vr) {
		convert_triangle(p0, p1, p2, cl,  0);
	} else {
		convert_triangle(p0, p1, p2, cl, -1);
		convert_triangle(p0, p1, p2, cl, +1);
	}
}
void draw_triangle(const float3& p0, const float3& p1, const float3& p2, const int c0, const int c1, const int c2, const bool translucent) { // points clockwise from above
	const float3 normal = cross(p1-p0, p2-p0);
	const int cl0 = shading(c0, p0, normal, translucent);
	const int cl1 = shading(c1, p1, normal, translucent);
	const int cl2 = shading(c2, p2, normal, translucent);
	if(!camera.vr) {
		convert_triangle_interpolated(p0, p1, p2, cl0, cl1, cl2,  0);
	} else {
		convert_triangle_interpolated(p0, p1, p2, cl0, cl1, cl2, -1);
		convert_triangle_interpolated(p0, p1, p2, cl0, cl1, cl2, +1);
	}
}
void draw_text(const float3& p, const float r, const string& s, const int color) {
	if(!camera.vr) {
		convert_text(p, s, r, color,  0);
	} else {
		convert_text(p, s, r, color, -1);
		convert_text(p, s, r, color, +1);
	}
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
		default: camera.input_key(key);
	}
}

#if defined(INTERACTIVE_GRAPHICS)
#if defined(_WIN32)

#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#include <windows.h>
HWND window;
HDC displayDC, memDC;
HBITMAP frameDC;

int key_windows(const int keycode) {
	switch(keycode) {
		case VK_SPACE     : return ' '; // space
		case VK_ESCAPE    : return  27; // escape
		case VK_UP        : return -38; // up arrow
		case VK_DOWN      : return  40; // down arrow
		case VK_LEFT      : return -37; // left arrow
		case VK_RIGHT     : return -39; // right arrow
		case VK_PRIOR     : return -33; // page up
		case VK_NEXT      : return -34; // page down
		case VK_HOME      : return -36; // pos1
		case VK_ADD       : case VK_OEM_PLUS : return '+'; // +
		case VK_SUBTRACT  : case VK_OEM_MINUS: return '-'; // -
		case VK_MULTIPLY  : return '*'; // *
		case VK_DIVIDE    : return '/'; // /
		case VK_OEM_COMMA : return ','; // ,
		case VK_OEM_PERIOD: return '.'; // .
		default: return keycode;
	}
}
void update_frame(const double frametime) {
	main_label(frametime);
	SetBitmapBits(frameDC, 4*(int)camera.width*(int)camera.height, camera.bitmap);
	BitBlt(displayDC, 0, 0, (int)camera.width, (int)camera.height, memDC, 0, 0, SRCCOPY); // copy back buffer to front buffer
	camera.clear_frame(); // clear frame
}
LRESULT CALLBACK WndProc(HWND window, UINT message, WPARAM wParam, LPARAM lParam) {
	if(message==WM_DESTROY) {
		running = false;
		PostQuitMessage(0);
		exit(0);
	} else if(message==WM_MOUSEMOVE) {
		camera.input_mouse_moved((int)LOWORD(lParam), (int)HIWORD(lParam));
		if(!camera.lockmouse) SetCursorPos((int)camera.width/2, (int)camera.height/2); // center cursor
	} else if(message==WM_MOUSEWHEEL) {
		if((short)HIWORD(wParam)>0) camera.input_scroll_up(); else camera.input_scroll_down();
	} else if(message==WM_LBUTTONDOWN||message==WM_MBUTTONDOWN||message==WM_RBUTTONDOWN) {
		if(!camera.lockmouse) {
			ShowCursor(true); // show cursor
		} else {
			ShowCursor(false); // hide cursor
			SetCursorPos(camera.width/2, camera.height/2); // reset cursor
		}
		camera.input_key('U');
	} else if(message==WM_KEYDOWN) {
		int key = key_windows((int)wParam);
		if(key=='U') {
			if(!camera.lockmouse) {
				ShowCursor(true); // show cursor
			} else {
				ShowCursor(false); // hide cursor
				SetCursorPos(camera.width/2, camera.height/2); // reset cursor
			}
		}
		camera.set_key_state(key, true);
		key_bindings(key);
	} else if(message==WM_KEYUP) {
		int key = key_windows((int)wParam);
		camera.set_key_state(key, false);
	}
	return DefWindowProc(window, message, wParam, lParam);
}
#ifdef GRAPHICS_CONSOLE
int main(int argc, char* argv[]) { // call WinMain from dummy main function in order to have an additional console window
	main_arguments = get_main_arguments(argc, argv);
	return WinMain(GetModuleHandle(0), 0, GetCommandLineA(), SW_SHOWMINIMIZED);
}
#endif // GRAPHICS_CONSOLE
INT WINAPI WinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE, _In_ PSTR, _In_ INT) {
	WNDCLASS wndClass = { 0, WndProc, 0, 0, hInstance, LoadIcon(0, IDI_APPLICATION), LoadCursor(0, IDC_ARROW), (HBRUSH)GetStockObject(BLACK_BRUSH), 0, TEXT("WindowClass") };
	RegisterClass(&wndClass);
	MONITORINFO mi = { sizeof(mi) };
	if(!GetMonitorInfo(MonitorFromWindow(window, MONITOR_DEFAULTTONEAREST), &mi)) return 1;
	const uint width  = (uint)(mi.rcMonitor.right-mi.rcMonitor.left); // get screen size, initialize variables
	const uint height = (uint)(mi.rcMonitor.bottom-mi.rcMonitor.top);
	ShowCursor(false); // hide cursor
	SetCursorPos(width/2, height/2);
	DEVMODE lpDevMode = { 0 }; // get monitor fps
	const uint fps_limit = (uint)EnumDisplaySettings(0, ENUM_CURRENT_SETTINGS, &lpDevMode)!=0 ? (uint)lpDevMode.dmDisplayFrequency : 60u; // find out screen refresh rate
	window = CreateWindow("WindowClass", WINDOW_NAME, WS_POPUP|WS_VISIBLE, mi.rcMonitor.left, mi.rcMonitor.top, width, height, 0, 0, hInstance, 0); // create fullscreen window
	displayDC = GetDC(window);
	memDC = CreateCompatibleDC(displayDC);
	frameDC = CreateCompatibleBitmap(displayDC, width, height); // initialize back buffer
	SelectObject(memDC, frameDC);

	camera = Camera(width, height, fps_limit);
	thread compute_thread(main_physics); // start main_physics() in a new thread

	Clock clock;
	double frametime = 1.0;
	MSG msg = { 0 };
	while(msg.message!=WM_QUIT) {
		while(PeekMessage(&msg, 0, 0, 0, PM_REMOVE)) {
			if(msg.message==WM_QUIT) break;
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
		// main loop ################################################################
		camera.update_state();
		main_graphics();
		update_frame(frametime);
		frametime = clock.stop();
		sleep(1.0/(double)camera.fps_limit-frametime);
		clock.start();
		// ##########################################################################
	}

	ReleaseDC(window, displayDC);
	compute_thread.join();
	return 0;
}

#elif defined(__linux__)||defined(__APPLE__)

#include "X11/include/X11/Xlib.h" // sources: libx11-dev, x11proto-dev, libxrandr-dev, libxrender-dev
#include "X11/include/X11/extensions/Xrandr.h" // for multi-monitor

Display* x11_display;
Window x11_window;
GC x11_gc;
XImage* x11_image;
std::atomic_bool updating_frame = true;

int key_linux(const int keycode) {
	switch(keycode) {
		case  38: return  'A'; case 42: return 'G'; case 58: return 'M'; case 39: return 'S'; case 52: return 'Y'; case 15: return '6';
		case  56: return  'B'; case 43: return 'H'; case 57: return 'N'; case 28: return 'T'; case 29: return 'Z'; case 16: return '7';
		case  54: return  'C'; case 31: return 'I'; case 32: return 'O'; case 30: return 'U'; case 10: return '1'; case 17: return '8';
		case  40: return  'D'; case 44: return 'J'; case 33: return 'P'; case 55: return 'V'; case 11: return '2'; case 18: return '9';
		case  26: return  'E'; case 45: return 'K'; case 24: return 'Q'; case 25: return 'W'; case 12: return '3'; case 14: return '5';
		case  41: return  'F'; case 46: return 'L'; case 27: return 'R'; case 53: return 'X'; case 13: return '4'; case 19: return '0';
		case  51: return  '#'; case 94: return '<'; case 35: return '+'; case 61: return '-'; case 59: return ','; case 60: return '.';
		case  65: return  ' '; // space
		case  36: return   10; // enter
		case  22: return    8; // backspace
		case   9: return   27; // escape
		case 107: return  127; // delete
		case  48: return  142; // Ä
		case  47: return  153; // Ö
		case  34: return  154; // Ü
		case  20: return  225; // ß
		case  49: return  248; // ^
		case  21: return  239; // ´
		case 106: return  -45; // insert
		case  77: return -144; // num lock
		case  66: return  -20; // caps lock
		case 115: return  -91; // windows key
		case 117: return  -93; // kontext menu key
		case  87: return  '1'; case  85: return  '6'; // numpad
		case  88: return  '2'; case  79: return  '7'; // numpad
		case  89: return  '3'; case  80: return  '8'; // numpad
		case  83: return  '4'; case  81: return  '9'; // numpad
		case  84: return  '5'; case  90: return  '0'; // numpad
		case  86: return  '+'; case  82: return  '-'; // numpad
		case  63: return  '*'; case 112: return  '/'; // numpad
		case  91: return  ','; case 108: return   10; // numpad
		case  98: return  -38; case 104: return  -40; // up/down arrow
		case 100: return  -37; case 102: return  -39; // left/right arrow
		case  99: return  -33; case 105: return  -34; // page up/down
		case  97: return  -36; case 103: return  -35; // pos1/end
		case  67: return -112; case  73: return -118; // F1/F7
		case  68: return -113; case  74: return -119; // F2/F8
		case  69: return -114; case  75: return -120; // F3/F9
		case  70: return -115; case  76: return -121; // F4/F10
		case  71: return -116; case  95: return -122; // F5/F11
		case  72: return -117; case  96: return -123; // F6/F12
		default: return 0;
	}
}
void update_frame(const double frametime) {
	main_label(frametime);
	updating_frame = true;
	XPutImage(x11_display, x11_window, x11_gc, x11_image, 0, 0, 0, 0, camera.width, camera.height);
	updating_frame = false;
	camera.clear_frame(); // clear frame
	if(!camera.lockmouse) XWarpPointer(x11_display, x11_window, x11_window, 0, 0, camera.width, camera.height, camera.width/2, camera.height/2); // center cursor
}
void hide_cursor() {
	const char zeroes[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	Pixmap x11_pixmap = XCreateBitmapFromData(x11_display, x11_window, zeroes, 8, 8);
	XColor black;
	black.red = black.green = black.blue = 0;
	Cursor x11_invisible_cursor = XCreatePixmapCursor(x11_display, x11_pixmap, x11_pixmap, &black, &black, 0, 0);
	XDefineCursor(x11_display, x11_window, x11_invisible_cursor);
	XFreeCursor(x11_display, x11_invisible_cursor);
	XFreePixmap(x11_display, x11_pixmap);
}
void show_cursor() {
	Cursor x11_cursor = XCreateFontCursor(x11_display, 68);
	XDefineCursor(x11_display, x11_window, x11_cursor);
	XFreeCursor(x11_display, x11_cursor);
}
void input_detection() {
	XEvent x11_event;
	while(running) {
		if(!updating_frame) {
			XNextEvent(x11_display, &x11_event);
			if(x11_event.type==MotionNotify) {
				camera.input_mouse_moved((int)x11_event.xmotion.x, (int)x11_event.xmotion.y);
			} else if(x11_event.type==ButtonPress) { // counterpart: ButtonRelease
				const int x11_button = x11_event.xbutton.button;
				if(x11_button==Button4) { // scroll up
					camera.input_scroll_up();
				} else if(x11_button==Button5) { // scroll down
					camera.input_scroll_down();
				} else {
					if(!camera.lockmouse) {
						show_cursor();
					} else {
						hide_cursor();
						XWarpPointer(x11_display, x11_window, x11_window, 0, 0, camera.width, camera.height, (int)camera.width/2, (int)camera.height/2); // reset cursor
					}
					camera.input_key('U');
				}
			} else if(x11_event.type==KeyPress) {
				const int key = key_linux((int)x11_event.xkey.keycode);
				if(key=='U') {
					if(!camera.lockmouse) {
						show_cursor();
					} else {
						hide_cursor();
						XWarpPointer(x11_display, x11_window, x11_window, 0, 0, camera.width, camera.height, (int)camera.width/2, (int)camera.height/2); // reset cursor
					}
				}
				camera.set_key_state(key, true);
				key_bindings(key);
			} else if(x11_event.type==KeyRelease) {
				const int key = key_linux((int)x11_event.xkey.keycode);
				camera.set_key_state(key, false);
			}
		} else {
			sleep(0.01);
		}
	}
}
int main(int argc, char* argv[]) {
	main_arguments = get_main_arguments(argc, argv);

	x11_display = XOpenDisplay(0);
	if(!x11_display) print_error("No X11 display available.");

	Window x11_root_window = DefaultRootWindow(x11_display);
	XRRScreenResources* x11_screen_resources = XRRGetScreenResources(x11_display, x11_root_window);
	XRROutputInfo* x11_output_info = XRRGetOutputInfo(x11_display, x11_screen_resources, XRRGetOutputPrimary(x11_display, x11_root_window));
	XRRCrtcInfo* x11_crtc_info = XRRGetCrtcInfo(x11_display, x11_screen_resources, x11_output_info->crtc);
	const uint width  = (uint)x11_crtc_info->width; // width and height of primary monitor
	const uint height = (uint)x11_crtc_info->height;
	const int window_offset_x = (int)x11_crtc_info->x; // offset of primary monitor in multi-monitor coordinates
	const int window_offset_y = (int)x11_crtc_info->y;
	XRRFreeCrtcInfo(x11_crtc_info);
	XRRFreeOutputInfo(x11_output_info);
	XRRFreeScreenResources(x11_screen_resources);

	camera = Camera(width, height, 60u);

	x11_window = XCreateWindow(x11_display, x11_root_window, window_offset_x, window_offset_y, width, height, 0, CopyFromParent, CopyFromParent, CopyFromParent, 0, 0);
	XSizeHints x11_size_hints = { PPosition|PSize, window_offset_x, window_offset_y, (int)width, (int)height };
	XSetNormalHints(x11_display, x11_window, &x11_size_hints); // place window on the primary monitor
	struct Hints { long flags=2l, functions=0l, decorations=0b0000000l, input_mode=0l, status=0l; } x11_hints; // decorations=maximize|minimize|menu|title|resize|border|all
	Atom x11_property = XInternAtom(x11_display, "_MOTIF_WM_HINTS", 0);
	XChangeProperty(x11_display, x11_window, x11_property, x11_property, 32, PropModeReplace, (unsigned char*)&x11_hints, 5); // remove window decorations
	XStoreName(x11_display, x11_window, WINDOW_NAME);
	XMapRaised(x11_display, x11_window);
	x11_gc = XCreateGC(x11_display, x11_window, 0, 0);
	x11_image = XCreateImage(x11_display, CopyFromParent, DefaultDepth(x11_display, DefaultScreen(x11_display)), ZPixmap, 0, (char*)camera.bitmap, width, height, 32, 0);
	XSelectInput(x11_display, x11_window, KeyPressMask|KeyReleaseMask|ButtonPressMask|ButtonReleaseMask|PointerMotionMask);
	hide_cursor();
	XWarpPointer(x11_display, None, x11_window, 0, 0, camera.width, camera.height, window_offset_x+(int)camera.width/2, window_offset_y+(int)camera.height/2); // catch cursor from anywhere on monitors

	thread compute_thread(main_physics); // start main_physics() in a new thread
	thread input_thread(input_detection);

	Clock clock;
	double frametime = 1.0;
	while(running) {
		// main loop ################################################################
		camera.update_state();
		main_graphics();
		update_frame(frametime);
		frametime = clock.stop();
		sleep(1.0/(double)camera.fps_limit-frametime);
		clock.start();
		// ##########################################################################
	}

	XFreeGC(x11_display, x11_gc);
	XDestroyWindow(x11_display, x11_window);
	XCloseDisplay(x11_display);
	compute_thread.join();
	input_thread.join();
	return 0;
}

#endif // Linux
#endif // INTERACTIVE_GRAPHICS

#ifdef INTERACTIVE_GRAPHICS_ASCII

uint last_textwidth=0u, last_textheight=0u;
uint fontwidth=8u, fontheight=16u;
void update_frame(const double frametime) {
	if(last_textwidth==0u&&last_textheight==0u) get_console_font_size(fontwidth, fontheight);
	uint textwidth=0u, textheight=0u;
	get_console_size(textwidth, textheight);
	textwidth = max(textwidth, 2u)-1u;
	textwidth = min(textwidth, textheight*camera.width*fontheight/(camera.height*fontwidth));
	textheight = min(textheight, textwidth*camera.height*fontwidth/(camera.width*fontheight));
	textwidth = max(textwidth, 1u);
	textheight = max(textheight, 1u);
	if(textwidth!=last_textwidth||textheight!=last_textheight) {
		clear_console();
		last_textwidth = textwidth;
		last_textheight = textheight;
	}
	show_console_cursor(false);
	const Image image(camera.width, camera.height, camera.bitmap);
	print_video_dither(&image, textwidth, textheight);
	print(alignr(textwidth, to_string(textwidth)+"x"+to_string(textheight)+" "+alignr(4, to_int(1.0/frametime))+"fps"));
	show_console_cursor(true);
	camera.clear_frame(); // clear frame
}
void input_detection() {
	while(running) {
		int key = key_press();
		key -= (key>96&&key<123)*32; // convert lower case to upper case
		key_bindings(key);
	}
}
int main(int argc, char* argv[]) {
	main_arguments = get_main_arguments(argc, argv);
	camera = Camera(384u, 216u, 60u); // width and height must be divisible by 8
	thread compute_thread(main_physics); // start main_physics() in a new thread
	thread input_thread(input_detection);
	Clock clock;
	double frametime = 1.0;
#ifdef UTILITIES_CONSOLE_DITHER_LOOKUP
	print_image_dither_initialize_lookup();
#endif // UTILITIES_CONSOLE_DITHER_LOOKUP
	clear_console();
	while(running) {
		// main loop ################################################################
		camera.update_state();
		main_graphics();
		update_frame(frametime);
		frametime = clock.stop();
		sleep(1.0/(double)camera.fps_limit-frametime);
		clock.start();
		// ##########################################################################
	}
	compute_thread.join();
	input_thread.join();
	return 0;
}

#endif // INTERACTIVE_GRAPHICS_ASCII

#if !defined(INTERACTIVE_GRAPHICS) && !defined(INTERACTIVE_GRAPHICS_ASCII)

int main(int argc, char* argv[]) {
	main_arguments = get_main_arguments(argc, argv);
	camera = Camera(GRAPHICS_FRAME_WIDTH, GRAPHICS_FRAME_HEIGHT, 60u); // width and height must be divisible by 8
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

#endif // no INTERACTIVE_GRAPHICS and no INTERACTIVE_GRAPHICS_ASCII
#endif // GRAPHICS