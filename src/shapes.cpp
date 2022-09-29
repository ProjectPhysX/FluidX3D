#include "shapes.hpp"
#include "lbm.hpp"

const float d = 0.8660254f; // sqrt(3)/2

bool sphere(const uint x, const uint y, const uint z, const float3& p, const float r) {
	const float3 t = float3(x, y, z)-p;
	return sq(t.x)+sq(t.y)+sq(t.z)<=sq(r);
}
bool ellipsoid(const uint x, const uint y, const uint z, const float3& p, const float3& r) {
	const float3 t = float3(x, y, z)-p;
	return sq(t.x/r.x)+sq(t.y/r.y)+sq(t.z/r.z)<=1.0f;
}
bool cube(const uint x, const uint y, const uint z, const float3& p, const float l) {
	const float3 t = float3(x, y, z)-p;
	return t.x>=-0.5f*l&&t.x<=0.5f*l && t.y>=-0.5f*l&&t.y<=0.5f*l && t.z>=-0.5f*l&&t.z<=0.5f*l;
}
bool cuboid(const uint x, const uint y, const uint z, const float3& p, const float3& l) {
	const float3 t = float3(x, y, z)-p;
	return t.x>=-0.5f*l.x&&t.x<=0.5f*l.x && t.y>=-0.5f*l.y&&t.y<=0.5f*l.y && t.z>=-0.5f*l.z&&t.z<=0.5f*l.z;
}
bool cylinder(const uint x, const uint y, const uint z, const float3& p, const float3& n, const float r) {
	const float3 t = float3(x, y, z)-p;
	const float sqnt = sq(dot(normalize(n), t));
	const float dist = sq(t.x)+sq(t.y)+sq(t.z)-sqnt;
	return dist<=sq(r) && sqnt<=sq(0.5f*length(n));
}
bool cone(const uint x, const uint y, const uint z, const float3& p, const float3& n, const float r1, const float r2) {
	const float3 t = float3(x, y, z)-p;
	const float nt = dot(normalize(n), t);
	const float l = length(n);
	const float r = r1+(r2-r1)*(nt+0.5f*l)/l;
	const float sqnt = sq(nt);
	const float dist = sq(t.x)+sq(t.y)+sq(t.z)-sqnt;
	return dist<=sq(r) && sqnt<=sq(0.5f*l);
}
bool pipe(const uint x, const uint y, const uint z, const float3& p, const float3& n, const float r) {
	const float3 t = float3(x, y, z)-p;
	const float sqnt = sq(dot(normalize(n), t));
	const float dist = sq(t.x)+sq(t.y)+sq(t.z)-sqnt;
	return dist>=sq(r-d)&&dist<=sq(r+d) && sqnt<=sq(0.5f*length(n));
}
bool conepipe(const uint x, const uint y, const uint z, const float3& p, const float3& n, const float r1, const float r2) {
	const float3 t = float3(x, y, z)-p;
	const float nt = dot(normalize(n), t);
	const float l = length(n);
	const float r = r1+(r2-r1)*(nt+0.5f*l)/l;
	const float sqnt = sq(nt);
	const float dist = sq(t.x)+sq(t.y)+sq(t.z)-sqnt;
	return dist>=sq(r-d)&&dist<=sq(r+d) && sqnt<=sq(0.5f*l);
}
bool triangle(const uint x, const uint y, const uint z, const float3& p0, const float3& p1, const float3& p2) {
	const float3 t = float3(x, y, z)-p0;
	const float3 u=p1-p0, v=p2-p0, normal=normalize(cross(u, v));
	const float tn = dot(t, normal);
	const float3 p = t-tn*normal; // closest point on triangle
	const float uu=dot(u, u), uv=dot(u, v), vv=dot(v, v), d=uu*vv-uv*uv, pu=dot(p, u), pv=dot(p, v);
	const float w0 = (vv*pu-uv*pv)/d; // barycentric coordinates
	const float w1 = (uu*pv-uv*pu)/d;
	const float w2 = 1.0f-w0-w1;
	const bool is_in_triangle = w0>=0.0f&&w0<=1.0f && w1>=0.0f&&w1<=1.0f && w2>=0.0f&&w2<=1.0f;
	const bool touches_plane = 2.0f*fabs(tn)<=fabs(normal.x)+fabs(normal.y)+fabs(normal.z);
	return is_in_triangle&&touches_plane;
}
bool plane(const uint x, const uint y, const uint z, const float3& p, const float3& n) {
	const float3 t = float3(x, y, z)-p;
	return dot(t, n)<=0.0f;
}
bool torus_x(const uint x, const uint y, const uint z, const float3& p, const float r, const float R) {
	const float3 t = float3(x, y, z)-p;
	const float sqtR = sq(t.x)+sq(t.y)+sq(t.z)+sq(R);
	const float sqr1 = sq(sqtR-sq(r-d));
	const float sqr2 = sq(sqtR-sq(r+d));
	const float sqr3 = 4.0f*sq(R)*(sq(t.y)+sq(t.z));
	return sqr3<=sqr1 && sqr3>=sqr2;
}
bool torus_y(const uint x, const uint y, const uint z, const float3& p, const float r, const float R) {
	const float3 t = float3(x, y, z)-p;
	const float sqtR = sq(t.x)+sq(t.y)+sq(t.z)+sq(R);
	const float sqr1 = sq(sqtR-sq(r-d));
	const float sqr2 = sq(sqtR-sq(r+d));
	const float sqr3 = 4.0f*sq(R)*(sq(t.x)+sq(t.z));
	return sqr3<=sqr1 && sqr3>=sqr2;
}
bool torus_z(const uint x, const uint y, const uint z, const float3& p, const float r, const float R) {
	const float3 t = float3(x, y, z)-p;
	const float sqtR = sq(t.x)+sq(t.y)+sq(t.z)+sq(R);
	const float sqr1 = sq(sqtR-sq(r-d));
	const float sqr2 = sq(sqtR-sq(r+d));
	const float sqr3 = 4.0f*sq(R)*(sq(t.x)+sq(t.y));
	return sqr3<=sqr1 && sqr3>=sqr2;
}

float sphere_plic(const uint x, const uint y, const uint z, const float3& p, const float r) { // sphere with PLIC fill levels returned
	const float3 t = float3(x, y, z)-p;
	const float rad = length(t);
	const float3 n = normalize(t); // sphere surface normal
	const float k = 0.5f*(fabs(n.x)+fabs(n.y)+fabs(n.z));
	/**/ if(rad> r+k) return -1.0f; // point is outside of sphere
	else if(rad<=r-k) return  1.0f; // point is completely inside of sphere
	else return plic_cube_inverse(clamp(r-rad, -k, k), n); // point is at the interface, return PLIC fill level
}
float ellipsoid_plic(const uint x, const uint y, const uint z, const float3& p, const float3& r) { // sphere with PLIC fill levels returned
	const float3 t = float3(x, y, z)-p;
	const float rad1 = (fabs(t.x)*r.x+fabs(t.y)*r.y+fabs(t.z)*r.z)/length(r);
	const float rad2 = rad1*sqrt(sq(t.x/r.x)+sq(t.y/r.y)+sq(t.z/r.z));
	const float3 n = normalize(float3(t.x/sq(r.x), t.y/sq(r.y), t.z/sq(r.z))); // ellipsoid surface normal
	const float k = 0.5f*(fabs(n.x)+fabs(n.y)+fabs(n.z));
	/**/ if(rad1-rad2<-k) return -1.0f; // point is outside of sphere
	else if(rad1-rad2> k) return  1.0f; // point is completely inside of sphere
	else return plic_cube_inverse(clamp(rad1-rad2, -k, k), n); // point is at the interface, return PLIC fill level
}
float cylinder_x_plic(const uint x, const uint y, const uint z, const float3& p, const float r, const float l) { // bar with PLIC fill levels returned
	const float3 t = float3(x, y, z)-p;
	const float rad = sqrt(sq(t.y)+sq(t.z));
	const float3 n = normalize(float3(0.0f, t.y, t.z)); // sphere surface normal
	const float k = 0.5f*(fabs(n.y)+fabs(n.z));
	/**/ if(rad> r+k || t.x< -l-0.5f||t.x> l+0.5f) return -1.0f; // point is outside of bar
	else if(rad<=r-k && t.x>=-l+0.5f&&t.x<=l-0.5f) return  1.0f; // point is completely inside of bar
	else if(t.x<-l+0.5f) return t.x+l+0.5f;
	else if(t.x> l-0.5f) return t.x-l-0.5f;
	else return plic_cube_inverse(clamp(r-rad, -k, k), n);
}
float cylinder_y_plic(const uint x, const uint y, const uint z, const float3& p, const float r, const float l) { // bar with PLIC fill levels returned
	const float3 t = float3(x, y, z)-p;
	const float rad = sqrt(sq(t.x)+sq(t.z));
	const float3 n = normalize(float3(t.x, 0.0f, t.z)); // sphere surface normal
	const float k = 0.5f*(fabs(n.x)+fabs(n.z));
	/**/ if(rad> r+k || t.y< -l-0.5f||t.y> l+0.5f) return -1.0f; // point is outside of bar
	else if(rad<=r-k && t.y>=-l+0.5f&&t.y<=l-0.5f) return  1.0f; // point is completely inside of bar
	else if(t.y<-l+0.5f) return t.y+l+0.5f;
	else if(t.y> l-0.5f) return t.y-l-0.5f;
	else return plic_cube_inverse(clamp(r-rad, -k, k), n);
}
float cylinder_z_plic(const uint x, const uint y, const uint z, const float3& p, const float r, const float l) { // bar with PLIC fill levels returned
	const float3 t = float3(x, y, z)-p;
	const float rad = sqrt(sq(t.x)+sq(t.y));
	const float3 n = normalize(float3(t.x, t.y, 0.0f)); // sphere surface normal
	const float k = 0.5f*(fabs(n.x)+fabs(n.y));
	/**/ if(rad> r+k || t.z< -l-0.5f||t.z> l+0.5f) return -1.0f; // point is outside of bar
	else if(rad<=r-k && t.z>=-l+0.5f&&t.z<=l-0.5f) return  1.0f; // point is completely inside of bar
	else if(t.z<-l+0.5f) return t.z+l+0.5f;
	else if(t.z> l-0.5f) return t.z-l-0.5f;
	else return plic_cube_inverse(clamp(r-rad, -k, k), n);
}
float plane_plic(const uint x, const uint y, const uint z, const float3& p, const float3& n) { // plane with PLIC fill levels returned
	const float3 t = float3(x, y, z)-p;
	const float3 nn = normalize(n);
	const float d0 = -dot(t, nn);
	const float dl = 0.5f*(fabs(nn.x)+fabs(nn.y)+fabs(nn.z));
	/**/ if(d0<=-dl) return -1.0f;
	else if(d0>= dl) return  1.0f;
	else return plic_cube_inverse(d0, nn);
}

void voxelize_line(LBM& lbm, const float3& p0, const float3& p1, const uchar flag) { // voxelize line
	const float xa=p0.x+0.5f, ya=p0.y+0.5f, za=p0.z+0.5f; // start point
	const float xb=p1.x+0.5f, yb=p1.y+0.5f, zb=p1.z+0.5f; // end point
	const int dx=(int)sign(xb-xa), dy=(int)sign(yb-ya), dz=(int)sign(zb-za); // fast ray-grid-traversal
	const float fxa=xa-floor(xa), fya=ya-floor(ya), fza=za-floor(za);
	int3 xyz = int3(floor(xa), floor(ya), floor(za));
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
		if(tmx>1.000001f&&tmy>1.000001f&&tmz>1.000001f) return; // terminate at end of ray
		if(xyz.x<0||xyz.y<0||xyz.z<0||xyz.x>=(int)lbm.get_Nx()||xyz.y>=(int)lbm.get_Ny()||xyz.z>=(int)lbm.get_Nz()) return; // terminate if out of box
		const uint n = lbm.index((uint)xyz.x, (uint)xyz.y, (uint)xyz.z);
		lbm.flags[n] = flag;
	}
}
void voxelize_triangle(LBM& lbm, const float3& p0, const float3& p1, const float3& p2, const uchar flag) { // voxelize triangle
	const float3 u=p1-p0, v=p2-p1, w=p0-p2;
	const float ul=length(u), vl=length(v), wl=length(w);
	for(float i=0.0f; i<ul; i+=0.25f) voxelize_line(lbm, p0+(i/ul)*u, p2, flag);
	for(float i=0.0f; i<vl; i+=0.25f) voxelize_line(lbm, p1+(i/vl)*v, p0, flag);
	for(float i=0.0f; i<wl; i+=0.25f) voxelize_line(lbm, p2+(i/wl)*w, p1, flag);
}
void voxelize_mesh_hull(LBM& lbm, const Mesh* mesh, const uchar flag) { // voxelize mesh
	for(uint i=0u; i<mesh->triangle_number; i++) voxelize_triangle(lbm, mesh->p0[i], mesh->p1[i], mesh->p2[i], flag); // voxelize mesh
}
void voxelize_stl_hull(LBM& lbm, const string& path, const float3& center, const float3x3& rotation, const float size, const uchar flag) { // read and voxelize binary .stl file
	const Mesh* mesh = read_stl(path, float3(lbm.get_Nx(), lbm.get_Ny(), lbm.get_Nz()), center, rotation, size);
	voxelize_mesh_hull(lbm, mesh, flag);
	delete mesh;
}
void voxelize_stl_hull(LBM& lbm, const string& path, const float3x3& rotation, const float size, const uchar flag) { // read and voxelize binary .stl file (place in box center)
	voxelize_stl_hull(lbm, path, lbm.center(), rotation, size, flag);
}
void voxelize_stl_hull(LBM& lbm, const string& path, const float3& center, const float size, const uchar flag) { // read and voxelize binary .stl file (no rotation)
	voxelize_stl_hull(lbm, path, center, float3x3(1.0f), size, flag);
}
void voxelize_stl_hull(LBM& lbm, const string& path, const float size, const uchar flag) { // read and voxelize binary .stl file (place in box center, no rotation)
	voxelize_stl_hull(lbm, path, lbm.center(), float3x3(1.0f), size, flag);
}