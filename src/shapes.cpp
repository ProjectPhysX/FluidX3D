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