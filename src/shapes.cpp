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
Mesh* read_stl(LBM& lbm, const string& path, const float3& center, const float3x3& rotation, const float size) { // read binary .stl file
	const string filename = create_file_extension(path, ".stl");
	std::ifstream file(filename, std::ios::in|std::ios::binary);
	if(file.fail()) print_error("File \""+filename+"\" does not exist!");
	file.seekg(0, std::ios::end);
	const uint filesize = (uint)file.tellg();
	file.seekg(0, std::ios::beg);
	uchar* data = new uchar[filesize];
	file.read((char*)data, filesize);
	file.close();
	if(filesize==0u) println("\rError: File \""+filename+"\" is corrupt!");
	const uint triangle_number = ((uint*)data)[20];
	uint counter = 84u;
	if(triangle_number>0u&&filesize==84u+50u*triangle_number) print_info("Loading \""+filename+"\" with "+to_string(triangle_number)+" triangles.");
	else print_error("File \""+filename+"\" is corrupt or unsupported! Only binary .stl files are supported.");
	Mesh* mesh = new Mesh(triangle_number, center);
	mesh->p0[0] = float3(0.0f); // to fix warning C6001
	for(uint i=0u; i<triangle_number; i++) {
		const float* triangle_data = (float*)(data+counter);
		counter += 50u;
		mesh->p0[i] = rotation*float3(triangle_data[ 3], triangle_data[ 4], triangle_data[ 5]); // read positions of triangle vertices and rotate them
		mesh->p1[i] = rotation*float3(triangle_data[ 6], triangle_data[ 7], triangle_data[ 8]);
		mesh->p2[i] = rotation*float3(triangle_data[ 9], triangle_data[10], triangle_data[11]);
	}
	delete[] data;
	float3 pmin=mesh->p0[0], pmax=pmin; // auto-rescale mesh
	for(uint i=0u; i<triangle_number; i++) {
		const float3 p0=mesh->p0[i], p1=mesh->p1[i], p2=mesh->p2[i];
		pmin.x = fmin(fmin(fmin(p0.x, p1.x), p2.x), pmin.x);
		pmin.y = fmin(fmin(fmin(p0.y, p1.y), p2.y), pmin.y);
		pmin.z = fmin(fmin(fmin(p0.z, p1.z), p2.z), pmin.z);
		pmax.x = fmax(fmax(fmax(p0.x, p1.x), p2.x), pmax.x);
		pmax.y = fmax(fmax(fmax(p0.y, p1.y), p2.y), pmax.y);
		pmax.z = fmax(fmax(fmax(p0.z, p1.z), p2.z), pmax.z);
	}
	const float3 offset = -0.5f*(pmin+pmax);
	float scale = 1.0f;
	if(size>0) { // rescale to specified size
		scale = size/fmax(fmax(pmax.x-pmin.x, pmax.y-pmin.y), pmax.z-pmin.z);
	} else { // auto-rescale to largest possible size
		const float scale_x = (float)lbm.get_Nx()/(pmax.x-pmin.x);
		const float scale_y = (float)lbm.get_Ny()/(pmax.y-pmin.y);
		const float scale_z = (float)lbm.get_Nz()/(pmax.z-pmin.z);
		scale = fmin(fmin(scale_x, scale_y), scale_z);
	}
	for(uint i=0u; i<triangle_number; i++) { // rescale mesh
		mesh->p0[i] = center+scale*(offset+mesh->p0[i]);
		mesh->p1[i] = center+scale*(offset+mesh->p1[i]);
		mesh->p2[i] = center+scale*(offset+mesh->p2[i]);
	}
	return mesh;
}
Mesh* read_stl(LBM& lbm, const string& path, const float3x3& rotation, const float size) { // read binary .stl file (place in box center)
	return read_stl(lbm, path, lbm.center(), rotation, size);
}
Mesh* read_stl(LBM& lbm, const string& path, const float3& center, const float size) { // read binary .stl file (no rotation)
	return read_stl(lbm, path, center, float3x3(1.0f), size);
}
Mesh* read_stl(LBM& lbm, const string& path, const float size) { // read binary .stl file (place in box center, no rotation)
	return read_stl(lbm, path, lbm.center(), float3x3(1.0f), size);
}
void voxelize_mesh(LBM& lbm, const Mesh* mesh, const uchar flag) { // voxelize mesh
	for(uint i=0u; i<mesh->triangle_number; i++) voxelize_triangle(lbm, mesh->p0[i], mesh->p1[i], mesh->p2[i], flag); // voxelize mesh
}
void voxelize_stl(LBM& lbm, const string& path, const float3& center, const float3x3& rotation, const float size, const uchar flag) { // read and voxelize binary .stl file
	const Mesh* mesh = read_stl(lbm, path, center, rotation, size);
	voxelize_mesh(lbm, mesh, flag);
	delete mesh;
}
void voxelize_stl(LBM& lbm, const string& path, const float3x3& rotation, const float size, const uchar flag) { // read and voxelize binary .stl file (place in box center)
	voxelize_stl(lbm, path, lbm.center(), rotation, size, flag);
}
void voxelize_stl(LBM& lbm, const string& path, const float3& center, const float size, const uchar flag) { // read and voxelize binary .stl file (no rotation)
	voxelize_stl(lbm, path, center, float3x3(1.0f), size, flag);
}
void voxelize_stl(LBM& lbm, const string& path, const float size, const uchar flag) { // read and voxelize binary .stl file (place in box center, no rotation)
	voxelize_stl(lbm, path, lbm.center(), float3x3(1.0f), size, flag);
}