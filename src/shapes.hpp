#pragma once

#include "defines.hpp"
#include "utilities.hpp"

bool sphere(const uint x, const uint y, const uint z, const float3& p, const float r);
bool ellipsoid(const uint x, const uint y, const uint z, const float3& p, const float3& r);
bool cube(const uint x, const uint y, const uint z, const float3& p, const float l);
bool cuboid(const uint x, const uint y, const uint z, const float3& p, const float3& l);
bool cylinder(const uint x, const uint y, const uint z, const float3& p, const float3& n, const float r);
bool cone(const uint x, const uint y, const uint z, const float3& p, const float3& n, const float r1, const float r2);
bool pipe(const uint x, const uint y, const uint z, const float3& p, const float3& n, const float r);
bool conepipe(const uint x, const uint y, const uint z, const float3& p, const float3& n, const float r1, const float r2);
bool triangle(const uint x, const uint y, const uint z, const float3& p0, const float3& p1, const float3& p2);
bool plane(const uint x, const uint y, const uint z, const float3& p, const float3& n);
bool torus_x(const uint x, const uint y, const uint z, const float3& p, const float r, const float R);
bool torus_y(const uint x, const uint y, const uint z, const float3& p, const float r, const float R);
bool torus_z(const uint x, const uint y, const uint z, const float3& p, const float r, const float R);

float sphere_plic(const uint x, const uint y, const uint z, const float3& p, const float r); // sphere with PLIC fill levels returned
float ellipsoid_plic(const uint x, const uint y, const uint z, const float3& p, const float3& r); // ellipsoid with PLIC fill levels returned
float cylinder_x_plic(const uint x, const uint y, const uint z, const float3& p, const float r, const float l); // bar with PLIC fill levels returned
float cylinder_y_plic(const uint x, const uint y, const uint z, const float3& p, const float r, const float l); // bar with PLIC fill levels returned
float cylinder_z_plic(const uint x, const uint y, const uint z, const float3& p, const float r, const float l); // bar with PLIC fill levels returned
float plane_plic(const uint x, const uint y, const uint z, const float3& p, const float3& n); // plane with PLIC fill levels returned

class LBM; // forward-declare LBM class
void voxelize_line(LBM& lbm, const float3& p0, const float3& p1, const uchar flag); // voxelize line
void voxelize_triangle(LBM& lbm, const float3& p0, const float3& p1, const float3& p2, const uchar flag); // voxelize triangle

struct Mesh { // triangle mesh
	uint triangle_number = 0u;
	float3 center;
	float3* p0;
	float3* p1;
	float3* p2;
	Mesh(const uint triangle_number, const float3& center) {
		this->triangle_number = triangle_number;
		this->center = center;
		this->p0 = new float3[triangle_number];
		this->p1 = new float3[triangle_number];
		this->p2 = new float3[triangle_number];
	}
	~Mesh() {
		delete[] p0;
		delete[] p1;
		delete[] p2;
	}
	void scale(const float scale) {
		for(uint i=0u; i<triangle_number; i++) {
			p0[i] = scale*(p0[i]-center)+center;
			p1[i] = scale*(p1[i]-center)+center;
			p2[i] = scale*(p2[i]-center)+center;
		}
	}
	void rotate(const float3x3& rotation) {
		for(uint i=0u; i<triangle_number; i++) {
			p0[i] = rotation*(p0[i]-center)+center;
			p1[i] = rotation*(p1[i]-center)+center;
			p2[i] = rotation*(p2[i]-center)+center;
		}
	}
};
Mesh* read_stl(LBM& lbm, const string& path, const float size=-1.0f); // read binary .stl file (place in box center, no rotation)
Mesh* read_stl(LBM& lbm, const string& path, const float3& center, const float size=-1.0f); // read binary .stl file (no rotation)
Mesh* read_stl(LBM& lbm, const string& path, const float3x3& rotation, const float size=-1.0f); // read binary .stl file (place in box center)
Mesh* read_stl(LBM& lbm, const string& path, const float3& center, const float3x3& rotation, const float size=-1.0f); // read binary .stl file
void voxelize_mesh(LBM& lbm, const Mesh* mesh, const uchar flag); // voxelize mesh

void voxelize_stl(LBM& lbm, const string& path, const float size=-1.0f, const uchar flag=TYPE_S); // read and voxelize binary .stl file (place in box center, no rotation)
void voxelize_stl(LBM& lbm, const string& path, const float3& center, const float size=-1.0f, const uchar flag=TYPE_S); // read and voxelize binary .stl file (no rotation)
void voxelize_stl(LBM& lbm, const string& path, const float3x3& rotation, const float size=-1.0f, const uchar flag=TYPE_S); // read and voxelize binary .stl file (place in box center)
void voxelize_stl(LBM& lbm, const string& path, const float3& center, const float3x3& rotation, const float size=-1.0f, const uchar flag=TYPE_S); // read and voxelize binary .stl file