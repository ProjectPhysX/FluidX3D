#include "setup.hpp"

/*void main_setup() { // 3D Taylor-Green vortices
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(128u, 128u, 128u, 0.01f);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		const float A = 0.25f;
		const uint periodicity = 1u;
		const float a=(float)Nx/(float)periodicity, b=(float)Ny/(float)periodicity, c=(float)Nz/(float)periodicity;
		const float fx = (float)x+0.5f-0.5f*(float)Nx;
		const float fy = (float)y+0.5f-0.5f*(float)Ny;
		const float fz = (float)z+0.5f-0.5f*(float)Nz;
		lbm.u.x[n] =  A*cos(2.0f*pif*fx/a)*sin(2.0f*pif*fy/b)*sin(2.0f*pif*fz/c);
		lbm.u.y[n] = -A*sin(2.0f*pif*fx/a)*cos(2.0f*pif*fy/b)*sin(2.0f*pif*fz/c);
		lbm.u.z[n] =  A*sin(2.0f*pif*fx/a)*sin(2.0f*pif*fy/b)*cos(2.0f*pif*fz/c);
		lbm.rho[n] = 1.0f-sq(A)*3.0f/4.0f*(cos(4.0f*pif*fx/a)+cos(4.0f*pif*fy/b));
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



/*void main_setup() { // 2D Taylor-Green vortices
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(1024u, 1024u, 1u, 0.01f);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		const float A = 0.2f;
		const uint periodicity = 5u;
		const float a=(float)Nx/(float)periodicity, b=(float)Ny/(float)periodicity;
		const float fx = (float)x+0.5f-0.5f*(float)Nx;
		const float fy = (float)y+0.5f-0.5f*(float)Ny;
		lbm.u.x[n] =  A*cos(2.0f*pif*fx/a)*sin(2.0f*pif*fy/b);
		lbm.u.y[n] = -A*sin(2.0f*pif*fx/a)*cos(2.0f*pif*fy/b);
		lbm.rho[n] = 1.0f-sq(A)*3.0f/4.0f*(cos(4.0f*pif*fx/a)+cos(4.0f*pif*fy/b));
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



/*void main_setup() { // Poiseuille flow validation
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint R = 63u; // channel radius (default: 63)
	const float umax = 0.1f; // maximum velocity in channel center (must be < 0.57735027f)
	const float tau = 1.0f; // relaxation time (must be > 0.5f), tau = nu*3+0.5
	const float nu = units.nu_from_tau(tau); // nu = (tau-0.5)/3
	const uint H = 2u*(R+1u);
#ifndef D2Q9
	LBM lbm(H, lcm(sq(H), WORKGROUP_SIZE)/sq(H), H, nu, 0.0f, units.f_from_u_Poiseuille_3D(umax, 1.0f, nu, R), 0.0f); // 3D
#else // D2Q9
	LBM lbm(lcm(H, WORKGROUP_SIZE)/H, H, 1u, nu, units.f_from_u_Poiseuille_2D(umax, 1.0f, nu, R), 0.0f, 0.0f); // 2D
#endif // D2Q9
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
#ifndef D2Q9
		if(!cylinder(x, y, z, lbm.center(), float3(0u, Ny, 0u), 0.5f*(float)min(Nx, Nz)-1.0f)) lbm.flags[n] = TYPE_S; // 3D
#else // D2Q9
		if(y==0u||y==Ny-1u) lbm.flags[n] = TYPE_S; // 2D
#endif // D2Q9
	}	// #########################################################################################################################################################################################
	double error_min = max_double;
	while(true) {
		lbm.run(1000);
		lbm.update_fields();
		lbm.u.read_from_device();
		string s;
		double error_dif=0.0, error_sum=0.0;
#ifndef D2Q9
		for(uint x=0u; x<Nx; x++) {
			for(uint y=Ny/2u; y<Ny/2u+1u; y++) {
				for(uint z=0; z<Nz; z++) {
					const uint n = x+(y+z*Ny)*Nx;
					const double r = (double)sqrt(sq(x+0.5f-0.5f*(float)Nx)+sq(z+0.5f-0.5f*(float)Nz)); // radius from channel center
					if(r<R) {
						const double unum = (double)sqrt(sq(lbm.u.x[n])+sq(lbm.u.y[n])+sq(lbm.u.z[n])); // numerical velocity
						const double uref = umax*(sq(R)-sq(r))/sq(R); // theoretical velocity profile u = G*(R^2-r^2)
						error_dif += sq(unum-uref); // L2 error (Krüger p. 138)
						error_sum += sq(uref);
						s += to_string(r)+" "+to_string(unum)+" "+to_string(uref)+"\n";
					}
				}
			}
		}
#else // D2Q9
		for(uint x=Nx/2u; x<Nx/2u+1u; x++) {
			for(uint y=1u; y<Ny-1u; y++) {
				const uint n = x+(y+0u*Ny)*Nx;
				const double r = (double)(y+0.5f-0.5f*(float)Ny); // radius from channel center
				const double unum = (double)sqrt(sq(lbm.u.x[n])+sq(lbm.u.y[n])); // numerical velocity
				const double uref = umax*(sq(R)-sq(r))/sq(R); // theoretical velocity profile u = G*(R^2-r^2)
				error_dif += sq(unum-uref); // L2 error (Krüger p. 138)
				error_sum += sq(uref);
				s += to_string(r)+" "+to_string(unum)+" "+to_string(uref)+"\n";
			}
		}
#endif // D2Q9
		if(sqrt(error_dif/error_sum)>=error_min) { // stop when error has converged
			print_info("Poiseuille flow error converged after "+to_string(lbm.get_t())+" steps to "+to_string(error_min)); // typical expected L2 errors: 2-5% (Krüger p. 256)
			wait();
			exit(0);
		}
		error_min = fmin(error_min, sqrt(error_dif/error_sum));
		print_info("Poiseuille flow error after t="+to_string(lbm.get_t())+" is "+to_string(error_min)); // typical expected L2 errors: 2-5% (Krüger p. 256)
	}
} /**/



/*void main_setup() { // cylinder in rectangular duct
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const float Re = 25000.0f;
	const float D = 64.0f;
	const float u = rsqrt(3.0f);
	const float w = D;
	const float h = 3.0f*D;
	const float nu = units.nu_from_Re(Re, D, u);
	const float f = units.f_from_u_rectangular_duct(w, D, 1.0f, nu, u);
	LBM lbm(to_uint(w), 12u*to_uint(D), to_uint(h), nu, 0.0f, f, 0.0f);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		lbm.u.y[n] = 0.1f*u;
		if(cylinder(x, y, z, float3(lbm.center().x, 2.0f*D, lbm.center().z), float3(Nx, 0u, 0u), 0.5f*D)) lbm.flags[n] = TYPE_S;
		if(x==0u||x==Nx-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // x and z non periodic
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



/*void main_setup() { // Taylor-Couette flow
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(96u, 96u, 192u, 0.04f);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(!cylinder(x, y, z, lbm.center(), float3(0u, 0u, Nz), (float)(Nx/2u-1u))) lbm.flags[n] = TYPE_S;
		if( cylinder(x, y, z, lbm.center(), float3(0u, 0u, Nz), (float)(Nx/4u   ))) {
			const float3 relative_position = lbm.relative_position(n);
			lbm.u.x[n] =  relative_position.y;
			lbm.u.y[n] = -relative_position.x;
			lbm.u.z[n] = (1.0f-random(2.0f))*0.001f;
			lbm.flags[n] = TYPE_S;
		}
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



/*void main_setup() { // lid-driven cavity
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 96u;
	const float Re = 1000.0f;
	const float u = 0.4f;
	LBM lbm(L, L, L, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(z==Nz-1) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



/*void main_setup() { // delta wing
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 128u;
	const float Re = 10000.0f;
	const float u = 0.1f;
	LBM lbm(L, 4u*L, L, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float3 offset = float3(lbm.center().x, 0.0f, lbm.center().z);
	const float3 p0 = offset+float3(  0*(int)L/64,  5*(int)L/64,  20*(int)L/64);
	const float3 p1 = offset+float3(-20*(int)L/64, 90*(int)L/64, -10*(int)L/64);
	const float3 p2 = offset+float3(+20*(int)L/64, 90*(int)L/64, -10*(int)L/64);
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(triangle(x, y, z, p0, p1, p2)) lbm.flags[n] = TYPE_S;
		else lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	//voxelize_triangle(lbm, p0, p1, p2, TYPE_S);
	lbm.run();
} /**/



void main_setup() { // giving it a shot
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 256u;
	const float Re = 1000000.0f;
	const float u = 0.1f;
	LBM lbm(L, L*3u, L/2u, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float size = 1.75f*(float)L;
	const float3 center = float3(lbm.center().x, 0.52f*size, lbm.center().z+0.03f*size);
	// #############################################################################################################################################################################################
	// QUICK ROTATION FOR OBJECT ALIGNMENT
	// const float3x3 rotation = float3x3(float3(1, 0, 0), radians(-0.0f))*float3x3(float3(0, 0, 1), radians(270.0f))*float3x3(float3(1, 0, 0), radians(0.0f));
	// #############################################################################################################################################################################################
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(-10.0f))*float3x3(float3(0, 0, 1), radians(90.0f))*float3x3(float3(1, 0, 0), radians(90.0f));
	lbm.voxelize_stl(get_exe_path()+"../stl/b52.stl", center, rotation, size);
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	lbm.run();
} /**/



/*void main_setup() { // Boeing 747
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 256u;
	const float Re = 1000000.0f;
	const float u = 0.1f;
	LBM lbm(L, L*2u, L/2u, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float size = 1.0f*(float)L;
	const float3 center = float3(lbm.center().x, 0.55f*size, lbm.center().z);
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(-15.0f));
	lbm.voxelize_stl(get_exe_path()+"../stl/747.stl", center, rotation, size); // https://www.thingiverse.com/thing:2772812/files
	const ulong N=lbm.get_N(); uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//lbm.graphics.set_camera_free(float3(1.0f*(float)Nx, -0.4f*(float)Ny, 2.0f*(float)Nz), -33.0f, 42.0f, 68.0f);
	//Clock clock;
	//lbm.run(0u);
	//while(lbm.get_t()<100000u) {
	//	lbm.graphics.write_frame_png();
	//	lbm.run(28u);
	//}
	//write_file(get_exe_path()+"time.txt", print_time(clock.stop()));
	lbm.run();
} /**/



/*void main_setup() { // Boeing 757
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 912u;
	const float Re = 100000.0f;
	const float u = 0.125f;
	LBM lbm(L, 2u*L, L/2u, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float size = 1.1f*(float)L;
	const float3 center = float3(lbm.center().x, 32.0f+0.5f*size, lbm.center().z);
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(75.0f));
	lbm.voxelize_stl(get_exe_path()+"../stl/757.stl", center, rotation, size); // https://www.thingiverse.com/thing:5091064/files
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//Clock clock;
	//lbm.run(0u);
	//while(lbm.get_t()<100000u) {
	//	lbm.graphics.set_camera_free(float3(1.0f*(float)Nx, -0.4f*(float)Ny, 2.0f*(float)Nz), -33.0f, 42.0f, 68.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/t/");
	//	lbm.graphics.set_camera_free(float3(0.5f*(float)Nx, -0.35f*(float)Ny, -0.7f*(float)Nz), -35.0f, -35.0f, 100.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/b/");
	//	lbm.graphics.set_camera_free(float3(0.0f*(float)Nx, 0.51f*(float)Ny, 0.75f*(float)Nz), 90.0f, 28.0f, 80.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/f/");
	//	lbm.graphics.set_camera_free(float3(0.6f*(float)Nx, -0.15f*(float)Ny, 0.06f*(float)Nz), 0.0f, 0.0f, 100.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/s/");
	//	lbm.run(28u);
	//}
	//write_file(get_exe_path()+"time.txt", print_time(clock.stop()));
	lbm.run();
} /**/



/*void main_setup() { // Star Wars X-wing
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 256u;
	const float Re = 100000.0f;
	const float u = 0.125f;
	LBM lbm(L, L*2u, L/2u, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float size = 1.0f*(float)L;
	const float3 center = float3(lbm.center().x, 32.0f+0.5f*size, lbm.center().z);
	const float3x3 rotation = float3x3(float3(0, 0, 1), radians(180.0f));
	voxelize_stl_hull(lbm, get_exe_path()+"../stl/X-wing.stl", center, rotation, size); // https://www.thingiverse.com/thing:353276/files
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//Clock clock;
	//lbm.run(0u);
	//while(lbm.get_t()<50000u) {
	//	lbm.graphics.set_camera_free(float3(1.0f*(float)Nx, -0.4f*(float)Ny, 2.0f*(float)Nz), -33.0f, 42.0f, 68.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/t/");
	//	lbm.graphics.set_camera_free(float3(0.5f*(float)Nx, -0.35f*(float)Ny, -0.7f*(float)Nz), -33.0f, -40.0f, 100.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/b/");
	//	lbm.graphics.set_camera_free(float3(0.0f*(float)Nx, 0.51f*(float)Ny, 0.75f*(float)Nz), 90.0f, 28.0f, 80.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/f/");
	//	lbm.graphics.set_camera_free(float3(0.7f*(float)Nx, -0.15f*(float)Ny, 0.06f*(float)Nz), 0.0f, 0.0f, 100.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/s/");
	//	lbm.run(28u);
	//}
	//write_file(get_exe_path()+"time.txt", print_time(clock.stop()));
	lbm.run();
} /**/



/*std::atomic_bool revoxelizing = false;
void revoxelize(LBM* lbm, Mesh* mesh) { // voxelize new frames in detached thread in parallel while LBM is running
	for(uint n=0u; n<lbm->get_N(); n++) lbm->flags[n] &= ~TYPE_S; // clear flags
	const float3x3 rotation = float3x3(float3(0.2f, 1.0f, 0.1f), radians(0.4032f)); // create rotation matrix to rotate mesh
	mesh->rotate(rotation); // rotate mesh
	voxelize_mesh_hull(*lbm, mesh, TYPE_S); // voxelize rotated mesh in lbm.flags
	revoxelizing = false; // indicate new voxelizer thread has finished
}
void main_setup() { // Star Wars TIE fighter
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 256u;
	const float Re = 100000.0f;
	const float u = 0.125f;
	LBM lbm(L, L*2u, L, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float size = 0.65f*(float)L;
	const float3 center = float3(lbm.center().x, 0.6f*size, lbm.center().z);
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(90.0f));
	Mesh* mesh = read_stl(get_exe_path()+"../stl/TIE-fighter.stl", lbm.size(), center, rotation, size); // https://www.thingiverse.com/thing:2919109/files
	voxelize_mesh_hull(lbm, mesh, TYPE_S);
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//Clock clock;
	lbm.run(0u);
	while(lbm.get_t()<50000u) {
	//	lbm.graphics.set_camera_free(float3(1.0f*(float)Nx, -0.4f*(float)Ny, 0.63f*(float)Nz), -33.0f, 33.0f, 80.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/t/");
	//	lbm.graphics.set_camera_free(float3(0.3f*(float)Nx, -1.5f*(float)Ny, -0.45f*(float)Nz), -83.0f, -10.0f, 40.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/b/");
	//	lbm.graphics.set_camera_free(float3(0.0f*(float)Nx, 0.57f*(float)Ny, 0.7f*(float)Nz), 90.0f, 29.5f, 80.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/f/");
	//	lbm.graphics.set_camera_free(float3(2.5f*(float)Nx, 0.0f*(float)Ny, 0.0f*(float)Nz), 0.0f, 0.0f, 50.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/s/");
		while(revoxelizing.load()) sleep(0.01f); // wait for voxelizer thread to finish
		lbm.flags.write_to_device(); // lbm.flags on host is finished, write to device now
		revoxelizing = true; // indicate new voxelizer thread is starting
		thread voxelizer(revoxelize, &lbm, mesh); // start new voxelizer thread
		voxelizer.detach(); // detatch voxelizer thread so LBM can run in parallel
		lbm.run(28u); // run LBM in parallel while CPU is voxelizing the next frame
	}
	//write_file(get_exe_path()+"time.txt", print_time(clock.stop()));
} /**/



/*void main_setup() { // Star Trek Enterprise-E
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 528u;
	const float Re = 100000.0f;
	const float u = 0.1f;
	LBM lbm(L, L*3u, L/2u, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float size = 2.5f*(float)L;
	const float3 center = float3(lbm.center().x, 16.0f+0.5f*size, lbm.center().z);
	const float3x3 rotation = float3x3(float3(0, 0, 1), radians(90.0f));
	lbm.voxelize_stl(get_exe_path()+"../stl/Enterprise-E.stl", center, rotation, size); // https://www.thingiverse.com/thing:1423364/files
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//lbm.graphics.set_camera_centered(-60.0f, 20.0f, 0.0f, 2.5f);
	//Clock clock;
	//lbm.run(0u);
	//while(lbm.get_t()<18000u) {
	//	lbm.graphics.write_frame_png();
	//	lbm.run(30u);
	//}
	//write_file(get_exe_path()+"time.txt", print_time(clock.stop()));
	lbm.run();
} /**/



/*void main_setup() { // F1 car
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 256u; // 2152u on 8x MI200
	const float kmh = 100.0f;
	const float si_u = kmh/3.6f;
	const float si_x = 2.0f;
	const float si_rho = 1.225f;
	const float si_nu = 1.48E-5f;
	const float Re = units.si_Re(si_x, si_u, si_nu);
	print_info("Re = "+to_string(Re));
	const float u = 0.08f;
	const float size = 1.6f*(float)L;
	units.set_m_kg_s(size*2.0f/5.5f, u, 1.0f, si_x, si_u, si_rho);
	const float nu = units.nu(si_nu);
	print_info("1s = "+to_string(units.t(1.0f)));
	LBM lbm(L, L*2u, L/2u, nu);
	// #############################################################################################################################################################################################
	const float3 center = float3(lbm.center().x, 0.525f*size, 0.116f*size);
	lbm.voxelize_stl(get_exe_path()+"../stl/Ferrari_SF71H_V5.stl", center, size); // https://www.thingiverse.com/thing:2990512/files
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==Nz-1u) lbm.flags[n] = TYPE_E;
		if(z==0u) lbm.flags[n] = TYPE_S;
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//Clock clock;
	//lbm.run(0u);
	//while(lbm.get_t()<=units.t(1.0f)) {
	//	lbm.graphics.set_camera_free(float3(0.779346f*(float)Nx, -0.315650f*(float)Ny, 0.329444f*(float)Nz), -27.0f, 19.0f, 100.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/a/");
	//	lbm.graphics.set_camera_free(float3(0.556877f*(float)Nx, 0.228191f*(float)Ny, 1.159613f*(float)Nz), 19.0f, 53.0f, 100.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/b/");
	//	lbm.graphics.set_camera_free(float3(0.220650f*(float)Nx, -0.589529f*(float)Ny, 0.085407f*(float)Nz), -72.0f, 21.0f, 86.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/c/");
	//	lbm.run(units.t(0.5f/600.0f)); // run LBM in parallel while CPU is voxelizing the next frame
	//}
	//write_file(get_exe_path()+"time.txt", print_time(clock.stop()));
	lbm.run();
} /**/



#ifdef SURFACE



/*void main_setup() { // hydraulic jump
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(96u, 352u, 96u, 0.007f, 0.0f, 0.0f, -0.0005f);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		const uint H2 = Nz*3u/5u;
		const uint H1 = Nz*2u/5u;
		const uint P1 = Ny*1u/20u;
		const uint P3 = Ny*3u/20u;
		if(z<H2) lbm.flags[n] = TYPE_F;
		if(y<P3&&z< H1) lbm.flags[n] = TYPE_S;
		if(y<P1&&z>=H2) lbm.flags[n] = TYPE_S;
		if(y==1u&&z>=H1&&z<H2) {
			lbm.flags[n] = TYPE_E;
			lbm.rho[n] = 1.55f;
		}
		if(y==Ny-2u) {
			lbm.flags[n] = TYPE_E;
			lbm.u.y[n] = 0.2f/5.0f;
			lbm.rho[n] = 0.99f;
		}
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



/*void main_setup() { // dam break
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(128u, 256u, 256u, 0.005f, 0.0f, 0.0f, -0.0002f, 0.0001f);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(z<Nz*6u/8u && y<Ny/8u) lbm.flags[n] = TYPE_F;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



/*void main_setup() { // swirl
	const uint D = 128u;
	const float f = 0.000001f;
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(D, D, D, 0.005f, 0.0f, 0.0f, -f);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		const uint H = D*5u/6u;
		const uint R = D/4u-1u;
		if(z<H) {
			lbm.flags[n] = TYPE_F;
			lbm.rho[n] = units.rho_hydrostatic(f, (float)z, (float)H);
		}
		if(cylinder(x, y, z, float3(lbm.center().x, lbm.center().y, 0.0f), float3(0u, 0u, 1u), (float)R)) {
			lbm.u.x[n] =  ((float)y+0.5f-0.5f*(float)Ny)/(float)R*0.25f;
			lbm.u.y[n] = -((float)x+0.5f-0.5f*(float)Nx)/(float)R*0.25f;
		}
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



/*void main_setup() { // liquid gallium on speaker
	const uint L = 128u;
	const float rho = 1.0f;
	const float u = 0.09f; // peak velocity of speaker membrane
	const float nu = 0.01f; // make smaller
	const float sigma = 0.005f; // make larger
	const float f = 0.0005f; // make smaller
	const float frequency = 0.01f; // amplitude = u/(2.0f*pif*frequency);
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(L, L, L*3u/4u, nu, 0.0f, 0.0f, -f, sigma);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(z<Nz/3u && x>0u&&x<Nx-1u&&y>0u&&y<Ny-1u&&z>0u&&z<Nz-1u) {
			lbm.rho[n] = units.rho_hydrostatic(f, (float)z, (float)(Nz/3u));
			lbm.u.x[n] = random_symmetric(1E-9f);
			lbm.u.y[n] = random_symmetric(1E-9f);
			lbm.u.z[n] = random_symmetric(1E-9f);
			lbm.flags[n] = TYPE_F;
		}
		if(z==0u) lbm.u.z[n] = 1E-16f;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}	// #########################################################################################################################################################################################
	lbm.run(0u);
	while(running) {
		const float uz = u*sin(2.0f*pif*frequency*(float)lbm.get_t());
		for(uint z=0u; z<1u; z++) {
			for(uint y=1u; y<Ny-1u; y++) {
				for(uint x=1u; x<Nx-1u; x++) {
					const uint n = x+(y+z*Ny)*Nx;
					lbm.u.z[n] = uz;
				}
			}
		}
		lbm.u.write_to_device(2u*N, Nx*Ny); // here faster than lbm.u.write_to_device_3d(1u, Nx-1u, 1u, Ny-1u, 0u, 1u, Nx, Ny, Nz, 2);
		lbm.run(1u);
	}
} /**/



/*void main_setup() { // beach
	const float f = 0.001f; // make smaller
	const float u = 0.12f; // peak velocity of speaker membrane
	const float frequency = 0.0007f; // amplitude = u/(2.0f*pif*frequency);
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(128u, 640u, 96u, 0.01f, 0.0f, 0.0f, -f);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		const uint H = Nz/2u;
		if(z<H) {
			lbm.flags[n] = TYPE_F;
			lbm.rho[n] = units.rho_hydrostatic(f, (float)z, (float)H);
		}
		if(plane(x, y, z, float3(lbm.center().x, 128.0f, 0.0f), float3(0.0f, -1.0f, 8.0f))) lbm.flags[n] = TYPE_S;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
		if(y==0u && x>0u&&x<Nx-1u&&z>0u&&z<Nz-1u) lbm.flags[n] = TYPE_E;
	}	// #########################################################################################################################################################################################
	lbm.run(0u);
	while(running) {
		lbm.u.read_from_device(N+Nx*Ny, N-2u*Nx*Ny); // read velocity first to prevent flickering of velocity field in graphics
		lbm.u.read_from_device(2u*N+Nx*Ny, N-2u*Nx*Ny);
		const float uy = u*sin(2.0f*pif*frequency*(float)lbm.get_t());
		const float uz = 0.5f*u*cos(2.0f*pif*frequency*(float)lbm.get_t());
		for(uint z=1u; z<Nz-1u; z++) {
			for(uint y=0u; y<1u; y++) {
				for(uint x=1u; x<Nx-1u; x++) {
					const uint n = x+(y+z*Ny)*Nx;
					lbm.u.y[n] = uy;
					lbm.u.z[n] = uz;
				}
			}
		}
		lbm.u.write_to_device_3d(1u, Nx-1u, 0u, 1u, 1u, Nz-1u, Nx, Ny, Nz, 1);
		lbm.u.write_to_device_3d(1u, Nx-1u, 0u, 1u, 1u, Nz-1u, Nx, Ny, Nz, 2);
		lbm.run(100u);
	}
} /**/



/*void main_setup() { // river
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(128u, 384u, 96u, 0.02f, 0.0f, -0.00007f, -0.0005f, 0.01f);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		const int R = 20, H = 32;
		if(z==0) lbm.flags[n] = TYPE_S;
		else if(z<H) {
			lbm.flags[n] = TYPE_F;
			lbm.u.y[n] = -0.1f;
		}
		if(cylinder(x, y, z, float3(Nx*2u/3u, Ny*2u/3u, Nz/2u)+0.5f, float3(0u, 0u, Nz), (float)R)) lbm.flags[n] = TYPE_S;
		if(cuboid(x, y, z, float3(Nx/3u, Ny/3u, Nz/2u)+0.5f, float3(2u*R, 2u*R, Nz))) lbm.flags[n] = TYPE_S;
		if(x==0u||x==Nx-1u) lbm.flags[n] = TYPE_S; // x non periodic
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



/*void main_setup() { // raindrop setup (#define D3Q19, SRT, VOLUME_FORCE, EQUILIBRIUM_BOUNDARIES, SURFACE)
	const int box_diameter = 256; // 936 for Quadro RTX 8000, 1064 for MI200
	float drop_diameter = box_diameter/5;
	const int select_drop_size = 12;
	const float alpha_sim = 20.0f;

	if(drop_diameter==-1.0f) drop_diameter = 0.1f*(float)box_diameter;
	const float scale = (float)box_diameter/(10.0f*drop_diameter); // 256.0f/400.0f; // 1.0f;

	// rain drop parameters from "Effects of Altitude on Maximum Raindrop Size and Fall Velocity as Limited by Collisional Breakup, Fig. 3" in SI-units
	float const si_nu = 1.0508E-6f; // kinematic shear viscosity [m^2/s] at 20°C and 35g/l salinity
	const float si_rho = 1024.8103f; // fluid density [kg/m^3] at 20°C and 35g/l salinity
	const float si_sigma = 73.81E-3f; // fluid surface tension [kg/s^2] at 20°C and 35g/l salinity
	const float si_g = 9.81f; // gravitational acceleration [m/s^2]
	const float alpha = alpha_sim; // impact angle [°], 0 = vertical

	//                            0        1        2        3        4        5        6        7        8        9       10       11       12       13 (13 is for validation)
	const float si_Ds[] = { 1.0E-3f, 1.5E-3f, 2.0E-3f, 2.5E-3f, 3.0E-3f, 3.5E-3f, 4.0E-3f, 4.5E-3f, 5.0E-3f, 5.5E-3f, 6.0E-3f, 6.5E-3f, 7.0E-3f, 4.1E-3f };
	const float si_us[] = {   4.50f,   5.80f,   6.80f,   7.55f,   8.10f,   8.45f,   8.80f,   9.05f,   9.20f,   9.30f,   9.40f,   9.45f,   9.55f,   7.21f };
	const float si_D = si_Ds[select_drop_size]; // drop diameter [m] (1-7mm)
	const float si_u = si_us[select_drop_size]; // impact velocity [m/s] (4.50-9.55m/s)
	const float si_H  =  4.0f*si_D*scale; // liquid pool height [m] (4*D is sufficient for deep pool)
	const float si_Lx = 10.0f*si_D*scale; // simulation box width [m]
	const float si_Lz =  8.5f*si_D*scale; // simulation box height [m]

	// determine a length, a velocity and the mean density in simulation units
	const float Lx = (float)box_diameter; // simulation box width
	const float u = 0.05f; // impact velocity
	const float rho = 1.0f; // density
	units.set_m_kg_s(Lx, u, rho, si_Lx, si_u, si_rho); // calculate 3 independent conversion factors (m, kg, s)
	print_info("D  = "+to_string(si_D, 6u));
	print_info("Re = "+to_string(units.si_Re(si_D, si_u, si_nu), 6u));
	print_info("We = "+to_string(units.si_We(si_D, si_u, si_rho, si_sigma), 6u));
	print_info("Fr = "+to_string(units.si_Fr(si_D, si_u, si_g), 6u));
	print_info("Ca = "+to_string(units.si_Ca(si_u, si_rho, si_nu, si_sigma), 6u));
	print_info("Bo = "+to_string(units.si_Bo(si_D, si_rho, si_g, si_sigma), 6u));
	print_info("10ms = "+to_string(units.t(0.01f))+" LBM steps");
	const float nu = units.nu(si_nu); // calculate values for remaining parameters in simulation units
	const float sigma = units.sigma(si_sigma);
	const float f = units.f(si_rho, si_g); // force per volume
	const float Lz = units.x(si_Lz);
	const float H = units.x(si_H);
	const float R = 0.5f*units.x(si_D); // drop radius
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(to_uint(Lx), to_uint(Lx), to_uint(Lz), nu, 0.0f, 0.0f, -f, sigma); // largest box size on Titan Xp with FP32: 384^2, FP16: 464^3
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		lbm.rho[n] = rho; // set density everywhere
		if(sphere(x, y, z, float3(0.5f*(float)Nx, 0.5f*(float)Ny-2.0f*R*tan(alpha*pif/180.0f), H+R+2.5f)+0.5f, R+2.0f)) {
			const float b = sphere_plic(x, y, z, float3(0.5f*(float)Nx, 0.5f*(float)Ny-2.0f*R*tan(alpha*pif/180.0f)+0.5f, H+R+2.5f), R);
			if(b!=-1.0f) {
				lbm.u.y[n] =  sin(alpha*pif/180.0f)*u;//+random_symmetric(0.1f); // break symmetry by initializing with noise
				lbm.u.z[n] = -cos(alpha*pif/180.0f)*u;//+random_symmetric(0.1f); // break symmetry by initializing with noise
				if(b==1.0f) {
					lbm.flags[n] = TYPE_F;
					lbm.phi[n] = 1.0f;
				} else {
					lbm.flags[n] = TYPE_I;
					lbm.phi[n] = b;
				}
			}
		}
		if(z==0) lbm.flags[n] = TYPE_S;
		else if(z==to_uint(H)) {
			lbm.flags[n] = TYPE_I;
			lbm.phi[n] = 0.5f; // not strictly necessary, but should be clearer (phi is automatically initialized to 0.5f for TYPE_I if not initialized)
		} else if((float)z<H) lbm.flags[n] = TYPE_F;
		else if((x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==Nz-1u)&&(float)z>H+R) { // make drops that hit the simulation box ceiling disappear
			lbm.rho[n] = 0.5f;
			lbm.flags[n] = TYPE_E;
		}
	}	// #########################################################################################################################################################################################
	lbm.run(0u);
	//key_1 = false; // turn off boundary
	//key_6 = true; // turn on surface raytracing
	//Clock clock;

	// image
	//lbm.run(units.t(0.0015f));
	//print_info("compute time: "+print_time(clock.stop()));
	//clock.start();
	//lbm.graphics.set_camera_centered(-30.0f, 20.0f, 100.0f, 1.0f);
	//lbm.graphics.write_frame_png();
	//print_info("render time: "+to_string(clock.stop(), 3u));

	// video
	//while(units.si_t(lbm.get_t())<=0.003f) {
	//	lbm.graphics.set_camera_centered(-30.0f, 20.0f, 100.0f, 1.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/new/");
	//	lbm.graphics.set_camera_centered(10.0f, 40.0f, 100.0f, 1.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/p/");
	//	lbm.graphics.set_camera_centered(0.0f, 0.0f, 45.0f, 1.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/o/");
	//	lbm.graphics.set_camera_centered(0.0f, 90.0f, 45.0f, 1.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/t/");
	//	lbm.run(units.t(0.004f/600u));
	//}
	//print_info("compute + render time: "+to_string(clock.stop(), 3u));

	lbm.run();
} /**/



/*void main_setup() { // bursting bubble setup
	const float d=64.0f;
	const float sigma=0.0003f;
	const float si_nu = 1E-6f; // kinematic shear viscosity (water) [m^2/s]
	const float si_rho = 1E3f; // density (water) [kg/m^3]
	const float si_sigma = 0.072f; // surface tension (water) [kg/s^2]
	const float si_d = 4E-3f; // bubble diameter [m]
	const float si_g = 9.81f; // gravitational acceleration [m/s^2]
	const float si_f = units.si_f_from_si_g(si_g, si_rho);
	const float si_rho_particles = si_rho;
	const float rho = 1.0f;
	const float m = si_d/d; // length si_x = x*[m]
	const float kg = si_rho/rho*cb(m); // density si_rho = rho*[kg/m^3]
	const float s = sqrt(sigma/si_sigma*kg); // velocity si_sigma = sigma*[kg/s^2]
	units.set_m_kg_s(m, kg, s);
	const float f = units.f(si_f);
	const float nu = units.nu(si_nu);
	const uint Lx = to_uint(4.0f*d);
	const uint Ly = to_uint(4.0f*d);
	const uint Lz = d==112.0f ? to_uint(2.5f*d) : to_uint(3.0f*d);
	const uint H = to_uint(2.0f*d);
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(Lx, Ly, Lz, nu, 0.0f, 0.0f, -f, sigma);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(z<H) { // pool
			lbm.flags[n] = TYPE_F;
		}
		const float r = 0.5f*d;
		if(sphere(x, y, z, float3(lbm.center().x, lbm.center().y, (float)H-0.5f*d), r+1.0f)) { // bubble
			const float b = clamp(sphere_plic(x, y, z, float3(lbm.center().x, lbm.center().y, (float)H-0.5f*d), r), 0.0f, 1.0f);
			if(b==1.0f) {
				lbm.flags[n] = TYPE_G;
				lbm.phi[n] = 0.0f;
			} else {
				lbm.flags[n] = TYPE_I;
				lbm.phi[n] = (1.0f-b);
			}
		}
		if(z==0) lbm.flags[n] = TYPE_S;
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



/*void main_setup() { // cube with changing gravity
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(96u, 96u, 96u, 0.02f, 0.0f, 0.0f, -0.001f, 0.001f);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(x<Nx*2u/3u&&y<Ny*2u/3u) lbm.flags[n] = TYPE_F;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S;
	}	// #########################################################################################################################################################################################
	while(running) {
		lbm.set_f(0.0f, 0.0f, -0.001f);
		lbm.run(2500u);
		lbm.set_f(0.0f, +0.001f, 0.0f);
		lbm.run(2500u);
		lbm.set_f(0.0f, 0.0f, +0.001f);
		lbm.run(2500u);
		lbm.set_f(0.0f, -0.001f, 0.0f);
		lbm.run(2000u);
		lbm.set_f(0.0f, 0.0f, 0.0f);
		lbm.run(3000u);
	}
} /**/



/*void main_setup() { // periodic faucet mass conservation test
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(96u, 192u, 128u, 0.02f, 0.0f, 0.0f, -0.001f);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(y>Ny*5/6) lbm.flags[n] = TYPE_F;
		const uint D = max(Nx, Nz);
		const uint r = D/6;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u) lbm.flags[n] = TYPE_S; // x and y non periodic
		if((z==0u||z==Nz-1u) && sq(x-Nx/2)+sq(y-Nx/2)>sq(r)) lbm.flags[n] = TYPE_S; // z non periodic
		if(y<=Nx/2+2*r && torus_x(x, y, z, float3(Nx/2, Nx/2+r, Nz)+0.5f, (float)r, (float)r)) lbm.flags[n] = TYPE_S;
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



/*void main_setup() { // two colliding droplets in force field
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(256u, 256u, 128u, 0.014f, 0.0f, 0.0f, 0.0f, 0.0001f);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(sphere(x, y, z, lbm.center()-float3(0u, 10u, 0u), 32.0f)) {
			lbm.flags[n] = TYPE_F;
			lbm.u.y[n] = 0.025f;
		}
		if(sphere(x, y, z, lbm.center()+float3(30u, 40u, 0u), 12.0f)) {
			lbm.flags[n] = TYPE_F;
			lbm.u.y[n] = -0.2f;
		}
		lbm.F.x[n] = -0.001f*lbm.relative_position(n).x;
		lbm.F.y[n] = -0.001f*lbm.relative_position(n).y;
		lbm.F.z[n] = -0.0005f*lbm.relative_position(n).z;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



#endif // SURFACE
#ifdef TEMPERATURE



/*void main_setup() { // Rayleigh-Benard convection
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(256u, 256u, 64u, 0.02f, 0.0f, 0.0f, -0.001f, 0.0f, 1.0f, 1.0f);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		lbm.u.x[n] = random_symmetric(0.015f);
		lbm.u.y[n] = random_symmetric(0.015f);
		lbm.u.z[n] = random_symmetric(0.015f);
		if(z==1u) {
			lbm.T[n] = 1.75f;
			lbm.flags[n] = TYPE_T;
		} else if(z==Nz-2u) {
			lbm.T[n] = 0.25f;
			lbm.flags[n] = TYPE_T;
		}
		//if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
		if(z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S;
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



/*void main_setup() { // TEMPERATURE test
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(32u, 196u, 60u, 0.02f, 0.0f, 0.0f, -0.001f, 0.0f, 1.0f, 1.0f);
	// #############################################################################################################################################################################################
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(y==1) {
			lbm.T[n] = 1.8f;
			lbm.flags[n] = TYPE_T;
		} else if(y==Ny-2) {
			lbm.T[n] = 0.3f;
			lbm.flags[n] = TYPE_T;
		}
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



#endif // TEMPERATURE
/* #else // BENCHMARK
#include "info.hpp"
void main_setup() { // benchmark
	uint mlups = 0u;
	{ // ######################################################## define simulation box size, viscosity and volume force ###########################################################################
		//LBM lbm( 32u,  32u,  32u, 1.0f);
		//LBM lbm( 48u,  48u,  48u, 1.0f);
		//LBM lbm( 64u,  64u,  64u, 1.0f);
		//LBM lbm( 96u,  96u,  96u, 1.0f);
		//LBM lbm(128u, 128u, 128u, 1.0f);
		//LBM lbm(192u, 192u, 192u, 1.0f);
		LBM lbm(256u, 256u, 256u, 1.0f);
		//LBM lbm(384u, 384u, 384u, 1.0f);
		//LBM lbm(464u, 464u, 464u, 1.0f);
		//LBM lbm(480u, 480u, 480u, 1.0f);
		//LBM lbm(512u, 512u, 512u, 1.0f);
		// #########################################################################################################################################################################################
		for(uint i=0u; i<1000u; i++) {
			lbm.run(10u);
			mlups = max(mlups, to_uint((double)lbm.get_N()*1E-6/info.dt_smooth));
		}
	} // make lbm object go out of scope to free its memory
	print_info("Peak MLUPs/s = "+to_string(mlups));
#if defined(_WIN32)
	wait();
#endif // Windows
} /**/
