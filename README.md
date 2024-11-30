# Easy to Run FluidX3D
<a href="https://youtu.be/5AzxwQpng0M"><img src="https://img.youtube.com/vi/5AzxwQpng0M/maxresdefault.jpg" alt="10 billion voxel Space Shuttle simulation" width="50%"></img></a><a href="https://youtu.be/3JNVBQyetMA"><img src="https://img.youtube.com/vi/3JNVBQyetMA/maxresdefault.jpg" alt="Star Wars X-wing simulation" width="50%"></img></a><br>

## Compute Features
<b>CFD model: Lattice-Boltzmann method</b>

  - Streaming (2/2): <p align="center"><i>f</i><sub>0</sub><sup>temp</sup>(<i>x</i>,<i>t</i>) = <i>f</i><sub>0</sub>(<i>x</i>, <i>t</i>)<br>
<i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>) = <i>f</i><sub>(<i>t</i>%2 ? <i>i</i> : (<i>i</i>%2 ? <i>i</i>+1 : <i>i</i>-1))</sub>(<i>i</i>%2 ? <i>x</i> : <i>x</i>-<i>e<sub>i</sub></i>, <i>t</i>) &nbsp; for &nbsp; <i>i</i> &isin; [1, <i>q</i>-1]</p>
  - Collision: <p align="center"><i>&rho;</i>(<i>x</i>,<i>t</i>) = (&Sigma;<sub><i>i</i></sub> <i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>)) + 1<br><br>
<i>u</i>(<i>x</i>,<i>t</i>) = <sup>1</sup>&#8725;<sub><i>&rho;</i>(<i>x</i>,<i>t</i>)</sub> &Sigma;<sub><i>i</i></sub> <i>c<sub>i</sub></i> <i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>)<br><br>
<i>f<sub>i</sub></i><sup>eq-shifted</sup>(<i>x</i>,<i>t</i>) = <i>w<sub>i</sub></i> <i>&rho;</i> · (<sup>(<i>u</i><sub>°</sub><i>c<sub>i</sub></i>)<sup>2</sup></sup>&#8725;<sub>(2<i>c</i><sup>4</sup>)</sub> - <sup>(<i>u</i><sub>°</sub><i>u</i>)</sup>&#8725;<sub>(2c<sup>2</sup>)</sub> + <sup>(<i>u</i><sub>°</sub><i>c<sub>i</sub></i>)</sup>&#8725;<sub><i>c</i><sup>2</sup></sub>) + <i>w<sub>i</sub></i> (<i>&rho;</i>-1)<br><br>
<i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>, <i>t</i>+&Delta;<i>t</i>) = <i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>) + <i>&Omega;<sub>i</sub></i>(<i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>), <i>f<sub>i</sub></i><sup>eq-shifted</sup>(<i>x</i>,<i>t</i>), <i>&tau;</i>)</p>
  - Streaming (1/2): <p align="center"><i>f</i><sub>0</sub>(<i>x</i>, <i>t</i>+&Delta;<i>t</i>) = <i>f</i><sub>0</sub><sup>temp</sup>(<i>x</i>, <i>t</i>+&Delta;<i>t</i>)<br>
<i>f</i><sub>(<i>t</i>%2 ? (<i>i</i>%2 ? <i>i</i>+1 : <i>i</i>-1) : <i>i</i>)</sub>(<i>i</i>%2 ? <i>x</i>+<i>e<sub>i</sub></i> : <i>x</i>, <i>t</i>+&Delta;<i>t</i>) = <i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>, <i>t</i>+&Delta;<i>t</i>) &nbsp; for &nbsp; <i>i</i> &isin; [1, <i>q</i>-1]</p>

---

### Velocity Sets:
  - D2Q9
  - D3Q15
  - D3Q19 (default)
  - D3Q27

### Collision Operators:
  - Single-relaxation-time (SRT/BGK) (default)
  - Two-relaxation-time (TRT)

### Only 8 flag bits per lattice point (can be used independently / at the same time):
  - `TYPE_S` (stationary or moving) solid boundaries
  - `TYPE_E` equilibrium boundaries (inflow/outflow)
  - `TYPE_T` temperature boundaries
  - `TYPE_F` free surface (fluid)
  - `TYPE_I` free surface (interface)
  - `TYPE_G` free surface (gas)
  - `TYPE_X` remaining for custom use or further extensions
  - `TYPE_Y` remaining for custom use or further extensions



## Setup
Clone this repository:
```
git clone https://github.com/morkev/FluidX3D-easy-run.git
```

Go to the `stl folder`, and add the `.stl` file (or compatible format) there.

Declare whatever file you what to simulate in `void main_setup()` located in `setup.cpp`:
```cpp
void main_setup() {
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
	lbm.voxelize_stl(get_exe_path()+"../stl/F117.stl", center, rotation, size); // replace F117.stl with your file
	const uint N=lbm.get_N(), Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(uint n=0u, x=0u, y=0u, z=0u; n<N; n++, lbm.coordinates(n, x, y, z)) {
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	lbm.run();
} /**/
```


## Compatibility
- Works in Windows, Linux with C++17
- Supports importing and voxelizing triangle meshes from binary `.stl` files
- Supports exporting volumetric data as binary `.vtk` files
- Supports exporting rendered frames as `.png`/`.qoi`/`.bmp` files; time-consuming image encoding is handled in parallel on the CPU while the simulation on GPU can continue without delay


## References
- Lehmann, M.: [Esoteric Pull and Esoteric Push: Two Simple In-Place Streaming Schemes for the Lattice Boltzmann Method on GPUs](https://doi.org/10.3390/computation10060092). Computation, 10, 92, (2022)
- Lehmann, M., Krause, M., Amati, G., Sega, M., Harting, J. and Gekle, S.: [Accuracy and performance of the lattice Boltzmann method with 64-bit, 32-bit, and customized 16-bit number formats](https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats). Phys. Rev. E 106, 015308, (2022)
- Lehmann, M.: [Combined scientific CFD simulation and interactive raytracing with OpenCL](https://www.researchgate.net/publication/360501260_Combined_scientific_CFD_simulation_and_interactive_raytracing_with_OpenCL). IWOCL'22: International Workshop on OpenCL, 3, 1-2, (2022)
- Lehmann, M., Oehlschlägel, L.M., Häusl, F., Held, A. and Gekle, S.: [Ejection of marine microplastics by raindrops: a computational and experimental study](https://doi.org/10.1186/s43591-021-00018-8). Micropl.&Nanopl. 1, 18, (2021)
- Lehmann, M.: [High Performance Free Surface LBM on GPUs](https://doi.org/10.15495/EPub_UBT_00005400). Master's thesis, (2019)
- Lehmann, M. and Gekle, S.: [Analytic Solution to the Piecewise Linear Interface Construction Problem and Its Application in Curvature Calculation for Volume-of-Fluid Simulation Codes](https://doi.org/10.3390/computation10020021). Computation, 10, 21, (2022)


## Contact
- This is a fork of FluidX3D, developed by Dr. Moritz Lehmann [moritz.lehmann@uni-bayreuth.de](mailto:moritz.lehmann@uni-bayreuth.de?subject=FluidX3D).
- This is a ready to run version for those that don't have a heavy coding background, modified by Kevin Mora [kmoragon@asu.edu](mailto:kmoragon@asu.edu?subject=FluidX3D).
