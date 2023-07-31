# FluidX3D Documentation - How to get started?



## 1. Download
[Download](https://github.com/ProjectPhysX/FluidX3D/archive/refs/heads/master.zip) and unzip the source code, or clone with `git clone https://github.com/ProjectPhysX/FluidX3D.git`.

<br>

## 2. Compiling the Source Code
- There is no "installation" of the FluidX3D software. Instead, you have to compile the source code yourself.
- I have made this as easy as possible and this documentation will guide you through it. Nontheless, some basic programming experience with C++ would be good, as all the setup scripts are written in C++.
- First, compile the code as-is; this is the standard FP32 benchmark test case. By default, the fastest installed GPU will be selected automatically. Compile time is about 10 seconds.

### Windows
- Download and install [Visual Studio Community](https://visualstudio.microsoft.com/de/vs/community/). In Visual Studio Installer, add:
  - Desktop development with C++
  - MSVC v142
  - Windows 10 SDK
- Open [`FluidX3D.sln`](FluidX3D.sln) in [Visual Studio Community](https://visualstudio.microsoft.com/de/vs/community/).
- Compile and run by clicking the <kbd>► Local Windows Debugger</kbd> button.
- To select a specific GPU, open Windows CMD in the `FluidX3D` folder (type `cmd` in File Explorer in the directory field and press <kbd>Enter</kbd>), then run `bin\FluidX3D.exe 0` to select device `0`. You can also select multiple GPUs with `bin\FluidX3D.exe 0 1 3 6` if the setup is [configured as multi-GPU](#the-lbm-class).

### Linux
- Compile and run with `chmod +x make.sh` and `./make.sh`.
- Compiling requires `C++17`, which is supported since `g++` version `8`. Check with `g++ --version`.
- If you use [`INTERACTIVE_GRAPHICS`](src/defines.hpp), change to the "[compile on Linux with X11](make.sh#L6)" command in [`make.sh`](make.sh#L6).
- To select a specific GPU, enter `./make.sh 0` to compile+run, or `bin/FluidX3D 0` to run on device `0`. You can also select multiple GPUs with `bin/FluidX3D 0 1 3 6` if the setup is [configured as multi-GPU](#the-lbm-class).

### macOS
- Select the "[compile on macOS](make.sh#L9)" command in [`make.sh`](make.sh#L9).
- Compile and run with `chmod +x make.sh` and `./make.sh`.

### Android
- Select the "[compile on Android](make.sh#L10)" command in [`make.sh`](make.sh#L10).
- Compile and run with `chmod +x make.sh` and `./make.sh`.

<br>

## 3. Go through Sample Setups
- Now open [`src/setup.cpp`](src/setup.cpp). In here are all the sample setups, each one being a `void main_setup() {...}` function block written in C++. Uncomment one of them, maybe start top-to-bottom.
- In the line where the `main_setup()` function starts, it says "required extensions in defines.hpp:", followed by a list of extensions in capital letters. Head over to [`src/defines.hpp`](src/defines.hpp) and comment `//#define BENCHMARK` with a `//`. Then, uncomment all of the extensions required for the setup by removing the `//` in front of the corresponding line.
- Finally, [compile](#2-compiling-the-source-code) and run the setup with the <kbd>► Local Windows Debugger</kbd> button (Windows) or `./make.sh` (Linux/macOS).
- Once the interactive graphics window opens, press key <kbd>P</kbd> to start/pause the simulation, and press <kbd>H</kbd> to show the help menu for keyboard controls and visualization settings.
- Go through some of the sample setups this way, get familiar with their code structure and test the graphics mode.

<br>

## 4. Keyboard/Mouse Controls for [`INTERACTIVE_GRAPHICS`](src/defines.hpp)
- <kbd>P</kbd>: start/pause the simulation
- <kbd>H</kbd>: show/hide help menu for keyboard controls and visualization settings
- <kbd>1</kbd>: flag wireframe / solid surface (and force vectors on solid cells or surface pressure if the extension is used)
- <kbd>2</kbd>: velocity field
- <kbd>3</kbd>: streamlines
- <kbd>4</kbd>: vorticity (velocity-colored Q-criterion isosurface)
- <kbd>5</kbd>: rasterized free surface
- <kbd>6</kbd>: raytraced free surface
- <kbd>7</kbd>: particles
- <kbd>T</kbd>: toggle slice visualization mode
- <kbd>Q</kbd>/<kbd>E</kbd>: move slice in slice visualization mode
- <kbd>Mouse</kbd> or <kbd>I</kbd>/<kbd>J</kbd>/<kbd>K</kbd>/<kbd>L</kbd>: rotate camera
- <kbd>Scrollwheel</kbd> or <kbd>+</kbd>/<kbd>-</kbd>: zoom (centered camera mode) or camera movement speed (free camera mode)
- <kbd>Mouseclick</kbd> or <kbd>U</kbd>: toggle rotation with <kbd>Mouse</kbd> and angle snap rotation with <kbd>I</kbd>/<kbd>J</kbd>/<kbd>K</kbd>/<kbd>L</kbd>
- <kbd>Y</kbd>/<kbd>X</kbd>: adjust camera field of view
- <kbd>G</kbd>: print current camera position/rotation in console as copy/paste command
- <kbd>R</kbd>: toggle camera autorotation
- <kbd>F</kbd>: toggle centered/free camera mode
- <kbd>W</kbd>/<kbd>A</kbd>/<kbd>S</kbd>/<kbd>D</kbd>/<kbd>Space</kbd>/<kbd>C</kbd>: move free camera
- <kbd>V</kbd>: toggle stereoscopic rendering for VR
- <kbd>B</kbd>: toggle VR-goggles/3D-TV mode for stereoscopic rendering
- <kbd>N</kbd>/<kbd>M</kbd>: adjust eye distance for stereoscopic rendering
- <kbd>Esc</kbd>/<kbd>Alt</kbd>+<kbd>F4</kbd>: quit

<br>

## 5. Writing your own Setups

### The LBM Class
- For initializing the simulation box, use call
  ```c
  LBM lbm(Nx, Ny, Nz, nu, ...);
  ```
  constructor. `Nx`/`Ny`/`Nz` is the grid resolution and `nu` is the kinematic shear viscosity in [LBM units](#unit-conversion).
- To use multiple GPUs, use
  ```c
  LBM lbm(Nx, Ny, Nz, Dx, Dy, Dz, nu, ...);
  ```
  with `Dx`/`Dy`/`Dz` indicating how many domains (GPUs) there are in each spatial direction. The product `Dx`×`Dy`×`Dz` is the total number of domains (GPUs).
- As long as the `lbm` object is in scope, you can access the memory. As soon as it goes out of scope, all memory associated with the current simulation is freed again.
- The grid resolution `Nx`/`Ny`/`Nz` ultimately determines the VRAM occupation. Quite often it's not obvious at which resolution you'll overshoot the VRAM capacity of the GPU(s). To aid with this, there is the function:
  ```c
  const uint3 lbm_N = resolution(float3(1.0f, 2.0f, 0.5f), 2000u);
  ```
  This takes as inputs the desired aspect ratio of the simulation box and the VRAM occupation in MB, and returns the grid resolution as a `uint3` with `.x`/`.y`/`.z` components. You can also directly feed the `uint3` into the LBM constructor as resolution: `LBM lbm(lbm_N, nu, ...);`

### Unit Conversion
- The LBM simulation uses a different unit system from SI units, where density `rho=1` and velocity `u≈0.001-0.1`, because floating-point arithmetic is most accurate close to `1`.
- To ease unit conversion from SI to LBM units and back, there is the [`units.hpp`](src/units.hpp) struct. By calling
  ```c
  units.set_m_kg_s(lbm_length, lbm_velocity, lbm_density=1, si_length, si_velocity, si_density);
  ```
  the base unit conversion factors [m], [kg], [s] are calculated and stored in the `units` struct. Thereafter, any of the conversion functions from [`src/units.hpp`](src/units.hpp) can be used to go from SI to LBM units and back, such as `lbm_nu = units.nu(si_nu)` to convert the kinematic viscosity from SI to LBM units.
- A good beginner example for this is the "[aerodynamics of a cow](src/setup.cpp)" setup.

### Initial and Boundary Conditions
- If not explicitly set, by default all cells have the default values `rho=1`, `u=0`, `flags=0`.
- The initial/boundary conditions of single grid cells are set in a parallelized loop that iterates over the entire grid:
  ```c
  const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
  	// ...
  });
  ```
  Within this loop, you can set the density, velocity and flags of each cell individually by assigning values to `lbm.rho[n]`, `lbm.u.x[n]`, `lbm.u.y[n]`, `lbm.u.z[n]` and `lbm.flags[n]`. The `n` here is the linearized 3D grid index, corresponding to an (`x`|`y`|`z`) position via the function `lbm.coordinates(n, x, y, z)`.
- For example, to set solid walls at the left and right sides of the simulation box, write
  ```c
  if(y==0u||y==Ny-1u) lbm.flags[n] = TYPE_S;
  ```
  within the loop.
- All box sides where no solid (`TYPE_S`) or other boundary type are set will remain periodic boundaries.
- Primitive geometry, such as spheres, ellipsoids, cubes, cuboids, cylinders, codes, pipes, triangles, inclined planes, or toruses can be set with the functions from [`shapes.hpp`](src/shapes.hpp). Example to insert a cylinder:
  ```c
  if(cylinder(x, y, z, lbm.center(), float3(axis_x, axis_z, axis_z), radius) lbm.flags[n] = TYPE_S;
  ```
- The non-moving no-slip mid-grid bounce-back boundaries (`TYPE_S`) are always available without further extensions. "No-slip bounce-back" refers to the property that the flow velocity directly at the boundary is 0 (no-slip condition). "Mid-grid" refers to the boundary being located exactly in the middle between the boundary cells and adjacent fluid cells.
- For inflow/outflow boundaries, you need to enable (uncomment) the [`EQUILIBRIUM_BOUNDARIES`](src/defines.hpp) extension. Then, for the specific inflow/outflow cells, set the flag `lbm.flags[n] = TYPE_E` and on the same cell specify either a density `lbm.rho[n]` unequal to `1` or a velocity `lbm.u.x[n]`/`lbm.u.y[n]`/`lbm.u.z[n]` unequal to `0`, or a combination of both. `TYPE_E` cells enforce the specified density/velocity value and absorb any incoming shockwaves.
- For moving solid boundaries, you need to enable (uncomment) the [`MOVING_BOUNDARIES`](src/defines.hpp) extension. Then, for the specific solid cells, set the flag `lbm.flags[n] = TYPE_S` and on the same cell specify a velocity `lbm.u.x[n]`/`lbm.u.y[n]`/`lbm.u.z[n]` unequal to `0`. `TYPE_S` cells reflect any incoming shockwaves.
- If strict mass conservation is required (for example flow through a linear pipe), use periodic boundaries (i.e. don't set any boundary type on the cells at these simulation box sides), and drive the flow with a volume force (equivalent to a pressure gradient). Therefore you need to enable (uncomment) the [`VOLUME_FORCE`](src/defines.hpp) extension, and in the [LBM constuctor](#the-lbm-class) set the force per volume (`fx`|`fy`|`fz`):
  ```c
  LBM lbm(Nx, Ny, Nz, nu, fx, fy, fz);
  ```
  These force per volume values should not exceed `0.001` in magnitude.

### Running the Simulation
- Call `lbm.run()` (without input parameter, it's infinite time steps) to initialize and execute the setup, or `lbm.run(time_steps)` to execute only a specific number of time steps.
- If you have a [more complicated simulation loop](#video-rendering) where you periodically compute time steps and render images for a video or export data, don't forget to place an `lbm.run(0u)` before that loop. This copies the initial/boundary conditions from CPU RAM to GPU VRAM and initializes the simulation on the GPU, without computing a time step. Without initialization, there is no data in VRAM yet for rendering.

### Loading .stl Files
- For more complex geometries, you can load `.stl` triangle meshes and voxelize them to the Cartesian simulation grid on the GPU(s).
- Create a `FluidX3D/stl/` folder next to the `FluidX3D/src/` folder and download the geometry from websites like [Thingiverse](https://www.thingiverse.com/), or create your own.
- Only binary `.stl` files are supported. For conversion from other formats or for splitting composite geometries like helicopter hull and rotors, I recommend [Microsoft 3D Builder](https://apps.microsoft.com/store/detail/3d-builder/9WZDNCRFJ3T6) on Windows or [Blender](https://www.blender.org/) on Windows/Linux.
- Load and voxelize simple `.stl` files directly with
  ```c
  lbm.voxelize_stl(get_exe_path()+"../stl/mesh.stl", center, rotation, size);
  ```
  This automatically repositions/rescales the mesh to the specified center. Use `lbm.center()` for the simulation box center, or add an offset with a `+float3(offset_x, offset_y, offset_z)`. You can generate and multiply together a rotation matrix like this (example: rotation around the z-axis by 180°, then around the x-axis by 90°):
  ```c
  float3x3 rotation = float3x3(float3(1, 0, 0), radians(90.0f))*float3x3(float3(0, 0, 1), radians(180.0f));
  ```
- To load composite geometries with several parts without automatic mesh repositioning/rescaling, use
  ```c
  Mesh* mesh_1 = read_stl(const string& path, const float scale=1.0f, const float3x3& rotation=float3x3(1.0f), const float3& offset=float3(0.0f)); // load mesh without automatic repositioning/rescaling
  Mesh* mesh_2 = read_stl(const string& path, const float scale=1.0f, const float3x3& rotation=float3x3(1.0f), const float3& offset=float3(0.0f));
  mesh_1->scale(const float scale); // manually scale meshes
  mesh_2->scale(const float scale);
  mesh_1->translate(const float3& translation); // manually reposition meshes
  mesh_2->translate(const float3& translation);
  lbm.voxelize_mesh_on_device(mesh_1); // voxelize meshes on GPU
  lbm.voxelize_mesh_on_device(mesh_2);
  
  ```
  to load the meshes from the `.stl` files, manually scale/reposition all parts of the mesh the same time, and finally voxelize them on the GPU.
- To aid with repositioning the mesh, there is `lbm.center()` for the center of the simulation box, as well as the min/max bounding-box coordinates of the mesh `mesh->pmin`/`mesh->pmax`, each a `float3` with (`x`|`y`|`z`) components.
- Rotating geometries have to be periodically revoxelized, about every 1-10 LBM time steps. In the main simulation loop in the [`main_setup()`](src/setup.cpp) function, first rotate the triangle mesh, then revoxelize on GPU, then compute a few LBM time steps:
  ```c
  const uint lbm_T = 100000u; // number of LBM time steps to simulate
  const uint lbm_dt = 4u; // number of LBM time steps between each mesh revoxelization
  lbm.run(0u); // initialize simulation
  while(lbm.get_t()<lbm_T) { // main simulation loop
  	mesh->rotate(float3x3(float3(0, 0, 1), lbm_omega*(float)lbm_dt)); // rotate the triangle mesh
  	lbm.voxelize_mesh_on_device(mesh, TYPE_S, center, float3(0.0f), float3(0.0f, 0.0f, lbm_omega)); // revoxelize the rotated triangle mesh, provide the instantaneous angular velocity vector for moving boundaries
  	lbm.run(lbm_dt); // run lbm_dt LBM time steps
  }
  ```
  Here `lbm_omega` is the angular velocity in radians per time step, `lbm_dt` is the number of simulated time steps between revoxelizations, and `float3(0.0f, 0.0f, lbm_omega)` is the instantaneous angular velocity as a vector along the axis of rotation. The largest displacement of the outermost cells should not exceed `1` cell between revoxelizations; set `lbm_omega = lbm_u/lbm_radius` accordingly.
- Have a look at the "[Cessna 172](src/setup.cpp)" and "[Bell 222](src/setup.cpp)" setups for some examples.

### Video Rendering
- For video rendering, disable (comment out) [`INTERACTIVE_GRAPHICS`](src/defines.hpp) and [`INTERACTIVE_GRAPHICS_ASCII`](src/defines.hpp) and enable (uncomment) [`GRAPHICS`](src/defines.hpp) in [`src/defines.hpp`](src/defines.hpp).
- Set the video resolution as [`GRAPHICS_FRAME_WIDTH`](src/defines.hpp)/[`GRAPHICS_FRAME_HEIGHT`](src/defines.hpp) and the background color as [`GRAPHICS_BACKGROUND_COLOR`](src/defines.hpp). You can also adjust the other [`GRAPHICS_...`](src/defines.hpp) options there, such as semi-transparent rendering mode, or adjust the color scale for velocity with [`GRAPHICS_U_MAX`](src/defines.hpp).
- A basic loop for rendering video in the [`main_setup()`](src/setup.cpp) function looks like this:
  ```c
  lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_Q_CRITERION; // set visualization modes, see all available visualization mode macros (VIZ_...) in defines.hpp
  const uint lbm_T = 10000u; // number of LBM time steps to simulate
  lbm.run(0u); // initialize simulation
  while(lbm.get_t()<lbm_T) { // main simulation loop
  	if(lbm.graphics.next_frame(lbm_T, 25.0f)) { // render enough frames for 25 seconds of 60fps video
  		lbm.graphics.set_camera_free(float3(2.5f*(float)Nx, 0.0f*(float)Ny, 0.0f*(float)Nz), 0.0f, 0.0f, 50.0f); // set camera to position 1
  		lbm.graphics.write_frame(get_exe_path()+"export/camera_angle_1/"); // export image from camera position 1
  		lbm.graphics.set_camera_centered(-40.0f, 20.0f, 78.0f, 1.25f); // set camera to position 2
  		lbm.graphics.write_frame(get_exe_path()+"export/camera_angle_2/"); // export image from camera position 2
  	}
  	lbm.run(1u); // run 1 LBM time step
  }
  ```
- To find suitable camera placement, run the simulation at low resolution in [`INTERACTIVE_GRAPHICS`](src/defines.hpp) mode, rotate/move the camera to the desired position, click the <kbd>Mouse</kbd> to disable mouse rotation, and press <kbd>G</kbd> to print the current camera settings as a copy-paste command in the console. <kbd>Alt</kbd>+<kbd>Tab</kbd> to the console and copy the camera placement command by selecting it with the mouse and right-clicking, then paste it into the [`main_setup()`](src/setup.cpp) function.
- The visualization mode(s) can be specified as `lbm.graphics.visualization_modes` with the [`VIS_...`](src/defines.hpp) macros. You can also set the `lbm.graphics.slice_mode` (`0`=no slice, `1`=x, `2`=y, `3`=z, `4`´=xz, `5`=xyz, `6`=yz, `7`=xy) and reposition the slices with `lbm.graphics.slice_x`/`lbm.graphics.slice_y`/`lbm.graphics.slice_z`.
- Exported frames will automatically be assigned the current simulation time step in their name, in the format `bin/export/image-123456789.png`.
- To convert the rendered `.png` images to video, use [FFmpeg](https://ffmpeg.org/):
  ```bash
  ffmpeg -framerate 60 -pattern_type glob -i "./bin/export/*/image-*.png" -c:v libx264 -pix_fmt yuv420p -b:v 24M "video.mp4"
  ```

### Data Export
- At any point in time, you can export volumetric data as binary `.vtk` files with:
  ```c
  lbm.rho.write_device_to_vtk();
  lbm.u.write_device_to_vtk();
  lbm.flags.write_device_to_vtk();
  lbm.phi.write_device_to_vtk(); // only for SURFACE extension
  lbm.T.write_device_to_vtk(); // only for TEMPERATURE extension
  lbm.write_mesh_to_vtk(const Mesh* mesh); // for exporting triangle meshes
  ```
- These functions first pull the data from the GPU(s) into CPU RAM, and then write it to the hard drive.
- Exported files will automatically be assigned the current simulation time step in their name, in the format `bin/export/u-123456789.vtk`.
- Be aware that these volumetric files can be gigantic in file size, tens of GigaByte for a single file.
- You can view/evaluate the `.vtk` files for example in [ParaView](https://www.paraview.org/).
- It is recommended to use the C++ functionality in the [`main_setup()`](src/setup.cpp) function directly to extract the data of interest and selectively only write that to the hard drive. Therefore, call `lbm.u.read_from_device()` to copy the data from the GPU(s) to CPU RAM, and then you can access it directly, for example
  ```c
  const float lbm_velocity_x = lbm.u.x[lbm.index(x, y, z)];
  ```
  to get the x-velocity at the position (`x`|`y`|`z`) in [LBM units](#unit-conversion).
- To [convert the velocity from LBM to SI units](#unit-conversion), use
  ```c
  const float si_velocity_x = units.si_u(lbm_velocity_x);
  ```
  after having done [unit conversion](#unit-conversion) with `units.set_m_kg_s(...)`.
- You can also export the `.stl` triangle meshes to binary `.vtk` files with:
  ```c
  lbm.write_mesh_to_vtk(const Mesh* mesh);
  ```

### Lift/Drag Forces
- Enable (uncomment) the [`FORCE_FIELD`](src/defines.hpp) extension. This extension allows computing boundary forces on every solid cell (`TYPE_S`) individually, as well as placing an individual volume force on every fluid cell (not used here).
- In the [`main_setup()`](src/setup.cpp) function's main simulation loop, alternatingly call:
  ```c
  lbm.run(lbm_dt); // run lbm_dt LBM time steps
  lbm.calculate_force_on_boundaries(); // compute boundary forces on GPU on all solid cells (TYPE_S)
  ```
  The latter computes the boundary forces on the GPU into the `lbm.F` field in VRAM.
- To copy `lbm.F` from GPU VRAM to CPU RAM, call:
  ```c
  lbm.F.read_from_device();
  ```
  You can then access the boundary forces at each individual cell with:
  ```c
  float lbm_force_x_n = lbm.F.x[lbm.index(x, y, z)];
  ```
- To sum over all the individual boundary cells that belong to the body, to get the total force on the body, first voxelize the body with
  ```c
  lbm.voxelize_mesh_on_device(mesh, TYPE_S|TYPE_X);
  ```
  with the additional `TYPE_X` flagging, and then call
  ```c
  const float3 lbm_force = lbm.calculate_force_on_object(TYPE_S|TYPE_X);
  ```
  to sum over all cells marked `TYPE_S|TYPE_X` that belong to the body. You can also use `TYPE_Y` for this.
- Finally, [convert from LBM to SI units](#unit-conversion) with
  ```c
  const float si_force_x = units.si_F(lbm_force.x);
  ```
  after having done [unit conversion](#unit-conversion) with `units.set_m_kg_s(...)`.
- See the "Ahmed body" setup for an example. Note that in the highly turbulent regime, computed body forces are too large by up to a factor 2, because even large resolution is not enough to fully capture the turbulent boundary layer. A wall function is needed, I'll scan literature on it.

<br>

## 6. Further LBM Extensions
By now you're already familiar with the [additional boundary types](#initial-and-boundary-conditions) through extensions [`VOLUME_FORCE`](src/defines.hpp), [`FORCE_FIELD`](src/defines.hpp), [`EQUILIBRIUM_BOUNDARIES`](src/defines.hpp), and [`MOVING_BOUNDARIES`](src/defines.hpp). The remaining available model extensions are briefly outlined here:

### [`SURFACE`](src/defines.hpp) Extension
- To simulate free water surfaces, enable (uncomment) the [`SURFACE`](src/defines.hpp) extension.
- All cells then get 3 additional flags: `TYPE_F` (fluid), `TYPE_I` (interface), and `TYPE_G` (gas). Fluid cells are computed with regular LBM. Interface cells account for the extra surface tension forces, if the surface tension coefficient `sigma` is set greater than `0` in the [LBM constructor](#the-lbm-class); the interface is always 1 cell layer thick. Gas cells are not simulated at all and are essentially treated as vacuum.
- If not set otherwise in the [initial conditions](#initial-and-boundary-conditions), all cells are initialized as `TYPE_G` by default. As initial conditions, set all cells that should be fluid to
  ```c
  lbm.flags[n] = TYPE_F;
  ```
  The interface layer will be automatically initialized during initialization with `lbm.run(0u)`.
- Addidionally to the 3 flags, each cell also gets assigned a fill level `lbm.phi[n]`: `1` for fluid cells (`TYPE_F`), `]0,1[` for interface cells (`TYPE_I`), and `0` for gas cells (`TYPE_G`). You can set this fill level at initialization, additionally to the cell flag. Do not forget to set the cell flag. If `lbm.phi[n]` is not set manually, it will automatically be initialized such that all fluid cells get `phi=1`, all interface cells get `phi=0.5`, and all gas clls get `phi=0` assigned.
- For a simple example, see the "[dam break](src/setup.cpp)" setup. A more advanced sample setup for free surfaces is the "[raindrop impact](src/setup.cpp)".

### [`TEMPERATURE`](src/defines.hpp) Extension
- With the [`TEMPERATURE`](src/defines.hpp) extension, FluidX3D can model thermal convection flows. This extension automatically also enables the [`VOLUME_FORCE`](src/defines.hpp) extension.
- In the [LBM constructor](#the-lbm-class), you then need to set the volume force (`fx`|`fy`|`fz`), the thermal diffusion coefficient `alpha`, and the thermal expansion coefficient `beta`, all [in LBM units](#unit-conversion):
  ```c
  LBM lbm(Nx, Ny, Nz, nu, fx, fy, fz, 0.0f, alpha, beta); // the "0.0f" is for the surface tension coefficient sigma which is not used here and has to remain 0
  ```
- With the extension, each grid cell gets an additional temperature `lbm.T[n]` ([in LBM units](#unit-conversion)) assigned. The default temperature in LBM units is `1`.
- To set temperature boundary conditions, use the flag `TYPE_T` and for the same cells assign a temperature unequal to `1`:
  ```c
  lbm.flags[n] = TYPE_T; // make the cell n a temperature boundary
  lbm.T[n] = 1.2f; // set this temperature boundary hotter than average
  ```
- See the "[Rayleigh-Benard convection](src/setup.cpp)" and "[thermal convection](src/setup.cpp)" setups for two examples.

### [`SUBGRID`](src/defines.hpp) Extension
- Fluid flow is characterized by the Reynolds number<p align="center"><i>Re</i> = <sup><i>x</i>·<i>u</i></sup>&#8725;<sub><i>nu</i></sub></p>with a characteristic length scale `x`, a characteristic velocity `u` and the kinematic shear viscosity `nu`. Larger length scale, larger velocity or smaller viscosity all mean larger Reynolds number.
- The Reynolds number is a unit-less number. A low value <i>Re</i> < 2300 means laminar flow, a high value <i>Re</i> > 2900 means turbulent flow. In between is a transitional regime.
- For very large Reynolds number <i>Re</i> > 100000, the LBM solver becomes [unstable](#7-suitable-parameters-and-simulation-instability), as tiny, very fast rotating vortices can be present in the flow field, and too fast velocity and shear rate makes the simulation blow up.
- To tackle this problem, there is subgrid models that model vortices smaller than single grid cells. This works by increasing the effective vscosity where the shear rate is large and lots of small eddies are assumed to be present. Coincidentally, locations of high shear rate and low viscosity cause instability, so increasing effective viscosity there keeps the simulation stable.
- The subgrid model in FLuidX3D is the Smagorinsky-Lilly model. You can enable it with the [`SUBGRID`](src/defines.hpp) extension.
- There is no additional performance cost for this extension.

### [`PARTICLES`](src/defines.hpp) Extension
- By default, the LBM is a grid-based simulation, so there are no particles.
- But the [`PARTICLES`](src/defines.hpp) extension allows to add particles to the simulation, either as passive tracers or as 2-way-coupled particles that can do floating/sedimentation.
- For passive tracers, only enable the [`PARTICLES`](src/defines.hpp) extension, and in the [LBM constructor](#the-lbm-class) simply add the particle count:
  ```c
  LBM lbm(Nx, Ny, Nz, nu, 50000u); // this will create 50000 particles
  ```
- Then, in [initialization](#initial-and-boundary-conditions), make a loop over all particles (outside of the initialization loop that iterates over all grid cells):
  ```c
  for(ulong n=0ull; n<lbm.particles->length(); n++) {
  	lbm.particles->x[n] = random_symmetric(0.5f*lbm.size().x); // this will palce the particles randomly anywhere in the simulation box
  	lbm.particles->y[n] = random_symmetric(0.5f*lbm.size().y);
  	lbm.particles->z[n] = random_symmetric(0.5f*lbm.size().z);
  }
  ```
- Note that the position (`0`|`0`|`0`) for particles corresponds to the simulation box center.
- For 2-way-coupled particles, additionally enable the [`VOLUME_FORCE`](src/defines.hpp) and [`FORCE_FIELD`](src/defines.hpp) extensions, and in the [LBM constructor](#the-lbm-class) add the particle density ([in LBM units](#unit-conversion)) unequal to `1`:
  ```c
  LBM lbm(Nx, Ny, Nz, nu, 50000u, 1.2f); // this will create 50000 particles that are more dense than the fluid and will sink to the bottom
  ```

<br>

## 7. Suitable Parameters and Simulation Instability
- Sometimes in the velocity field or streamlines visualization, you will see fuzzyness, or something that looks like a rapidly growing white crystal, blowing up from a certain point and filling the entire simulation box. This is instability, i.e. when velocities turn `NaN` or `Inf`.
- Often times, the cause of instability is an unfortunate choice of unsuitable parameters:
  - too high/low density `rho` (ideally should be very close to `1` at all times)
  - too high velocity `u` (must never exceed `0.57` anywhere in the box, ideally should be somewhere around `0.1`, but can be as small as `0.001`)
  - too low kinematic shear viscosity `nu` (ideally close to `1/6`, becomes unstable when it's very very close to `0` (then enable the [`SUBGRID`](src/defines.hpp) extension), and should not exceed `3`)
  - too high force per volume (`fx`|`fy`|`fz`) (should not exceed `0.001` in magnitude)
  - too high surface tension coefficient `sigma` (should not exceed `0.1`)
- The best parametrization for LBM simulations is an art in itself and needs some practice.