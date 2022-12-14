# FluidX3D

The fastest and most memory efficient lattice Boltzmann CFD software, running on any GPU via [OpenCL](https://github.com/ProjectPhysX/OpenCL-Wrapper "OpenCL-Wrapper").

<a href="https://youtu.be/5AzxwQpng0M"><img src="https://img.youtube.com/vi/5AzxwQpng0M/maxresdefault.jpg" alt="10 billion voxel Space Shuttle simulation" width="50%"></img></a><a href="https://youtu.be/VadLwt9OqMo"><img src="https://img.youtube.com/vi/VadLwt9OqMo/maxresdefault.jpg" alt="1 billion voxel raindrop simulation" width="50%"></img></a><br>
<a href="https://youtu.be/NQPgumd3Ei8"><img src="https://img.youtube.com/vi/NQPgumd3Ei8/maxresdefault.jpg" alt="Hydraulic jump simulation" width="50%"></img></a><a href="https://youtu.be/3JNVBQyetMA"><img src="https://img.youtube.com/vi/3JNVBQyetMA/maxresdefault.jpg" alt="Star Wars X-wing simulation" width="50%"></img></a>



## Compute Features

- CFD model: lattice Boltzmann method (LBM)
  - streaming (part 2/2): <p align="center"><i>f</i><sub>0</sub><sup>temp</sup>(<i>x</i>,<i>t</i>) = <i>f</i><sub>0</sub>(<i>x</i>, <i>t</i>)<br>
<i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>) = <i>f</i><sub>(<i>t</i>%2 ? <i>i</i> : (<i>i</i>%2 ? <i>i</i>+1 : <i>i</i>-1))</sub>(<i>i</i>%2 ? <i>x</i> : <i>x</i>-<i>e<sub>i</sub></i>, <i>t</i>) &nbsp; for &nbsp; <i>i</i> &isin; [1, <i>q</i>-1]</p>
  - collision: <p align="center"><i>&rho;</i>(<i>x</i>,<i>t</i>) = (&Sigma;<sub><i>i</i></sub> <i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>)) + 1<br><br>
<i>u</i>(<i>x</i>,<i>t</i>) = <sup>1</sup>&#8725;<sub><i>&rho;</i>(<i>x</i>,<i>t</i>)</sub> &Sigma;<sub><i>i</i></sub> <i>c<sub>i</sub></i> <i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>)<br><br>
<i>f<sub>i</sub></i><sup>eq-shifted</sup>(<i>x</i>,<i>t</i>) = <i>w<sub>i</sub></i> <i>&rho;</i> · (<sup>(<i>u</i><sub>°</sub><i>c<sub>i</sub></i>)<sup>2</sup></sup>&#8725;<sub>(2<i>c</i><sup>4</sup>)</sub> - <sup>(<i>u</i><sub>°</sub><i>u</i>)</sup>&#8725;<sub>(2c<sup>2</sup>)</sub> + <sup>(<i>u</i><sub>°</sub><i>c<sub>i</sub></i>)</sup>&#8725;<sub><i>c</i><sup>2</sup></sub>) + <i>w<sub>i</sub></i> (<i>&rho;</i>-1)<br><br>
<i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>, <i>t</i>+&Delta;<i>t</i>) = <i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>) + <i>&Omega;<sub>i</sub></i>(<i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>), <i>f<sub>i</sub></i><sup>eq-shifted</sup>(<i>x</i>,<i>t</i>), <i>&tau;</i>)</p>
  - streaming (part 1/2): <p align="center"><i>f</i><sub>0</sub>(<i>x</i>, <i>t</i>+&Delta;<i>t</i>) = <i>f</i><sub>0</sub><sup>temp</sup>(<i>x</i>, <i>t</i>+&Delta;<i>t</i>)<br>
<i>f</i><sub>(<i>t</i>%2 ? (<i>i</i>%2 ? <i>i</i>+1 : <i>i</i>-1) : <i>i</i>)</sub>(<i>i</i>%2 ? <i>x</i>+<i>e<sub>i</sub></i> : <i>x</i>, <i>t</i>+&Delta;<i>t</i>) = <i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>, <i>t</i>+&Delta;<i>t</i>) &nbsp; for &nbsp; <i>i</i> &isin; [1, <i>q</i>-1]</p>

<!-- markdown equations don't render properly in mobile browser

  - streaming (part 2/2):
$$j=0\\ \textrm{for}\\ i=0$$
$$j=t\\%2\\ ?\\ i\\ :\\ (i\\%2\\ ?\\ i+1\\ :\\ i-1)\\ \textrm{for}\\ i\in[1,q-1]$$
$$f_i^\textrm{temp}(\vec{x},t)=f_j(i\\%2\\ ?\\ \vec{x}\\ :\\ \vec{x}-\vec{e}_i,\\ t)$$
  - collision:
$$\rho(\vec{x},t)=\left(\sum_i f_i^\textrm{temp}(\vec{x},t)\right)+1$$
$$\vec{u}(\vec{x},t)=\frac{1}{\rho(\vec{x},t)}\sum_i\vec{c}_i f_i^\textrm{temp}(\vec{x},t)$$
$$f_i^\textrm{eq-shifted}(\vec{x},t)=w_i \rho \cdot\left(\frac{(\vec{u} _{^{^\circ}}\vec{c}_i)^2}{2 c^4}-\frac{\vec{u} _{^{^\circ}}\vec{u}}{2 c^2}+\frac{\vec{u} _{^{^\circ}}\vec{c}_i}{c^2}\right)+w_i (\rho-1)$$
$$f_i^\textrm{temp}(\vec{x},\\ t+\Delta t)=f_i^\textrm{temp}(\vec{x},t)+\Omega_i(f_i^\textrm{temp}(\vec{x},t),\\ f_i^\textrm{eq-shifted}(\vec{x},t),\\ \tau)$$
  - streaming (part 1/2):
$$j=0\\ \textrm{for}\\ i=0$$
$$j=t\\%2\\ ?\\ (i\\%2\\ ?\\ i+1\\ :\\ i-1)\\ :\\ i\\ \textrm{for}\\ i\in[1,q-1]$$
$$f_j(i\\%2\\ ?\\ \vec{x}+\vec{e}_i\\ :\\ \vec{x},\\ t+\Delta t)=f_i^\textrm{temp}(\vec{x},\\ t+\Delta t)$$

 -->

- peak performance on most GPUs (datacenter/gaming/professional/laptop), validated with roofline model
- up to 4.29 billion (2³²) lattice points, or 1624³ resolution, on a single GPU (if it has 225 GB memory)
- optimized to minimize memory demand to 55 Bytes/node (~⅙ (~⅓) of conventional FP64 (FP32) LBM solvers)
  - in-place streaming with [Esoteric-Pull](https://doi.org/10.3390/computation10060092): almost cuts memory demand in half and slightly increases performance due to implicit bounce-back boundaries; offers optimal memory access patterns for single-node in-place streaming
  - [decoupled arithmetic precision (FP32) and memory precision (FP32 or FP16S or FP16C)](https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats): all arithmetic is done in FP32 for compatibility on all hardware, but LBM density distribution functions in memory can be compressed to FP16S or FP16C: almost cuts memory demand in half again and almost doubles performance, without impacting overall accuracy for most setups
- [DDF-shifting](https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats) and other algebraic optimization to minimize round-off error
- updating density and velocity in the `stream_collide()` kernel is optional (higher performance if disabled)
- velocity sets:
  - D2Q9
  - D3Q15
  - D3Q19 (default)
  - D3Q27
- collision operators
  - single-relaxation-time (SRT/BGK) (default)
  - two-relaxation-time (TRT)
- only 8 flag bits per lattice point (can be used independently / at the same time):
  - `TYPE_S` (stationary or moving) solid boundaries
  - `TYPE_E` equilibrium boundaries (inflow/outflow)
  - `TYPE_T` temperature boundaries
  - `TYPE_F` free surface (fluid)
  - `TYPE_I` free surface (interface)
  - `TYPE_G` free surface (gas)
  - `TYPE_X` remaining for custom use or further extensions
  - `TYPE_Y` remaining for custom use or further extensions



## Optional Compute Extensions

- [boundary types](https://doi.org/10.15495/EPub_UBT_00005400)
  - stationary mid-grid bounce-back boundaries (stationary solid boundaries)
  - moving mid-grid bounce-back boundaries (moving solid boundaries)
  - equilibrium boundaries (non-reflective inflow/outflow)
  - temperature boundaries (fixed temperature)
- global force per volume (Guo forcing), can be modified on-the-fly
- local force per volume (force field)
  - optional computation of forces from the fluid on solid boundaries
- state-of-the-art [free surface LBM](https://doi.org/10.3390/computation10060092) (FSLBM) implementation:
  - [volume-of-fluid model](https://doi.org/10.15495/EPub_UBT_00005400)
  - [fully analytic PLIC](https://doi.org/10.3390/computation10020021) for efficient curvature calculation
  - improved mass conservation
  - ultra efficient implementation with only [4 kernels](https://doi.org/10.3390/computation10060092) additionally to `stream_collide()` kernel
- thermal LBM to simulate thermal convection
  - D3Q7 subgrid for thermal DDFs
  - in-place streaming with [Esoteric-Pull](https://doi.org/10.3390/computation10060092) for thermal DDFs
  - optional [FP16S or FP16C compression](https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats) for thermal DDFs with [DDF-shifting](https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats)
- Smagorinsky-Lilly subgrid turbulence LES model to keep simulations with very large Reynolds number stable <details><summary>equations</summary><p align="center"><i>&Pi;<sub>&alpha;&beta;</sub></i> = &Sigma;<sub><i>i</i></sub> <i>e<sub>i&alpha;</sub></i> <i>e<sub>i&beta;</sub></i> (<i>f<sub>i</sub></i> - <i>f<sub>i</sub></i><sup>eq-shifted</sup>)<br><br>
Q = &Sigma;<sub><i>&alpha;&beta;</i></sub> <i>&Pi;<sub>&alpha;&beta;</sub></i><sup>2</sup><br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;______________________<br>
&tau; = &frac12; (&tau;<sub>0</sub> + &radic; &tau;<sub>0</sub><sup>2</sup> + <sup>(16&radic;2)</sup>&#8725;<sub>(<i>3&pi;</i><sup>2</sup>)</sub> <sup>&radic;Q</sup>&#8725;<sub><i>&rho;</i></sub> )
</p></details>




## Graphics Features

- on Windows and Linux: real time [interactive rasterization and raytracing graphics](https://www.researchgate.net/publication/360501260_Combined_scientific_CFD_simulation_and_interactive_raytracing_with_OpenCL)
- on Windows and Linux (even in WSL and/or in a remote console through SSH): real time interactive ASCII console graphics
- with interactive graphics mode disabled, image resolution can be as large as VRAM allows for (132 Megapixel (16K) and above)
- (interacitive) visualization modes:
  - flags (and force vectors on solid boundary nodes if the extension is used)
  - velocity field
  - streamlines
  - velocity-colored Q-criterion isosurface
  - rasterized free surface with [marching-cubes](http://paulbourke.net/geometry/polygonise/)
  - [raytraced free surface](https://www.researchgate.net/publication/360501260_Combined_scientific_CFD_simulation_and_interactive_raytracing_with_OpenCL) with fast ray-grid traversal and marching-cubes, either 1-4 rays/pixel or 1-10 rays/pixel



## How to get started?

1. Check the settings and extensions in [`src/defines.hpp`](src/defines.hpp) by uncommenting corresponding lines.
2. Write a C++ setup skript as `main_setup()` function in [`src/setup.cpp`](src/setup.cpp) (get inspiration from existing setups):
   - For unit conversion, use the `units` struct.
   - For initializing the box, use call `LBM lbm(...);` constructor.
   - Set the initial condition in a loop that iterates over the entire lattice by writing to `lbm.rho[n]`/`lbm.u.x[n]`/`lbm.u.y[n]`/`lbm.u.z[n]`/`lbm.flags[n]`.
   - Call `lbm.run()` to initialize and execute the setup (infinite time steps) or `lbm.run(time_steps)` to execute only a specific number of time steps.
   - As long as the `lbm` object is in scope, you can access the memory. As soon as it goes out of scope, all memory associated to the current simulation is freed again.
3. When done with the setup, hit compile+run on Windows in Visual Studio Community or execute `./make.sh` on Linux; this will automatically select the fastest installed GPU. Alternatively, you can add the device ID as the first command-line argument, for example `./make.sh 2` to compile+run on device 2, or `bin/FluidX3D.exe 2` to run the executable on device 2. Compile time for the entire code is about 10 seconds. If you use `INTERACTIVE_GRAPHICS` on Linux, change to the "compile on Linux with X11" command in `make.sh`.
4. With `INTERACTIVE_GRAPHICS`/`INTERACTIVE_GRAPHICS_ASCII` enabled, press the <kbd>P</kbd> key to start/pause the simulation.

   Toggle rendering modes with the keyboard:
   - <kbd>1</kbd>: flags (and force vectors on solid boundary nodes if the extension is used)
   - <kbd>2</kbd>: velocity field
   - <kbd>3</kbd>: streamlines
   - <kbd>4</kbd>: vorticity / velocity-colored Q-criterion isosurface
   - <kbd>5</kbd>: rasterized free surface
   - <kbd>6</kbd>: raytraced free surface

   Camera movement:
   - <kbd>Mouse</kbd> or <kbd>I</kbd>/<kbd>J</kbd>/<kbd>K</kbd>/<kbd>L</kbd>: rotate camera
   - <kbd>Scrollwheel</kbd> or <kbd>,</kbd>/<kbd>.</kbd>: zoom (centered camera mode) or adjust camera movement speed (free camera mode)
   - <kbd>U</kbd>: toggle rotation with <kbd>Mouse</kbd> and angle snap rotation with <kbd>I</kbd>/<kbd>J</kbd>/<kbd>K</kbd>/<kbd>L</kbd>
   - <kbd>Y</kbd>/<kbd>X</kbd>: adjust camera field of view
   - <kbd>R</kbd>: toggle camera autorotation
   - <kbd>F</kbd>: toggle centered/free camera mode
   - <kbd>W</kbd>/<kbd>A</kbd>/<kbd>S</kbd>/<kbd>D</kbd>/<kbd>Space</kbd>/<kbd>C</kbd>: move free camera
   - <kbd>V</kbd>: toggle stereoscopic rendering for VR
   - <kbd>B</kbd>: toggle VR-goggles/3D-TV mode for stereoscopic rendering
   - <kbd>N</kbd>/<kbd>M</kbd>: adjust eye distance for stereoscopic rendering
   - <kbd>H</kbd>: print current camera position/rotation in console as copy/paste command
   - <kbd>Esc</kbd>/<kbd>Alt</kbd>+<kbd>F4</kbd>: quit



## Compatibility

- works in Windows, Linux and Android with C++17
- runs on any hardware that supports OpenCL 1.2, from any vendor (Nvidia, AMD, Intel, ...):
  - world's fastest datacenter GPUs like H100, A100, MI250(X), MI210, MI100, V100(S), P100, ...
  - gaming GPUs (desktop or laptop)
  - "professional"/workstation GPUs
  - integrated GPUs
  - Xeon Phi
  - CPUs
  - even smartphone ARM GPUs
- supports importing and voxelizing triangle meshes from binary `.stl` files
- supports exporting volumetric data as binary `.vtk` files
- supports exporting rendered frames as `.png`/`.qoi`/`.bmp` files; time-consuming image encoding is handled in parallel on the CPU while the simulation on GPU can continue without delay



## Benchmarks

Here are [performance benchmarks](https://doi.org/10.3390/computation10060092) on various hardware in MLUPs/s, or how many million lattice points are updated per second. The settings used for the benchmark are D3Q19 SRT with no extensions enabled (only LBM with implicit mid-grid bounce-back boundaries) and the setup consists of an empty cubic box with sufficient size (typically 256³). Without extensions, a single lattice point requires:
- a memory capacity of 93 (FP32/FP32) or 55 (FP32/FP16) Bytes
- a memory bandwidth of 153 (FP32/FP32) or 77 (FP32/FP16) Bytes per time step
- 326 (FP32/FP32) or 364 (FP32/FP16S) or 1241 (FP32/FP16C) FLOPs per time step (FP32+INT32 operations counted combined)

In consequence, the arithmetic intensity of this implementation is 2.13 (FP32/FP32) or 4.73 (FP32/FP16S) or 16.12 (FP32/FP16C) FLOPs/Byte. So performance is only limited by memory bandwidth.

| Device                        | FP32<br>[TFLOPs/s] | Mem<br>[GB] | BW<br>[GB/s] | FP32/FP32<br>[MLUPs/s] | FP32/FP16S<br>[MLUPs/s] | FP32/FP16C<br>[MLUPs/s] |
| :---------------------------- | -----------------: | ----------: | -----------: | ---------------------: | ----------------------: | ----------------------: |
| AMD Instinct MI250 (1 GCD)    |              45.26 |          64 |         1638 |             5638 (53%) |              9030 (42%) |              8506 (40%) |
| AMD Radeon VII                |              13.83 |          16 |         1024 |             4898 (73%) |              7778 (58%) |              5256 (40%) |
| Nvidia H100 PCIe 80GB         |              51.01 |          80 |         2000 |            11128 (85%) |             20624 (79%) |             13862 (53%) |
| Nvidia A100 SXM4 80GB         |              19.49 |          80 |         2039 |            10228 (77%) |             18448 (70%) |             11197 (42%) |
| Nvidia A100 SXM4 40GB         |              19.49 |          40 |         1555 |             8522 (84%) |             16013 (79%) |             11251 (56%) |
| Nvidia A100 PCIe 40GB         |              19.49 |          40 |         1555 |             8526 (84%) |             16035 (79%) |             11088 (55%) |
| Nvidia Tesla V100 16GB        |              14.13 |          16 |          900 |             5128 (87%) |             10325 (88%) |              7683 (66%) |
| Nvidia Quadro GV100           |              16.66 |          32 |          870 |             3442 (61%) |              6641 (59%) |              5863 (52%) |
| Nvidia Tesla P100 16GB        |               9.52 |          16 |          732 |             3295 (69%) |              5950 (63%) |              4176 (44%) |
| Nvidia Tesla P100 12GB        |               9.52 |          12 |          549 |             2427 (68%) |              4141 (58%) |              3999 (56%) |
| Nvidia Tesla K40m             |               4.29 |          12 |          288 |             1131 (60%) |              1868 (50%) |               912 (24%) |
| Nvidia Tesla K80  (1 GPU)     |               4.11 |          12 |          240 |              916 (58%) |              1642 (53%) |               943 (30%) |
| Nvidia Tesla K20c             |               3.52 |           5 |          208 |              861 (63%) |              1507 (56%) |               720 (27%) |
| AMD Radeon RX 6900 XT         |              20.63 |          16 |          512 |             1968 (59%) |              4227 (64%) |              4207 (63%) |
| AMD Radeon RX 5700 XT         |               9.75 |           8 |          448 |             1368 (47%) |              3253 (56%) |              3049 (52%) |
| AMD Radeon RX Vega 64         |              13.35 |           8 |          484 |             1875 (59%) |              2878 (46%) |              3227 (51%) |
| AMD Radeon RX 580 4GB         |               6.50 |           4 |          256 |              946 (57%) |              1848 (56%) |              1577 (47%) |
| AMD Radeon HD 7850            |               1.84 |           2 |          154 |              112 (11%) |               120 ( 6%) |               635 (32%) |
| Intel Arc A770 LE             |              19.66 |          16 |          560 |             2741 (75%) |              4591 (63%) |              4626 (64%) |
| Intel Arc A750 LE             |              17.20 |           8 |          512 |             2625 (78%) |              4184 (63%) |              4238 (64%) |
| Nvidia GeForce RTX 4090       |              82.58 |          24 |         1008 |             5624 (85%) |             11091 (85%) |             11496 (88%) |
| Nvidia GeForce RTX 3090 Ti    |              40.00 |          24 |         1008 |             5717 (87%) |             10956 (84%) |             10400 (79%) |
| Nvidia GeForce RTX 3090       |              39.05 |          24 |          936 |             5418 (89%) |             10732 (88%) |             10215 (84%) |
| Nvidia GeForce RTX 3080 Ti    |              37.17 |          12 |          912 |             5202 (87%) |              9832 (87%) |              9347 (79%) |
| Nvidia GeForce RTX 3080       |              29.77 |          10 |          760 |             4230 (85%) |              8118 (82%) |              7714 (78%) |
| Nvidia GeForce RTX 3070       |              20.31 |           8 |          448 |             2578 (88%) |              5096 (88%) |              5060 (87%) |
| Nvidia GeForce RTX 3060 Ti    |              16.49 |           8 |          448 |             2644 (90%) |              5129 (88%) |              4718 (81%) |
| Nvidia RTX A5000M             |              16.59 |          16 |          448 |             2228 (76%) |              4461 (77%) |              3662 (63%) |
| Nvidia GeForce RTX 3060       |              13.17 |          12 |          360 |             2108 (90%) |              4070 (87%) |              3566 (76%) |
| Nvidia GeForce RTX 3060M      |              10.94 |           6 |          336 |             2019 (92%) |              4012 (92%) |              3572 (82%) |
| Nvidia GeForce RTX 3050M      |               7.13 |           4 |          192 |             1180 (94%) |              2339 (94%) |              2016 (81%) |
| Nvidia Quadro RTX 6000        |              16.31 |          24 |          672 |             3307 (75%) |              6836 (78%) |              6879 (79%) |
| Nvidia Quadro RTX 8000 Pass.  |              14.93 |          48 |          624 |             2591 (64%) |              5408 (67%) |              5607 (69%) |
| Nvidia GeForce RTX 2080 Ti    |              13.45 |          11 |          616 |             3194 (79%) |              6700 (84%) |              6853 (86%) |
| Nvidia GeForce RTX 2080 Sup.  |              11.34 |           8 |          496 |             2434 (75%) |              5284 (82%) |              5087 (79%) |
| Nvidia Quadro RTX 5000        |              11.15 |          16 |          448 |             2341 (80%) |              4766 (82%) |              4773 (82%) |
| Nvidia GeForce RTX 2060 Sup.  |               7.18 |           8 |          448 |             2503 (85%) |              5035 (87%) |              4463 (77%) |
| Nvidia Quadro RTX 4000        |               7.12 |           8 |          416 |             2284 (84%) |              4584 (85%) |              4062 (75%) |
| Nvidia GeForce RTX 2060 KO    |               6.74 |           6 |          336 |             1643 (75%) |              3376 (77%) |              3266 (75%) |
| Nvidia GeForce RTX 2060       |               6.74 |           6 |          336 |             1681 (77%) |              3604 (83%) |              3571 (82%) |
| Nvidia GeForce GTX 1660 Sup.  |               5.03 |           6 |          336 |             1696 (77%) |              3551 (81%) |              3040 (70%) |
| Nvidia Tesla T4               |               8.14 |          15 |          300 |             1356 (69%) |              2869 (74%) |              2887 (74%) |
| Nvidia GeForce GTX 1660 Ti    |               5.48 |           6 |          288 |             1467 (78%) |              3041 (81%) |              3019 (81%) |
| Nvidia GeForce GTX 1660       |               5.07 |           6 |          192 |             1016 (81%) |              1924 (77%) |              1992 (80%) |
| Nvidia GeForce GTX 1650M      |               3.20 |           4 |          128 |              706 (84%) |              1214 (73%) |              1400 (84%) |
| Nvidia Titan Xp               |              12.15 |          12 |          548 |             2919 (82%) |              5495 (77%) |              5375 (76%) |
| Nvidia GeForce GTX 1080 Ti    |              12.06 |          11 |          484 |             2631 (83%) |              4837 (77%) |              4877 (78%) |
| Nvidia GeForce GTX 1080       |               9.78 |           8 |          320 |             1623 (78%) |              3100 (75%) |              3182 (77%) |
| Nvidia GeForce GTX 1060M      |               4.44 |           6 |          192 |              983 (78%) |              1882 (75%) |              1803 (72%) |
| Nvidia GeForce GTX 1050M Ti   |               2.49 |           4 |          112 |              631 (86%) |              1224 (84%) |              1115 (77%) |
| Nvidia Quadro P1000           |               1.89 |           4 |           82 |              426 (79%) |               839 (79%) |               778 (73%) |
| Nvidia Quadro M4000           |               2.57 |           8 |          192 |              899 (72%) |              1519 (61%) |              1050 (42%) |
| Nvidia Tesla M60 (1 GPU)      |               4.82 |           8 |          160 |              853 (82%) |              1571 (76%) |              1557 (75%) |
| Nvidia GeForce GTX 960M       |               1.51 |           4 |           80 |              442 (84%) |               872 (84%) |               627 (60%) |
| Nvidia Quadro K2000           |               0.73 |           2 |           64 |              312 (75%) |               444 (53%) |               171 (21%) |
| Nvidia GeForce GT 630 (OEM)   |               0.46 |           2 |           29 |              151 (81%) |               185 (50%) |                78 (21%) |
| Nvidia Quadro NVS 290         |               0.03 |       0.256 |            6 |                1 ( 2%) |                 1 ( 1%) |                 1 ( 1%) |
| Apple M1 Pro GPU 16C 16GB     |               4.10 |          11 |          200 |             1204 (92%) |              2329 (90%) |              1855 (71%) |
| AMD Radeon Vega 8 (4750G)     |               2.15 |          27 |           57 |              263 (71%) |               511 (70%) |               501 (68%) |
| AMD Radeon Vega 8 (3500U)     |               1.23 |           7 |           38 |              157 (63%) |               282 (57%) |               288 (58%) |
| Intel UHD Graphics 630        |               0.46 |           7 |           51 |              151 (45%) |               301 (45%) |               187 (28%) |
| Intel HD Graphics 5500        |               0.35 |           3 |           26 |               75 (45%) |               192 (58%) |               108 (32%) |
| Intel HD Graphics 4600        |               0.38 |           2 |           26 |              105 (63%) |               115 (35%) |                34 (10%) |
| Samsung ARM Mali-G72 MP18     |               0.24 |           4 |           29 |               14 ( 7%) |                17 ( 5%) |                12 ( 3%) |
| 2x AMD EPYC 9654              |              29.49 |        1536 |          922 |             1381 (23%) |              1814 (15%) |              1801 (15%) |
| Intel Xeon Phi 7210           |               5.32 |         192 |          102 |              415 (62%) |               193 (15%) |               223 (17%) |
| 4x Intel Xeon E5-4620 v4      |               2.69 |         512 |          273 |              460 (26%) |               275 ( 8%) |               239 ( 7%) |
| 2x Intel Xeon E5-2630 v4      |               1.41 |          64 |          137 |              264 (30%) |               146 ( 8%) |               129 ( 7%) |
| 2x Intel Xeon E5-2623 v4      |               0.67 |          64 |          137 |              125 (14%) |                66 ( 4%) |                59 ( 3%) |
| 2x Intel Xeon E5-2680 v3      |               1.92 |          64 |          137 |              209 (23%) |               305 (17%) |               281 (16%) |
| Intel Core i9-10980XE         |               3.23 |         128 |           94 |              286 (47%) |               251 (21%) |               223 (18%) |
| Intel Core i5-9600            |               0.60 |          16 |           43 |              146 (52%) |               127 (23%) |               147 (27%) |
| Intel Core i7-8700K           |               0.71 |          16 |           51 |              152 (45%) |               134 (20%) |               116 (17%) |
| Intel Core i7-7700HQ          |               0.36 |          12 |           38 |               81 (32%) |                82 (16%) |               108 (22%) |
| Intel Core i7-4770            |               0.44 |          16 |           26 |              104 (62%) |                69 (21%) |                59 (18%) |
| Intel Core i7-4720HQ          |               0.33 |          16 |           26 |               58 (35%) |                13 ( 4%) |                47 (14%) |



## Maximum Grid Resolution for D3Q19 LBM

| Memory | FP32/FP32 | FP32/FP16 |
| -----: | --------: | --------: |
|   1 GB |      224³ |      264³ |
|   2 GB |      280³ |      336³ |
|   3 GB |      320³ |      384³ |
|   4 GB |      352³ |      424³ |
|   6 GB |      404³ |      484³ |
|   8 GB |      448³ |      532³ |
|  10 GB |      480³ |      572³ |
|  11 GB |      496³ |      592³ |
|  12 GB |      512³ |      608³ |
|  16 GB |      564³ |      672³ |
|  24 GB |      644³ |      768³ |
|  32 GB |      708³ |      848³ |
|  40 GB |      764³ |      912³ |
|  48 GB |      812³ |      968³ |
|  64 GB |      896³ |     1068³ |
|  80 GB |      964³ |     1148³ |
|  96 GB |     1024³ |     1220³ |
| 128 GB |     1128³ |     1344³ |
| 192 GB |     1292³ |     1540³ |
| 256 GB |     1420³ |     1624³ |
| 384 GB |     1624³ |     1624³ |



## FAQs

### General

- <details><summary>What physical model does FluidX3D use?</summary><br>FluidX3D implements the lattice Boltzmann method, a type of direct numerical simulation (DNS), the most accurate type of fluid simulation, but also the most computationally challenging. Optional extension models include volume force (Guo forcing), free surface (<a href="https://doi.org/10.3390/computation10060092">volume-of-fluid</a> and <a href="https://doi.org/10.3390/computation10020021">PLIC</a>), a temperature model and Smagorinsky-Lilly subgrid turbulence model.<br><br></details>

- <details><summary>FluidX3D only uses FP32 or even FP32/FP16, in contrast to FP64. Are simulation results physically accurate?</summary><br>Yes, in all but extreme edge cases. The code has been specially optimized to minimize arithmetic round-off errors and make the most out of lower precision. With these optimizations, accuracy in most cases is indistinguishable from FP64 double-precision, even with FP32/FP16 mixed-precision. Details can be found in <a href="https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats">this paper</a>.<br><br></details>

- <details><summary>Why is the simulation box size limited to 2³² grid points?</summary><br>The 32-bit unsigned integer grid index will overflow above this number. Using 64-bit index calculation would slow the simulation down by ~20%, as 64-bit uint is calculated on special function units and not the regular GPU cores. 2³² grid points with FP32/FP16 mixed-precision is equivalent to 225GB memory and single GPUs currently are only at 80GB, so it should be fine for a while to come.<br><br></details>

- <details><summary>Comparted to the benchmark numbers stated <a href="https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats">here</a>, efficiency seems much lower but performance is slightly better for most devices. How can this be?</summary><br>In that paper, the One-Step-Pull swap algorithm is implemented, using only misaligned reads and coalesced writes. On almost all GPUs, the performance penalty for misaligned writes is much larger than for misaligned reads, and sometimes there is almost no penalty for misaligned reads at all. Because of this, One-Step-Pull runs at peak bandwidth and thus peak efficiency.<br>Here, a different swap algorithm termed <a href="https://doi.org/10.3390/computation10060092">Esoteric-Pull</a> is used, a type of in-place streaming. This makes the LBM require much less memory (93 vs. 169 (FP32/FP32) or 55 vs. 93 (FP32/FP16) Bytes/node for D3Q19), and also less memory bandwidth (153 vs. 171 (FP32/FP32) or 77 vs. 95 (FP32/FP16) Bytes/node per time step for D3Q19) due to so-called implicit bounce-back boundaries. However memory access now is half coalesced and half misaligned for both reads and writes, so memory access efficiency is lower. For overall performance, these two effects approximately cancel out. The benefit of Esoteric-Pull - being able to simulate domains twice as large with the same amount of memory - clearly outweights the cost of slightly lower memory access efficiency, especially since performance is not reduced overall.<br><br></details>

- <details><summary>Why don't you use CUDA? Wouldn't that be more efficient?</summary><br>No, that is a wrong myth. OpenCL is exactly as efficient as CUDA on Nvidia GPUs if optimized properly. <a href="https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats">Here</a> I did roofline model and analyzed OpenCL performance on various hardware. OpenCL efficiency on modern Nvidia GPUs can be 100% with the right memory access pattern, so CUDA can't possibly be any more efficient. Without any performance advantage, there is no reason to use proprietary CUDA over OpenCL, since OpenCL is compatible with a lot more hardware.<br><br></details>

- <details><summary>Why no multi-relaxation-time (MRT) collision operator?</summary><br>The idea of MRT is to linearly transform the DDFs into "moment space" by matrix multiplication and relax these moments individually, promising better stability and accuracy. In practice, in the vast majority of cases, it has zero or even negative effects on stability and accuracy, and simple SRT is much superior. Apart from the kinematic shear viscosity and conserved terms, the remaining moments are non-physical quantities and their tuning is a blackbox. Although MRT can be implemented in an efficient manner with only a single matrix-vector multiplication in registers, leading to identical performance compared to SRT by remaining bandwidth-bound, storing the matrices vastly elongates and over-complicates the code for no real benefit.<br><br></details>

- <details><summary>How about multi-GPU?</summary><br>Multi-GPU would allow for larger simulation domains, but is currently not yet supported in FluidX3D. It makes code development and maintaining much more effort and many extensions close to impossible, for example ray-tracing graphics. With the in-place streaming implementation and FP16 memory compression combined, FluidX3D already only requires 1/3 the memory of existing LBM GPU implementations, so large domain size is already possible on a single GPU. Nevertheless, I will look into multi-GPU support in the future, so stay tuned.<br><br></details>

### Hardware

- <details><summary>I'm on a budget and have only a cheap computer. Can I run FluidX3D on my toaster PC/laptop?</summary><br>Absolutely. Today even the most inexpensive hardware, like integrated GPUs or entry-level gaming GPUs, support OpenCL. You might be a bit more limited on memory capacity and grid resolution, but you should be good to go. I've tested FluidX3D on very old and inexpensive hardware and even on my Samsung S9+ smartphone, and it runs just fine, although admittedly a bit slower.<br><br></details>

- <details><summary>I don't have an expensive workstation GPU, but only a gaming GPU. Will performance suffer?</summary><br>No. Efficiency on gaming GPUs is exactly as good as on their "professional"/workstation counterparts. Performance often is even better as gaming GPUs have higher boost clocks.<br><br></details>

- <details><summary>Do I need a GPU with ECC memory?</summary><br>No. Gaming GPUs work just fine. Some Nvidia GPUs automatically reduce memory clocks for compute applications to almost entirely eliminate memory errors.<br><br></details>

- <details><summary>My GPU does not support CUDA. Can I still use FluidX3D?</summary><br>Yes. FluidX3D uses OpenCL 1.2 and not CUDA, so it runs on any GPU from any vendor since around 2012.<br><br></details>

- <details><summary>I don't have a dedicated graphics card at all. Can I still run FluidX3D on my PC/laptop?</summary><br>Yes. FluidX3D also runs on all integrated GPUs since around 2012, and also on CPUs.<br><br></details>

- <details><summary>I need more memory than my GPU can offer. Can I run FluidX3D on my CPU as well?</summary><br>Yes. You only need to install the <a href="https://www.intel.com/content/www/us/en/developer/articles/tool/opencl-drivers.html">OpenCL Runtime for Intel CPUs</a>.<br><br></details>

- <details><summary>In the benchmarks you list some very expensive hardware. How do you get access to that?</summary><br>I'm a scientist (PhD candidate in computational physics) and I use FluidX3D for my research, so I have access to BZHPC, SuperMUC-NG and JURECA-DC supercomputers.<br><br></details>

### Graphics

- <details><summary>I don't have an RTX/DXR GPU that supports raytracing. Can I still use raytracing graphics in FluidX3D?</summary><br>Yes, and at full performance. FluidX3D does not use a bounding volume hierarchy (BVH) to accelerate raytracing, but fast ray-grid traversal instead, implemented directly in OpenCL C. This is much faster than BVH for moving isosurfaces in the LBM grid (~N vs. ~N²+log(N) runtime; LBM itself is ~N³), and it does not require any dedicated raytracing hardware. Raytracing in FluidX3D runs on any GPU that supports OpenCL 1.2.<br><br></details>

- <details><summary>I have a datacenter/mining GPU without any video output or graphics hardware. Can FluidX3D still render simulation results?</summary><br>Yes. FluidX3D does all rendering (rasterization and raytracing) in OpenCL C, so no display output and no graphics features like OpenGL/Vulkan/DirectX are required. Rendering is just another form of compute after all. Rendered frames are passed to the CPU over PCIe and then the CPU can either draw them on screen through dedicated/integrated graphics or write them to the hard drive.<br><br></details>

- <details><summary>I'm running FluidX3D on a remote (super-)computer and only have an SSH terminal. Can I still use graphics somehow?</summary><br>Yes, either directly as interactive ASCII graphics in the terminal or by storing rendered frames on the hard drive and then copying them over via `scp -r user@server.url:"~/path/to/images/folder" .`.<br><br></details>

### Licensing

- <details><summary>I want to learn about programming/software/physics/engineering. Can I use FluidX3D for free?</summary><br>Yes. Anyone can use FluidX3D for free for public research, education or personal use. Use by scientists, students and hobbyists is free of charge and well encouraged.<br><br></details>

- <details><summary>I am a scientist/teacher with a paid position at a public institution. Can I use FluidX3D for my research/teaching?</summary><br>Yes, you can use FluidX3D free of charge. This is considered research/education, not commercial use. To give credit, the <a href="https://github.com/ProjectPhysX/FluidX3D#references">references</a> listed below should be cited. If you publish data/results generated by altered source versions, the altered source code must be published as well.<br><br></details>

- <details><summary>I work at a company in CFD/consulting/R&D or related fields. Can I use FluidX3D commercially?</summary><br>No. Commercial use is not allowed with the current license.<br><br></details>

- <details><summary>Is FluidX3D open-source?</summary><br>No. "Open-source" as a technical term is defined as freely available without any restriction on use, but I am not comfortable with that. I have written FluidX3D in my spare time and no one should milk it for profits while I remain uncompensated, especially considering what other CFD software sells for. The technical term for the type of license I choose is "source-available no-cost non-commercial". The source code is freely available, and you are free to use, to alter and to redistribute it, as long as you do not sell it or make a profit from derived products/services, and as long as you do not use it for any military purposes (see the <a href="https://github.com/ProjectPhysX/FluidX3D/blob/master/LICENSE.md">license</a> for details).<br><br></details>

- <details><summary>Will FluidX3D at some point be available with a commercial license?</summary><br>Maybe I will add the option for a second, commercial license later on. If you are interested in commercial use, let me know. For non-commercial use in science and education, FluidX3D is and will always be free.<br><br></details>



## External Code/Libraries/Images used in FluidX3D

- [OpenCL-Headers](https://github.com/KhronosGroup/OpenCL-Headers) for GPU parallelization ([Khronos Group](https://www.khronos.org/opencl/))
- [Win32 API](https://learn.microsoft.com/en-us/windows/win32/api/winbase/) for interactive graphics in Windows ([Microsoft](https://www.microsoft.com/))
- [X11/Xlib](https://www.x.org/releases/current/doc/libX11/libX11/libX11.html) for interactive graphics in Linux ([The Open Group](https://www.x.org/releases/current/doc/libX11/libX11/libX11.html))
- [marching-cubes tables](http://paulbourke.net/geometry/polygonise/) for isosurface generation on GPU ([Paul Bourke](http://paulbourke.net/geometry/))
- [`src/lodepng.cpp`](https://github.com/lvandeve/lodepng/blob/master/lodepng.cpp) and [`src/lodepng.hpp`](https://github.com/lvandeve/lodepng/blob/master/lodepng.h) for `.png` encoding and decoding ([Lode Vandevenne](https://lodev.org/))
- [SimplexNoise](https://weber.itn.liu.se/~stegu/simplexnoise/SimplexNoise.java) class in [`src/utilities.hpp`](https://github.com/ProjectPhysX/FluidX3D/blob/master/src/utilities.hpp) for generating continuous noise in 2D/3D/4D space ([Stefan Gustavson](https://github.com/stegu))
- [`skybox/skybox8k.png`](https://www.hdri-hub.com/hdri-skies-aviation-aerospace) for free surface raytracing ([HDRI Hub](https://www.hdri-hub.com/))



## References

- Lehmann, M.: [Esoteric Pull and Esoteric Push: Two Simple In-Place Streaming Schemes for the Lattice Boltzmann Method on GPUs](https://doi.org/10.3390/computation10060092). Computation, 10, 92, (2022)
- Lehmann, M., Krause, M., Amati, G., Sega, M., Harting, J. and Gekle, S.: [Accuracy and performance of the lattice Boltzmann method with 64-bit, 32-bit, and customized 16-bit number formats](https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats). Phys. Rev. E 106, 015308, (2022)
- Lehmann, M.: [Combined scientific CFD simulation and interactive raytracing with OpenCL](https://www.researchgate.net/publication/360501260_Combined_scientific_CFD_simulation_and_interactive_raytracing_with_OpenCL). IWOCL'22: International Workshop on OpenCL, 3, 1-2, (2022)
- Lehmann, M., Oehlschlägel, L.M., Häusl, F., Held, A. and Gekle, S.: [Ejection of marine microplastics by raindrops: a computational and experimental study](https://doi.org/10.1186/s43591-021-00018-8). Micropl.&Nanopl. 1, 18, (2021)
- Lehmann, M.: [High Performance Free Surface LBM on GPUs](https://doi.org/10.15495/EPub_UBT_00005400). Master's thesis, (2019)
- Lehmann, M. and Gekle, S.: [Analytic Solution to the Piecewise Linear Interface Construction Problem and Its Application in Curvature Calculation for Volume-of-Fluid Simulation Codes](https://doi.org/10.3390/computation10020021). Computation, 10, 21, (2022)



## Contact

- For any questions, feedback or other inquiries, don't hesitate to contact me at [moritz.lehmann@uni-bayreuth.de](mailto:moritz.lehmann@uni-bayreuth.de?subject=FluidX3D).
- Updates will be posted on Twitter via [@FluidX3D](https://twitter.com/FluidX3D) and [@ProjectPhysX](https://twitter.com/ProjectPhysX), under the hashtag [#FluidX3D](https://twitter.com/hashtag/FluidX3D?src=hashtag_click&f=live) or on my [YouTube channel](https://www.youtube.com/c/ProjectPhysX).