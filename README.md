# FluidX3D

The fastest and most memory efficient lattice Boltzmann CFD software, running on all GPUs via [OpenCL](https://github.com/ProjectPhysX/OpenCL-Wrapper "OpenCL-Wrapper").

<a href="https://youtu.be/-MkRBeQkLk8"><img src="https://img.youtube.com/vi/o3TPN142HxM/maxresdefault.jpg" width="50%"></img></a><a href="https://youtu.be/oC6U1M0Fsug"><img src="https://img.youtube.com/vi/oC6U1M0Fsug/maxresdefault.jpg" width="50%"></img></a><br>
<a href="https://youtu.be/XOfXHgP4jnQ"><img src="https://img.youtube.com/vi/XOfXHgP4jnQ/maxresdefault.jpg" width="50%"></img></a><a href="https://youtu.be/3JNVBQyetMA"><img src="https://img.youtube.com/vi/3JNVBQyetMA/maxresdefault.jpg" width="50%"></img></a>
(click on images to show videos on YouTube)


<details><summary>Update History</summary>

- v1.0 (04.08.2022)
  - initial release
- v1.1 (29.09.2022)
  - added solid voxelization on GPU (slow algorithm)
  - added tool to print current camera position (key_H)
  - minor bug fix (workaround for Intel iGPU driver bug with triangle rendering)
- v1.2 (24.10.2022)
  - added functions to compute force/torque on objects
  - added function to translate Mesh
  - added Stokes drag validation setup
- v1.3 (10.11.2022)
  - added unit conversion functions for torque
  - `FORCE_FIELD` and `VOLUME_FORCE` can now be used independently
  - minor bug fix (workaround for AMD legacy driver bug with binary number literals)
- v1.4 (14.12.2022)
  - added interactive graphics mode on Linux with X11
  - fixed streamline visualization bug in 2D
- v2.0 (09.01.2023)
  - added (cross-vendor) multi-GPU support on a single node (PC/laptop/server)
- v2.1 (15.01.2023)
  - made solid voxelization on GPU lightning fast (new algorithm, from minutes to milliseconds)
- v2.2 (20.01.2023)
  - added option to voxelize moving/rotating geometry on GPU, with automatic velocity initialization for each grid point based on center of rotation, linear velocity and rotational velocity
  - cells that are converted from solid->fluid during re-voxelization now have their DDFs properly initialized
  - added option to not auto-scale mesh during `read_stl(...)`, with negative `size` parameter
  - added kernel for solid boundary rendering with marching-cubes
- v2.3 (30.01.2023)
  - added particles with immersed-boundary method (either passive or 2-way-coupled, only supported with single-GPU)
  - minor optimization to GPU voxelization algorithm (workgroup threads outside mesh bounding-box return after ray-mesh intersections have been found)
  - displayed GPU memory allocation size is now fully accurate
  - fixed bug in `write_line()` function in `src/utilities.hpp`
  - removed `.exe` file extension for Linux/macOS
- v2.4 (11.03.2023)
  - added a help menu with key H that shows keyboard/mouse controls, visualization settings and simulation stats
  - improvements to keyboard/mouse control (+/- for zoom, mouseclick frees/locks cursor)
  - added suggestion of largest possible grid resolution if resolution is set larger than memory allows
  - minor optimizations in multi-GPU communication (insignificant performance difference)
  - fixed bug in temperature equilibrium function for temperature extension
  - fixed erroneous double literal for Intel iGPUs in skybox color functions
  - fixed bug in make.sh where multi-GPU device IDs would not get forwarded to the executable
  - minor bug fixes in graphics engine (free cursor not centered during rotation, labels in VR mode)
  - fixed bug in LBM::voxelize_stl() size parameter standard initialization
- v2.5 (11.04.2023)
  - implemented light absorption in fluid for raytracing graphics (no performance impact)
  - improved raytracing framerate when camera is inside fluid
  - fixed skybox pole flickering artifacts
  - fixed bug where moving objects during re-voxelization would leave an erroneous trail of solid grid cells behind
- v2.6 (16.04.2023)
  - patched OpenCL issues of Intel Arc GPUs: now VRAM allocations >4GB are possible and correct VRAM capacity is reported
- v2.7 (29.05.2023)
  - added slice visualization (key <kbd>2</kbd> / key <kbd>3</kbd> modes, then switch through slice modes with key <kbd>T</kbd>, move slice with keys <kbd>Q</kbd>/<kbd>E</kbd>)
  - made flag wireframe / solid surface visualization kernels toggleable with key <kbd>1</kbd>
  - added surface pressure visualization (key <kbd>1</kbd> when `FORCE_FIELD` is enabled and `lbm.calculate_force_on_boundaries();` is called)
  - added binary `.vtk` export function for meshes with `lbm.write_mesh_to_vtk(Mesh* mesh);`
  - added `time_step_multiplicator` for `integrate_particles()` function in PARTICLES extension
  - made correction of wrong memory reporting on Intel Arc more robust
  - fixed bug in `write_file()` template functions
  - reverted back to separate `cl::Context` for each OpenCL device, as the shared Context otherwise would allocate extra VRAM on all other unused Nvidia GPUs
  - removed Debug and x86 configurations from Visual Studio solution file (one less complication for compiling)
  - fixed bug that particles could get too close to walls and get stuck, or leave the fluid phase (added boundary force)

</details>


## Compute Features - Getting the Memory Problem under Control

- <details><summary>CFD model: lattice Boltzmann method (LBM)</summary>

  - streaming (part 2/2)<p align="center"><i>f</i><sub>0</sub><sup>temp</sup>(<i>x</i>,<i>t</i>) = <i>f</i><sub>0</sub>(<i>x</i>, <i>t</i>)<br><i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>) = <i>f</i><sub>(<i>t</i>%2 ? <i>i</i> : (<i>i</i>%2 ? <i>i</i>+1 : <i>i</i>-1))</sub>(<i>i</i>%2 ? <i>x</i> : <i>x</i>-<i>e<sub>i</sub></i>, <i>t</i>) &nbsp; for &nbsp; <i>i</i> &isin; [1, <i>q</i>-1]</p>
  - collision<p align="center"><i>&rho;</i>(<i>x</i>,<i>t</i>) = (&Sigma;<sub><i>i</i></sub> <i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>)) + 1<br><br><i>u</i>(<i>x</i>,<i>t</i>) = <sup>1</sup>&#8725;<sub><i>&rho;</i>(<i>x</i>,<i>t</i>)</sub> &Sigma;<sub><i>i</i></sub> <i>c<sub>i</sub></i> <i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>)<br><br><i>f<sub>i</sub></i><sup>eq-shifted</sup>(<i>x</i>,<i>t</i>) = <i>w<sub>i</sub></i> <i>&rho;</i> · (<sup>(<i>u</i><sub>°</sub><i>c<sub>i</sub></i>)<sup>2</sup></sup>&#8725;<sub>(2<i>c</i><sup>4</sup>)</sub> - <sup>(<i>u</i><sub>°</sub><i>u</i>)</sup>&#8725;<sub>(2c<sup>2</sup>)</sub> + <sup>(<i>u</i><sub>°</sub><i>c<sub>i</sub></i>)</sup>&#8725;<sub><i>c</i><sup>2</sup></sub>) + <i>w<sub>i</sub></i> (<i>&rho;</i>-1)<br><br><i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>, <i>t</i>+&Delta;<i>t</i>) = <i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>) + <i>&Omega;<sub>i</sub></i>(<i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>,<i>t</i>), <i>f<sub>i</sub></i><sup>eq-shifted</sup>(<i>x</i>,<i>t</i>), <i>&tau;</i>)</p>
  - streaming (part 1/2)<p align="center"><i>f</i><sub>0</sub>(<i>x</i>, <i>t</i>+&Delta;<i>t</i>) = <i>f</i><sub>0</sub><sup>temp</sup>(<i>x</i>, <i>t</i>+&Delta;<i>t</i>)<br><i>f</i><sub>(<i>t</i>%2 ? (<i>i</i>%2 ? <i>i</i>+1 : <i>i</i>-1) : <i>i</i>)</sub>(<i>i</i>%2 ? <i>x</i>+<i>e<sub>i</sub></i> : <i>x</i>, <i>t</i>+&Delta;<i>t</i>) = <i>f<sub>i</sub></i><sup>temp</sup>(<i>x</i>, <i>t</i>+&Delta;<i>t</i>) &nbsp; for &nbsp; <i>i</i> &isin; [1, <i>q</i>-1]</p>
  - velocity sets: D2Q9, D3Q15, D3Q19 (default), D3Q27
  - collision operators: single-relaxation-time (SRT/BGK) (default), two-relaxation-time (TRT)

  </details>

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

- <details><summary>optimized to minimize VRAM footprint to 1/6 of other LBM codes</summary>

  - traditional LBM (D3Q19) with FP64 requires ~344 Bytes/cell<br>
    - 🟧🟧🟧🟧🟧🟧🟧🟧🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟨🟨🟨🟨🟨🟨🟨🟨🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥🟥<br>(density 🟧, velocity 🟦, flags 🟨, 2 copies of DDFs 🟩/🟥; each square = 1 Byte)
    - allows for 3 Million cells per 1 GB VRAM
  - FluidX3D (D3Q19) requires only 55 Bytes/cell with [Esoteric-Pull](https://doi.org/10.3390/computation10060092)+[FP16](https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats)<br>
    - 🟧🟧🟧🟧🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟨🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩<br>(density 🟧, velocity 🟦, flags 🟨, DDFs 🟩; each square = 1 Byte)
    - allows for 19 Million cells per 1 GB VRAM
    - in-place streaming with [Esoteric-Pull](https://doi.org/10.3390/computation10060092): eliminates redundant copy `B` of density distribution functions (DDFs) in memory; almost cuts memory demand in half and slightly increases performance due to implicit bounce-back boundaries; offers optimal memory access patterns for single-cell in-place streaming
    - [decoupled arithmetic precision (FP32) and memory precision (FP32 or FP16S or FP16C)](https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats): all arithmetic is done in FP32 for compatibility on all hardware, but DDFs in memory can be compressed to FP16S or FP16C: almost cuts memory demand in half again and almost doubles performance, without impacting overall accuracy for most setups
    - <details><summary>only 8 flag bits per lattice point (can be used independently / at the same time)</summary>

      - `TYPE_S` (stationary or moving) solid boundaries
      - `TYPE_E` equilibrium boundaries (inflow/outflow)
      - `TYPE_T` temperature boundaries
      - `TYPE_F` free surface (fluid)
      - `TYPE_I` free surface (interface)
      - `TYPE_G` free surface (gas)
      - `TYPE_X` remaining for custom use or further extensions
      - `TYPE_Y` remaining for custom use or further extensions

      </details>
  - large cost saving: comparison of maximum single-GPU grid resolution for D3Q19 LBM

    | GPU&nbsp;VRAM&nbsp;capacity      | 1&nbsp;GB | 2&nbsp;GB | 3&nbsp;GB | 4&nbsp;GB | 6&nbsp;GB | 8&nbsp;GB | 10&nbsp;GB | 11&nbsp;GB | 12&nbsp;GB | 16&nbsp;GB | 20&nbsp;GB | 24&nbsp;GB | 32&nbsp;GB | 40&nbsp;GB | 48&nbsp;GB | 64&nbsp;GB | 80&nbsp;GB | 94&nbsp;GB | 128&nbsp;GB | 192&nbsp;GB | 256&nbsp;GB |
    | :------------------------------- | --------: | --------: | --------: | --------: | --------: | --------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ----------: | ----------: | ----------: |
    | approximate&nbsp;GPU&nbsp;price  | $25<br>GT&nbsp;210 | $25<br>GTX&nbsp;950 | $12<br>GTX&nbsp;1060 | $50<br>GT&nbsp;730 | $35<br>GTX&nbsp;1060 | $70<br>RX&nbsp;470 | $500<br>RTX&nbsp;3080 | $240<br>GTX&nbsp;1080&nbsp;Ti | $75<br>Tesla&nbsp;M40 | $75<br>Instinct&nbsp;MI25 | $900<br>RX&nbsp;7900&nbsp;XT | $205<br>Tesla&nbsp;P40 | $600<br>Instinct&nbsp;MI60 | $5500<br>A100 | $2400<br>RTX&nbsp;8000 | $10k<br>Instinct&nbsp;MI210 | $11k<br>A100 | >$40k<br>H100&nbsp;NVL | ?<br>Max&nbsp;Series&nbsp;1550 | - | - |
    | traditional&nbsp;LBM&nbsp;(FP64) |      144³ |      182³ |      208³ |      230³ |      262³ |      288³ |       312³ |       322³ |       330³ |       364³ |       392³ |       418³ |       460³ |       494³ |       526³ |       578³ |       624³ |       658³ |        730³ |        836³ |        920³ |
    | FluidX3D&nbsp;(FP32/FP32)        |      224³ |      282³ |      322³ |      354³ |      406³ |      448³ |       482³ |       498³ |       512³ |       564³ |       608³ |       646³ |       710³ |       766³ |       814³ |       896³ |       966³ |      1018³ |       1130³ |       1292³ |       1422³ |
    | FluidX3D&nbsp;(FP32/FP16)        |      266³ |      336³ |      384³ |      424³ |      484³ |      534³ |       574³ |       594³ |       610³ |       672³ |       724³ |       770³ |       848³ |       912³ |       970³ |      1068³ |      1150³ |      1214³ |       1346³ |       1540³ |       1624³ |

  </details>
- <details><summary>cross-vendor multi-GPU support on a single PC/laptop/server</summary>

  - domain decomposition allows pooling VRAM from multiple GPUs for much larger grid resolution
  - each domain (GPU) can hold up to 4.29 billion (2³², 1624³) lattice points (225 GB memory)
  - GPUs don't have to be identical (not even from the same vendor), but similar VRAM capacity/bandwidth is recommended
  - domain communication architecture (simplified)
    ```diff
    ++   .-----------------------------------------------------------------.   ++
    ++   |                              GPU 0                              |   ++
    ++   |                          LBM Domain 0                           |   ++
    ++   '-----------------------------------------------------------------'   ++
    ++              |                 selective                /|\             ++
    ++             \|/               in-VRAM copy               |              ++
    ++        .-------------------------------------------------------.        ++
    ++        |               GPU 0 - Transfer Buffer 0               |        ++
    ++        '-------------------------------------------------------'        ++
    !!                            |     PCIe     /|\                           !!
    !!                           \|/    copy      |                            !!
    @@        .-------------------------.   .-------------------------.        @@
    @@        | CPU - Transfer Buffer 0 |   | CPU - Transfer Buffer 1 |        @@
    @@        '-------------------------'\ /'-------------------------'        @@
    @@                           pointer  X   swap                             @@
    @@        .-------------------------./ \.-------------------------.        @@
    @@        | CPU - Transfer Buffer 1 |   | CPU - Transfer Buffer 0 |        @@
    @@        '-------------------------'   '-------------------------'        @@
    !!                           /|\    PCIe      |                            !!
    !!                            |     copy     \|/                           !!
    ++        .-------------------------------------------------------.        ++
    ++        |               GPU 1 - Transfer Buffer 1               |        ++
    ++        '-------------------------------------------------------'        ++
    ++             /|\                selective                 |              ++
    ++              |                in-VRAM copy              \|/             ++
    ++   .-----------------------------------------------------------------.   ++
    ++   |                              GPU 1                              |   ++
    ++   |                          LBM Domain 1                           |   ++
    ++   '-----------------------------------------------------------------'   ++
    ##                                    |                                    ##
    ##                      domain synchronization barrier                     ##
    ##                                    |                                    ##
    ||   -------------------------------------------------------------> time   ||
    ```
  - domain communication architecture (detailed)
    ```diff
    ++   .-----------------------------------------------------------------.   ++
    ++   |                              GPU 0                              |   ++
    ++   |                          LBM Domain 0                           |   ++
    ++   '-----------------------------------------------------------------'   ++
    ++     |  selective in- /|\  |  selective in- /|\  |  selective in- /|\    ++
    ++    \|/ VRAM copy (X)  |  \|/ VRAM copy (Y)  |  \|/ VRAM copy (Z)  |     ++
    ++   .---------------------.---------------------.---------------------.   ++
    ++   |    GPU 0 - TB 0X+   |    GPU 0 - TB 0Y+   |    GPU 0 - TB 0Z+   |   ++
    ++   |    GPU 0 - TB 0X-   |    GPU 0 - TB 0Y-   |    GPU 0 - TB 0Z-   |   ++
    ++   '---------------------'---------------------'---------------------'   ++
    !!          | PCIe /|\            | PCIe /|\            | PCIe /|\         !!
    !!         \|/ copy |            \|/ copy |            \|/ copy |          !!
    @@   .---------. .---------.---------. .---------.---------. .---------.   @@
    @@   | CPU 0X+ | | CPU 1X- | CPU 0Y+ | | CPU 3Y- | CPU 0Z+ | | CPU 5Z- |   @@
    @@   | CPU 0X- | | CPU 2X+ | CPU 0Y- | | CPU 4Y+ | CPU 0Z- | | CPU 6Z+ |   @@
    @@   '---------\ /---------'---------\ /---------'---------\ /---------'   @@
    @@      pointer X swap (X)    pointer X swap (Y)    pointer X swap (Z)     @@
    @@   .---------/ \---------.---------/ \---------.---------/ \---------.   @@
    @@   | CPU 1X- | | CPU 0X+ | CPU 3Y- | | CPU 0Y+ | CPU 5Z- | | CPU 0Z+ |   @@
    @@   | CPU 2X+ | | CPU 0X- | CPU 4Y+ | | CPU 0Y- | CPU 6Z+ | | CPU 0Z- |   @@
    @@   '---------' '---------'---------' '---------'---------' '---------'   @@
    !!         /|\ PCIe |            /|\ PCIe |            /|\ PCIe |          !!
    !!          | copy \|/            | copy \|/            | copy \|/         !!
    ++   .--------------------..---------------------..--------------------.   ++
    ++   |   GPU 1 - TB 1X-   ||    GPU 3 - TB 3Y-   ||   GPU 5 - TB 5Z-   |   ++
    ++   :====================::=====================::====================:   ++
    ++   |   GPU 2 - TB 2X+   ||    GPU 4 - TB 4Y+   ||   GPU 6 - TB 6Z+   |   ++
    ++   '--------------------''---------------------''--------------------'   ++
    ++    /|\ selective in-  |  /|\ selective in-  |  /|\ selective in-  |     ++
    ++     |  VRAM copy (X) \|/  |  VRAM copy (Y) \|/  |  VRAM copy (Z) \|/    ++
    ++   .--------------------..---------------------..--------------------.   ++
    ++   |        GPU 1       ||        GPU 3        ||        GPU 5       |   ++
    ++   |    LBM Domain 1    ||    LBM Domain 3     ||    LBM Domain 5    |   ++
    ++   :====================::=====================::====================:   ++
    ++   |        GPU 2       ||        GPU 4        ||        GPU 6       |   ++
    ++   |    LBM Domain 2    ||    LBM Domain 4     ||    LBM Domain 6    |   ++
    ++   '--------------------''---------------------''--------------------'   ++
    ##              |                     |                     |              ##
    ##              |      domain synchronization barriers      |              ##
    ##              |                     |                     |              ##
    ||   -------------------------------------------------------------> time   ||
    ```

  </details>
- [peak performance on GPUs](#single-gpu-benchmarks) (datacenter/gaming/professional/laptop), validated with roofline model
- [DDF-shifting](https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats) and other algebraic optimization to minimize round-off error

- <details><summary>powerful model extensions</summary>

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
  - Smagorinsky-Lilly subgrid turbulence LES model to keep simulations with very large Reynolds number stable
    <p align="center"><i>&Pi;<sub>&alpha;&beta;</sub></i> = &Sigma;<sub><i>i</i></sub> <i>e<sub>i&alpha;</sub></i> <i>e<sub>i&beta;</sub></i> (<i>f<sub>i</sub></i>   - <i>f<sub>i</sub></i><sup>eq-shifted</sup>)<br><br>Q = &Sigma;<sub><i>&alpha;&beta;</i></sub>   <i>&Pi;<sub>&alpha;&beta;</sub></i><sup>2</sup><br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;______________________<br>&tau; = &frac12; (&tau;<sub>0</sub> + &radic; &tau;<sub>0</sub><sup>2</sup> + <sup>(16&radic;2)</sup>&#8725;<sub>(<i>3&pi;</i><sup>2</sup>)</sub> <sup>&radic;Q</sup>&#8725;<sub><i>&rho;</i></sub> )</p>
  - particles with immersed-boundary method (either passive or 2-way-coupled, only supported with single-GPU)

  </details>



## Solving the Visualization Problem

- FluidX3D can do simulations so large that storing the volumetric data for later rendering becomes unmanageable (like 120GB for a single frame, hundreds of TeraByte for a video)
- instead, FluidX3D allows [rendering raw simulation data directly in VRAM](https://www.researchgate.net/publication/360501260_Combined_scientific_CFD_simulation_and_interactive_raytracing_with_OpenCL), so no large volumetric files have to be exported to the hard disk (see my [technical talk](https://youtu.be/pD8JWAZ2f8o))
- the rendering is so fast that it works interactively in real time for both rasterization and raytracing
- if no monitor is available (like on a remote Linux server), there is an ASCII rendering mode to interactively visualize the simulation in the terminal (even in WSL and/or through SSH)
- rendering is fully multi-GPU-parallelized via seamless domain decomposition rasterization
- with interactive graphics mode disabled, image resolution can be as large as VRAM allows for (4K/8K/16K and above)
- (interacitive) visualization modes:
  - flag wireframe / solid surface (and force vectors on solid cells or surface pressure if the extension is used)
  - velocity field (with slice mode)
  - streamlines (with slice mode)
  - velocity-colored Q-criterion isosurface
  - rasterized free surface with [marching-cubes](http://paulbourke.net/geometry/polygonise/)
  - [raytraced free surface](https://www.researchgate.net/publication/360501260_Combined_scientific_CFD_simulation_and_interactive_raytracing_with_OpenCL) with fast ray-grid traversal and marching-cubes, either 1-4 rays/pixel or 1-10 rays/pixel



## Solving the Compatibility Problem

- FluidX3D is written in OpenCL 1.2, so it runs on any hardware from any vendor (Nvidia, AMD, Intel, ...):
  - world's fastest datacenter GPUs, like H100, A100, MI250(X), MI210, MI100, V100(S), P100, ...
  - gaming GPUs (desktop or laptop), like Nvidia GeForce, AMD Radeon, Intel Arc
  - professional/workstation GPUs, like Nvidia Quadro, AMD Radeon Pro / FirePro
  - integrated GPUs
  - Intel Xeon Phi (requires installation of the [Intel OpenCL CPU Runtime ("oclcpuexp")](https://github.com/intel/llvm/releases?q=oneAPI+DPC%2B%2B+Compiler))
  - Intel/AMD CPUs (requires installation of the [Intel OpenCL CPU Runtime ("oclcpuexp")](https://github.com/intel/llvm/releases?q=oneAPI+DPC%2B%2B+Compiler))
  - even smartphone ARM GPUs
- supports parallelization across multiple GPUs on a single PC/laptop/server with PCIe communication, no SLI/Crossfire/NVLink/InfinityFabric or MPI installation required; the GPUs don't even have to be from the same vendor, but similar memory capacity and bandwidth is recommended
- works in Windows and Linux with C++17, with limited support also for MacOS and Android
- supports importing and voxelizing triangle meshes from binary `.stl` files, with fast GPU voxelization
- supports exporting volumetric data as binary `.vtk` files with `lbm.<field>.write_device_to_vtk();`
- supports exporting rendered frames as `.png`/`.qoi`/`.bmp` files with `lbm.graphics.write_frame();`, encoding is handled in parallel on the CPU while the simulation on GPU can continue without delay



## How to get started?

1. Check the settings and extensions in [`src/defines.hpp`](src/defines.hpp) by uncommenting corresponding lines.
2. Write a C++ setup skript as `main_setup()` function in [`src/setup.cpp`](src/setup.cpp) (get inspiration from existing setups):
   - For unit conversion, use the `units` struct.
   - For initializing the box, use call `LBM lbm(Nx, Ny, Nz, nu, ...);` constructor. To use multiple GPUs, use `LBM lbm(Nx, Ny, Nz, Dx, Dy, Dz, nu, ...);`, with `Dx`/`Dy`/`Dz` indicating how many domains (GPUs) there are in each spatial direction.
   - Set the initial condition in a loop that iterates over the entire lattice by writing to `lbm.rho[n]`/`lbm.u.x[n]`/`lbm.u.y[n]`/`lbm.u.z[n]`/`lbm.flags[n]`.
   - Call `lbm.run();` to initialize and execute the setup (infinite time steps) or `lbm.run(time_steps);` to execute only a specific number of time steps.
   - As long as the `lbm` object is in scope, you can access the memory. As soon as it goes out of scope, all memory associated to the current simulation is freed again.
3. On Windows in Visual Studio Community click compile+run, or on Linux run `chmod +x make.sh` and `./make.sh`; this will automatically select the fastest installed GPU(s). Alternatively, you can add the device ID(s) as command-line arguments, for example `./make.sh 2` to compile+run on device 2, or `bin/FluidX3D 1 3` to run the executable on devices 1 and 3. Compile time for the entire code is about 10 seconds. If you use `INTERACTIVE_GRAPHICS` on Linux, change to the "compile on Linux with X11" command in `make.sh`.
4. Keyboard/mouse controls with `INTERACTIVE_GRAPHICS`/`INTERACTIVE_GRAPHICS_ASCII` enabled:
   - <kbd>P</kbd>: start/pause the simulation
   - <kbd>H</kbd>: show/hide help
   - <kbd>1</kbd>: flag wireframe / solid surface (and force vectors on solid cells or surface pressure if the extension is used)
   - <kbd>2</kbd>: velocity field
   - <kbd>3</kbd>: streamlines
   - <kbd>4</kbd>: vorticity / velocity-colored Q-criterion isosurface
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



## Single-GPU/CPU Benchmarks

Here are [performance benchmarks](https://doi.org/10.3390/computation10060092) on various hardware in MLUPs/s, or how many million lattice points are updated per second. The settings used for the benchmark are D3Q19 SRT with no extensions enabled (only LBM with implicit mid-grid bounce-back boundaries) and the setup consists of an empty cubic box with sufficient size (typically 256³). Without extensions, a single lattice point requires:
- a memory capacity of 93 (FP32/FP32) or 55 (FP32/FP16) Bytes
- a memory bandwidth of 153 (FP32/FP32) or 77 (FP32/FP16) Bytes per time step
- 363 (FP32/FP32) or 406 (FP32/FP16S) or 1275 (FP32/FP16C) FLOPs per time step (FP32+INT32 operations counted combined)

In consequence, the arithmetic intensity of this implementation is 2.37 (FP32/FP32) or 5.27 (FP32/FP16S) or 16.56 (FP32/FP16C) FLOPs/Byte. So performance is only limited by memory bandwidth.

If your GPU/CPU is not on the list yet, you can report your benchmarks [here](https://github.com/ProjectPhysX/FluidX3D/issues/8).

Colors: 🔴 AMD, 🔵 Intel, 🟢 Nvidia, 🟣 Apple, 🟡 Samsung

| Device                                           | FP32<br>[TFlops/s] | Mem<br>[GB] | BW<br>[GB/s] | FP32/FP32<br>[MLUPs/s] | FP32/FP16S<br>[MLUPs/s] | FP32/FP16C<br>[MLUPs/s] |
| :----------------------------------------------- | -----------------: | ----------: | -----------: | ---------------------: | ----------------------: | ----------------------: |
|                                                  |                    |             |              |                        |                         |                         |
| 🔴&nbsp;Instinct&nbsp;MI250&nbsp;(1&nbsp;GCD)    |              45.26 |          64 |         1638 |             5638 (53%) |              9030 (42%) |              8506 (40%) |
| 🔴&nbsp;Instinct&nbsp;MI100                      |              46.14 |          32 |         1228 |             5093 (63%) |              8133 (51%) |              8542 (54%) |
| 🔴&nbsp;Instinct&nbsp;MI60                       |              14.75 |          32 |         1024 |             3570 (53%) |              5047 (38%) |              5111 (38%) |
| 🔴&nbsp;Radeon&nbsp;VII                          |              13.83 |          16 |         1024 |             4898 (73%) |              7778 (58%) |              5256 (40%) |
| 🟢&nbsp;H100&nbsp;PCIe&nbsp;80GB                 |              51.01 |          80 |         2000 |       11128&nbsp;(85%) |             20624 (79%) |             13862 (53%) |
| 🟢&nbsp;A100&nbsp;SXM4&nbsp;80GB                 |              19.49 |          80 |         2039 |       10228&nbsp;(77%) |             18448 (70%) |             11197 (42%) |
| 🟢&nbsp;A100&nbsp;SXM4&nbsp;40GB                 |              19.49 |          40 |         1555 |             8522 (84%) |             16013 (79%) |             11251 (56%) |
| 🟢&nbsp;A100&nbsp;PCIe&nbsp;40GB                 |              19.49 |          40 |         1555 |             8526 (84%) |             16035 (79%) |             11088 (55%) |
| 🟢&nbsp;Tesla&nbsp;V100&nbsp;SXM2&nbsp;32GB      |              15.67 |          32 |          900 |             4471 (76%) |              8947 (77%) |              7217 (62%) |
| 🟢&nbsp;Tesla&nbsp;V100&nbsp;PCIe&nbsp;16GB      |              14.13 |          16 |          900 |             5128 (87%) |             10325 (88%) |              7683 (66%) |
| 🟢&nbsp;Quadro&nbsp;GV100                        |              16.66 |          32 |          870 |             3442 (61%) |              6641 (59%) |              5863 (52%) |
| 🟢&nbsp;Titan&nbsp;V                             |              14.90 |          12 |          653 |             3601 (84%) |              7253 (86%) |              6957 (82%) |
| 🟢&nbsp;Tesla&nbsp;P100&nbsp;16GB                |               9.52 |          16 |          732 |             3295 (69%) |              5950 (63%) |              4176 (44%) |
| 🟢&nbsp;Tesla&nbsp;P100&nbsp;12GB                |               9.52 |          12 |          549 |             2427 (68%) |              4141 (58%) |              3999 (56%) |
| 🟢&nbsp;GeForce&nbsp;GTX&nbsp;TITAN              |               4.71 |           6 |          288 |             1460 (77%) |              2500 (67%) |              1113 (30%) |
| 🟢&nbsp;Tesla&nbsp;K40m                          |               4.29 |          12 |          288 |             1131 (60%) |              1868 (50%) |               912 (24%) |
| 🟢&nbsp;Tesla&nbsp;K80&nbsp;(1&nbsp;GPU)         |               4.11 |          12 |          240 |              916 (58%) |              1642 (53%) |               943 (30%) |
| 🟢&nbsp;Tesla&nbsp;K20c                          |               3.52 |           5 |          208 |              861 (63%) |              1507 (56%) |               720 (27%) |
|                                                  |                    |             |              |                        |                         |                         |
| 🔴&nbsp;Radeon&nbsp;RX&nbsp;7900&nbsp;XTX        |              61.44 |          24 |          960 |             3665 (58%) |              7644 (61%) |              7716 (62%) |
| 🔴&nbsp;Radeon&nbsp;RX&nbsp;7900&nbsp;XT         |              51.61 |          20 |          800 |             3013 (58%) |              5856 (56%) |              5986 (58%) |
| 🔴&nbsp;Radeon&nbsp;RX&nbsp;7600                 |              21.75 |           8 |          288 |             1250 (66%) |              2561 (68%) |              2512 (67%) |
| 🔴&nbsp;Radeon&nbsp;RX&nbsp;6900&nbsp;XT         |              23.04 |          16 |          512 |             1968 (59%) |              4227 (64%) |              4207 (63%) |
| 🔴&nbsp;Radeon&nbsp;RX&nbsp;6800&nbsp;XT         |              20.74 |          16 |          512 |             2008 (60%) |              4241 (64%) |              4224 (64%) |
| 🔴&nbsp;Radeon&nbsp;Pro&nbsp;W6800               |              17.83 |          32 |          512 |             1620 (48%) |              3361 (51%) |              3180 (48%) |
| 🔴&nbsp;Radeon&nbsp;RX&nbsp;6700M                |              10.60 |          10 |          320 |             1194 (57%) |              2388 (57%) |              2429 (58%) |
| 🔴&nbsp;Radeon&nbsp;RX&nbsp;5700&nbsp;XT         |               9.75 |           8 |          448 |             1368 (47%) |              3253 (56%) |              3049 (52%) |
| 🔴&nbsp;Radeon&nbsp;RX&nbsp;Vega&nbsp;64         |              13.35 |           8 |          484 |             1875 (59%) |              2878 (46%) |              3227 (51%) |
| 🔴&nbsp;Radeon&nbsp;RX&nbsp;580&nbsp;4GB         |               6.50 |           4 |          256 |              946 (57%) |              1848 (56%) |              1577 (47%) |
| 🔴&nbsp;Radeon&nbsp;R9&nbsp;390X                 |               5.91 |           8 |          384 |             1733 (69%) |              2217 (44%) |              1722 (35%) |
| 🔴&nbsp;Radeon&nbsp;HD&nbsp;7850                 |               1.84 |           2 |          154 |              112 (11%) |               120 ( 6%) |               635 (32%) |
| 🔵&nbsp;Arc&nbsp;A770&nbsp;LE                    |              19.66 |          16 |          560 |             2741 (75%) |              4591 (63%) |              4626 (64%) |
| 🔵&nbsp;Arc&nbsp;A750&nbsp;LE                    |              17.20 |           8 |          512 |             2625 (78%) |              4184 (63%) |              4238 (64%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;4090               |              82.58 |          24 |         1008 |             5624 (85%) |             11091 (85%) |             11496 (88%) |
| 🟢&nbsp;RTX&nbsp;6000&nbsp;Ada                   |              91.10 |          48 |          960 |             4997 (80%) |             10249 (82%) |             10293 (83%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;4080               |              55.45 |          16 |          717 |             3914 (84%) |              7626 (82%) |              7933 (85%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;3090&nbsp;Ti       |              40.00 |          24 |         1008 |             5717 (87%) |             10956 (84%) |             10400 (79%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;3090               |              39.05 |          24 |          936 |             5418 (89%) |             10732 (88%) |             10215 (84%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;3080&nbsp;Ti       |              37.17 |          12 |          912 |             5202 (87%) |              9832 (87%) |              9347 (79%) |
| 🟢&nbsp;RTX&nbsp;A6000                           |              40.00 |          48 |          768 |             4421 (88%) |              8814 (88%) |              8533 (86%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;3080               |              29.77 |          10 |          760 |             4230 (85%) |              8118 (82%) |              7714 (78%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;3070               |              20.31 |           8 |          448 |             2578 (88%) |              5096 (88%) |              5060 (87%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;3060&nbsp;Ti       |              16.49 |           8 |          448 |             2644 (90%) |              5129 (88%) |              4718 (81%) |
| 🟢&nbsp;RTX&nbsp;A5000M                          |              16.59 |          16 |          448 |             2228 (76%) |              4461 (77%) |              3662 (63%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;3060               |              13.17 |          12 |          360 |             2108 (90%) |              4070 (87%) |              3566 (76%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;3060M              |              10.94 |           6 |          336 |             2019 (92%) |              4012 (92%) |              3572 (82%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;3050M              |               7.13 |           4 |          192 |             1180 (94%) |              2339 (94%) |              2016 (81%) |
| 🟢&nbsp;Titan&nbsp;RTX                           |              16.31 |          24 |          672 |             3471 (79%) |              7456 (85%) |              7554 (87%) |
| 🟢&nbsp;Quadro&nbsp;RTX&nbsp;6000                |              16.31 |          24 |          672 |             3307 (75%) |              6836 (78%) |              6879 (79%) |
| 🟢&nbsp;Quadro&nbsp;RTX&nbsp;8000&nbsp;Pass.     |              14.93 |          48 |          624 |             2591 (64%) |              5408 (67%) |              5607 (69%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;2080&nbsp;Ti       |              13.45 |          11 |          616 |             3194 (79%) |              6700 (84%) |              6853 (86%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;2080&nbsp;Sup.     |              11.34 |           8 |          496 |             2434 (75%) |              5284 (82%) |              5087 (79%) |
| 🟢&nbsp;Quadro&nbsp;RTX&nbsp;5000                |              11.15 |          16 |          448 |             2341 (80%) |              4766 (82%) |              4773 (82%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;2060&nbsp;Sup.     |               7.18 |           8 |          448 |             2503 (85%) |              5035 (87%) |              4463 (77%) |
| 🟢&nbsp;Quadro&nbsp;RTX&nbsp;4000                |               7.12 |           8 |          416 |             2284 (84%) |              4584 (85%) |              4062 (75%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;2060&nbsp;KO       |               6.74 |           6 |          336 |             1643 (75%) |              3376 (77%) |              3266 (75%) |
| 🟢&nbsp;GeForce&nbsp;RTX&nbsp;2060               |               6.74 |           6 |          336 |             1681 (77%) |              3604 (83%) |              3571 (82%) |
| 🟢&nbsp;GeForce&nbsp;GTX&nbsp;1660&nbsp;Sup.     |               5.03 |           6 |          336 |             1696 (77%) |              3551 (81%) |              3040 (70%) |
| 🟢&nbsp;Tesla&nbsp;T4                            |               8.14 |          15 |          300 |             1356 (69%) |              2869 (74%) |              2887 (74%) |
| 🟢&nbsp;GeForce&nbsp;GTX&nbsp;1660&nbsp;Ti       |               5.48 |           6 |          288 |             1467 (78%) |              3041 (81%) |              3019 (81%) |
| 🟢&nbsp;GeForce&nbsp;GTX&nbsp;1660               |               5.07 |           6 |          192 |             1016 (81%) |              1924 (77%) |              1992 (80%) |
| 🟢&nbsp;GeForce&nbsp;GTX&nbsp;1650M              |               3.20 |           4 |          128 |              706 (84%) |              1214 (73%) |              1400 (84%) |
| 🟢&nbsp;Titan&nbsp;Xp                            |              12.15 |          12 |          548 |             2919 (82%) |              5495 (77%) |              5375 (76%) |
| 🟢&nbsp;GeForce&nbsp;GTX&nbsp;1080&nbsp;Ti       |              12.06 |          11 |          484 |             2631 (83%) |              4837 (77%) |              4877 (78%) |
| 🟢&nbsp;GeForce&nbsp;GTX&nbsp;1080               |               9.78 |           8 |          320 |             1623 (78%) |              3100 (75%) |              3182 (77%) |
| 🟢&nbsp;GeForce&nbsp;GTX&nbsp;1060M              |               4.44 |           6 |          192 |              983 (78%) |              1882 (75%) |              1803 (72%) |
| 🟢&nbsp;GeForce&nbsp;GTX&nbsp;1050M Ti           |               2.49 |           4 |          112 |              631 (86%) |              1224 (84%) |              1115 (77%) |
| 🟢&nbsp;Quadro&nbsp;P1000                        |               1.89 |           4 |           82 |              426 (79%) |               839 (79%) |               778 (73%) |
| 🟢&nbsp;GeForce&nbsp;GTX&nbsp;970                |               4.17 |           4 |          224 |              980 (67%) |              1721 (59%) |              1623 (56%) |
| 🟢&nbsp;Quadro&nbsp;M4000                        |               2.57 |           8 |          192 |              899 (72%) |              1519 (61%) |              1050 (42%) |
| 🟢&nbsp;Tesla&nbsp;M60&nbsp;(1&nbsp;GPU)         |               4.82 |           8 |          160 |              853 (82%) |              1571 (76%) |              1557 (75%) |
| 🟢&nbsp;GeForce&nbsp;GTX&nbsp;960M               |               1.51 |           4 |           80 |              442 (84%) |               872 (84%) |               627 (60%) |
| 🟢&nbsp;GeForce&nbsp;GTX&nbsp;770                |               3.33 |           2 |          224 |              800 (55%) |              1215 (42%) |               876 (30%) |
| 🟢&nbsp;GeForce&nbsp;GTX&nbsp;680&nbsp;4GB       |               3.33 |           4 |          192 |              783 (62%) |              1274 (51%) |               814 (33%) |
| 🟢&nbsp;Quadro&nbsp;K2000                        |               0.73 |           2 |           64 |              312 (75%) |               444 (53%) |               171 (21%) |
| 🟢&nbsp;GeForce&nbsp;GT&nbsp;630&nbsp;(OEM)      |               0.46 |           2 |           29 |              151 (81%) |               185 (50%) |                78 (21%) |
| 🟢&nbsp;Quadro&nbsp;NVS&nbsp;290                 |               0.03 |       0.256 |            6 |                1 ( 2%) |                 1 ( 1%) |                 1 ( 1%) |
|                                                  |                    |             |              |                        |                         |                         |
| 🟣&nbsp;M2&nbsp;Max&nbsp;GPU&nbsp;38C&nbsp;32GB  |               9.73 |          22 |          400 |             2405 (92%) |              4641 (89%) |              2444 (47%) |
| 🟣&nbsp;M1&nbsp;Max&nbsp;GPU&nbsp;24C&nbsp;32GB  |               6.14 |          22 |          400 |             2369 (91%) |              4496 (87%) |              2777 (53%) |
| 🟣&nbsp;M1&nbsp;Pro&nbsp;GPU&nbsp;16C&nbsp;16GB  |               4.10 |          11 |          200 |             1204 (92%) |              2329 (90%) |              1855 (71%) |
| 🔴&nbsp;Radeon&nbsp;Vega&nbsp;8&nbsp;(4750G)     |               2.15 |          27 |           57 |              263 (71%) |               511 (70%) |               501 (68%) |
| 🔴&nbsp;Radeon&nbsp;Vega&nbsp;8&nbsp;(3500U)     |               1.23 |           7 |           38 |              157 (63%) |               282 (57%) |               288 (58%) |
| 🔵&nbsp;UHD&nbsp;Graphics&nbsp;Xe&nbsp;32EUs     |               0.74 |          25 |           51 |              128 (38%) |               245 (37%) |               216 (32%) |
| 🔵&nbsp;UHD&nbsp;Graphics&nbsp;630               |               0.46 |           7 |           51 |              151 (45%) |               301 (45%) |               187 (28%) |
| 🔵&nbsp;UHD&nbsp;Graphics&nbsp;P630              |               0.46 |          51 |           42 |              177 (65%) |               288 (53%) |               137 (25%) |
| 🔵&nbsp;HD&nbsp;Graphics&nbsp;5500               |               0.35 |           3 |           26 |               75 (45%) |               192 (58%) |               108 (32%) |
| 🔵&nbsp;HD&nbsp;Graphics&nbsp;4600               |               0.38 |           2 |           26 |              105 (63%) |               115 (35%) |                34 (10%) |
| 🟡&nbsp;ARM&nbsp;Mali-G72&nbsp;MP18              |               0.24 |           4 |           29 |               14 ( 7%) |                17 ( 5%) |                12 ( 3%) |
|                                                  |                    |             |              |                        |                         |                         |
| 🔴&nbsp;2x&nbsp;EPYC&nbsp;9654                   |              29.49 |        1536 |          922 |             1381 (23%) |              1814 (15%) |              1801 (15%) |
| 🔵&nbsp;2x&nbsp;Xeon&nbsp;CPU&nbsp;Max&nbsp;9480 |              13.62 |         256 |          614 |             2037 (51%) |              1520 (19%) |              1464 (18%) |
| 🔵&nbsp;2x&nbsp;Xeon&nbsp;Platinum&nbsp;8480+    |              14.34 |         512 |          614 |             2162 (54%) |              1845 (23%) |              1884 (24%) |
| 🔵&nbsp;2x&nbsp;Xeon&nbsp;Platinum&nbsp;8380     |              11.78 |        2048 |          410 |             1410 (53%) |              1159 (22%) |              1298 (24%) |
| 🔵&nbsp;2x&nbsp;Xeon&nbsp;Platinum&nbsp;8358     |              10.65 |         256 |          410 |             1285 (48%) |              1007 (19%) |              1120 (21%) |
| 🔵&nbsp;2x&nbsp;Xeon&nbsp;Platinum&nbsp;8256     |               1.95 |        1536 |          282 |              396 (22%) |               158 ( 4%) |               175 ( 5%) |
| 🔵&nbsp;2x&nbsp;Xeon&nbsp;Platinum&nbsp;8153     |               4.10 |         384 |          256 |              691 (41%) |               290 ( 9%) |               328 (10%) |
| 🔵&nbsp;2x&nbsp;Xeon&nbsp;Gold&nbsp;6128         |               2.61 |         192 |          256 |              254 (15%) |               185 ( 6%) |               193 ( 6%) |
| 🔵&nbsp;Xeon&nbsp;Phi&nbsp;7210                  |               5.32 |         192 |          102 |              415 (62%) |               193 (15%) |               223 (17%) |
| 🔵&nbsp;4x&nbsp;Xeon&nbsp;E5-4620&nbsp;v4        |               2.69 |         512 |          273 |              460 (26%) |               275 ( 8%) |               239 ( 7%) |
| 🔵&nbsp;2x&nbsp;Xeon&nbsp;E5-2630&nbsp;v4        |               1.41 |          64 |          137 |              264 (30%) |               146 ( 8%) |               129 ( 7%) |
| 🔵&nbsp;2x&nbsp;Xeon&nbsp;E5-2623&nbsp;v4        |               0.67 |          64 |          137 |              125 (14%) |                66 ( 4%) |                59 ( 3%) |
| 🔵&nbsp;2x&nbsp;Xeon&nbsp;E5-2680&nbsp;v3        |               1.92 |          64 |          137 |              209 (23%) |               305 (17%) |               281 (16%) |
| 🔵&nbsp;Core&nbsp;i9-11900KB                     |               0.84 |          32 |           51 |              109 (33%) |               195 (29%) |               208 (31%) |
| 🔵&nbsp;Core&nbsp;i9-10980XE                     |               3.23 |         128 |           94 |              286 (47%) |               251 (21%) |               223 (18%) |
| 🔵&nbsp;Core&nbsp;i5-9600                        |               0.60 |          16 |           43 |              146 (52%) |               127 (23%) |               147 (27%) |
| 🔵&nbsp;Core&nbsp;i7-8700K                       |               0.71 |          16 |           51 |              152 (45%) |               134 (20%) |               116 (17%) |
| 🔵&nbsp;Xeon&nbsp;E-2176G                        |               0.71 |          64 |           42 |              201 (74%) |               136 (25%) |               148 (27%) |
| 🔵&nbsp;Core&nbsp;i7-7700HQ                      |               0.36 |          12 |           38 |               81 (32%) |                82 (16%) |               108 (22%) |
| 🔵&nbsp;Core&nbsp;i7-4770                        |               0.44 |          16 |           26 |              104 (62%) |                69 (21%) |                59 (18%) |
| 🔵&nbsp;Core&nbsp;i7-4720HQ                      |               0.33 |          16 |           26 |               58 (35%) |                13 ( 4%) |                47 (14%) |



## Multi-GPU Benchmarks

Multi-GPU benchmarks are done at the largest possible grid resolution with a cubic domain, and either 2x1x1, 2x2x1 or 2x2x2 of these cubic domains together. The percentages in brackets are single-GPU roofline model efficiency, and the multiplicator numbers in brackets are scaling factors relative to benchmarked single-GPU performance.

Colors: 🔴 AMD, 🔵 Intel, 🟢 Nvidia, 🟣 Apple, 🟡 Samsung

| Device                                                          | FP32<br>[TFlops/s] | Mem<br>[GB] | BW<br>[GB/s] | FP32/FP32<br>[MLUPs/s] | FP32/FP16S<br>[MLUPs/s] | FP32/FP16C<br>[MLUPs/s] |
| :-------------------------------------------------------------- | -----------------: | ----------: | -----------: | ---------------------: | ----------------------: | ----------------------: |
|                                                                 |                    |             |              |                        |                         |                         |
| 🔴&nbsp;1x&nbsp;Instinct&nbsp;MI250&nbsp;(1&nbsp;GCD)           |              45.26 |          64 |         1638 |             5638 (53%) |              9030 (42%) |              8506 (40%) |
| 🔴&nbsp;1x&nbsp;Instinct&nbsp;MI250&nbsp;(2&nbsp;GCD)           |              90.52 |         128 |         3277 |            9460 (1.7x) |            14313 (1.6x) |            17338 (2.0x) |
| 🔴&nbsp;2x&nbsp;Instinct&nbsp;MI250&nbsp;(4&nbsp;GCD)           |             181.04 |         256 |         6554 |      16925&nbsp;(3.0x) |            29163 (3.2x) |            29627 (3.5x) |
| 🔴&nbsp;4x&nbsp;Instinct&nbsp;MI250&nbsp;(8&nbsp;GCD)           |             362.08 |         512 |        13107 |      27350&nbsp;(4.9x) |            52258 (5.8x) |            53521 (6.3x) |
|                                                                 |                    |             |              |                        |                         |                         |
| 🔴&nbsp;1x&nbsp;Radeon&nbsp;VII                                 |              13.83 |          16 |         1024 |             4898 (73%) |              7778 (58%) |              5256 (40%) |
| 🔴&nbsp;2x&nbsp;Radeon&nbsp;VII                                 |              27.66 |          32 |         2048 |            8113 (1.7x) |            15591 (2.0x) |            10352 (2.0x) |
| 🔴&nbsp;4x&nbsp;Radeon&nbsp;VII                                 |              55.32 |          64 |         4096 |      12911&nbsp;(2.6x) |            24273 (3.1x) |            17080 (3.2x) |
| 🔴&nbsp;8x&nbsp;Radeon&nbsp;VII                                 |             110.64 |         128 |         8192 |      21946&nbsp;(4.5x) |            30826 (4.0x) |            24572 (4.7x) |
|                                                                 |                    |             |              |                        |                         |                         |
| 🟢&nbsp;1x&nbsp;A100&nbsp;SXM4&nbsp;40GB                        |              19.49 |          40 |         1555 |             8543 (84%) |             15917 (79%) |              8748 (43%) |
| 🟢&nbsp;2x&nbsp;A100&nbsp;SXM4&nbsp;40GB                        |              38.98 |          80 |         3110 |      14311&nbsp;(1.7x) |            23707 (1.5x) |            15512 (1.8x) |
| 🟢&nbsp;4x&nbsp;A100&nbsp;SXM4&nbsp;40GB                        |              77.96 |         160 |         6220 |      23411&nbsp;(2.7x) |            42400 (2.7x) |            29017 (3.3x) |
| 🟢&nbsp;8x&nbsp;A100&nbsp;SXM4&nbsp;40GB                        |             155.92 |         320 |        12440 |      37619&nbsp;(4.4x) |            72965 (4.6x) |            63009 (7.2x) |
|                                                                 |                    |             |              |                        |                         |                         |
| 🟢&nbsp;1x&nbsp;A100&nbsp;SXM4&nbsp;40GB                        |              19.49 |          40 |         1555 |             8522 (84%) |             16013 (79%) |             11251 (56%) |
| 🟢&nbsp;2x&nbsp;A100&nbsp;SXM4&nbsp;40GB                        |              38.98 |          80 |         3110 |      13629&nbsp;(1.6x) |            24620 (1.5x) |            18850 (1.7x) |
| 🟢&nbsp;4x&nbsp;A100&nbsp;SXM4&nbsp;40GB                        |              77.96 |         160 |         6220 |      17978&nbsp;(2.1x) |            30604 (1.9x) |            30627 (2.7x) |
|                                                                 |                    |             |              |                        |                         |                         |
| 🟢&nbsp;1x&nbsp;Tesla&nbsp;V100&nbsp;SXM2&nbsp;32GB             |              15.67 |          32 |          900 |             4471 (76%) |              8947 (77%) |              7217 (62%) |
| 🟢&nbsp;2x&nbsp;Tesla&nbsp;V100&nbsp;SXM2&nbsp;32GB             |              31.34 |          64 |         1800 |            7953 (1.8x) |            15469 (1.7x) |            12932 (1.8x) |
| 🟢&nbsp;4x&nbsp;Tesla&nbsp;V100&nbsp;SXM2&nbsp;32GB             |              62.68 |         128 |         3600 |      13135&nbsp;(2.9x) |            26527 (3.0x) |            22686 (3.1x) |
|                                                                 |                    |             |              |                        |                         |                         |
| 🟢&nbsp;1x&nbsp;Tesla&nbsp;K40m                                 |               4.29 |          12 |          288 |             1131 (60%) |              1868 (50%) |               912 (24%) |
| 🟢&nbsp;2x&nbsp;Tesla&nbsp;K40m                                 |               8.58 |          24 |          577 |            1971 (1.7x) |             3300 (1.8x) |             1801 (2.0x) |
| 🟢&nbsp;3x&nbsp;K40m&nbsp;+&nbsp;1x&nbsp;Titan&nbsp;Xp          |              17.16 |          48 |         1154 |            3117 (2.8x) |             5174 (2.8x) |             3127 (3.4x) |
|                                                                 |                    |             |              |                        |                         |                         |
| 🟢&nbsp;1x&nbsp;RTX&nbsp;A6000                                  |              40.00 |          48 |          768 |             4421 (88%) |              8814 (88%) |              8533 (86%) |
| 🟢&nbsp;2x&nbsp;RTX&nbsp;A6000                                  |              80.00 |          96 |         1536 |            8041 (1.8x) |            15026 (1.7x) |            14795 (1.7x) |
| 🟢&nbsp;4x&nbsp;RTX&nbsp;A6000                                  |             160.00 |         192 |         3072 |      14314&nbsp;(3.2x) |            27915 (3.2x) |            27227 (3.2x) |
| 🟢&nbsp;8x&nbsp;RTX&nbsp;A6000                                  |             320.00 |         384 |         6144 |      19311&nbsp;(4.4x) |            40063 (4.5x) |            39004 (4.6x) |
|                                                                 |                    |             |              |                        |                         |                         |
| 🟢&nbsp;1x&nbsp;Quadro&nbsp;RTX&nbsp;8000&nbsp;Pa.              |              14.93 |          48 |          624 |             2591 (64%) |              5408 (67%) |              5607 (69%) |
| 🟢&nbsp;2x&nbsp;Quadro&nbsp;RTX&nbsp;8000&nbsp;Pa.              |              29.86 |          96 |         1248 |            4767 (1.8x) |             9607 (1.8x) |            10214 (1.8x) |
|                                                                 |                    |             |              |                        |                         |                         |
| 🟢&nbsp;1x&nbsp;GeForce&nbsp;RTX&nbsp;2080&nbsp;Ti              |              13.45 |          11 |          616 |             3194 (79%) |              6700 (84%) |              6853 (86%) |
| 🟢&nbsp;2x&nbsp;GeForce&nbsp;RTX&nbsp;2080&nbsp;Ti              |              26.90 |          22 |         1232 |            5085 (1.6x) |            10770 (1.6x) |            10922 (1.6x) |
| 🟢&nbsp;4x&nbsp;GeForce&nbsp;RTX&nbsp;2080&nbsp;Ti              |              53.80 |          44 |         2464 |            9117 (2.9x) |            18415 (2.7x) |            18598 (2.7x) |
| 🟢&nbsp;7x&nbsp;2080&nbsp;Ti&nbsp;+&nbsp;1x&nbsp;A100&nbsp;40GB |             107.60 |          88 |         4928 |      16146&nbsp;(5.1x) |            33732 (5.0x) |            33857 (4.9x) |



## FAQs

### General

- <details><summary>What physical model does FluidX3D use?</summary><br>FluidX3D implements the lattice Boltzmann method, a type of direct numerical simulation (DNS), the most accurate type of fluid simulation, but also the most computationally challenging. Optional extension models include volume force (Guo forcing), free surface (<a href="https://doi.org/10.3390/computation10060092">volume-of-fluid</a> and <a href="https://doi.org/10.3390/computation10020021">PLIC</a>), a temperature model and Smagorinsky-Lilly subgrid turbulence model.<br><br></details>

- <details><summary>FluidX3D only uses FP32 or even FP32/FP16, in contrast to FP64. Are simulation results physically accurate?</summary><br>Yes, in all but extreme edge cases. The code has been specially optimized to minimize arithmetic round-off errors and make the most out of lower precision. With these optimizations, accuracy in most cases is indistinguishable from FP64 double-precision, even with FP32/FP16 mixed-precision. Details can be found in <a href="https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats">this paper</a>.<br><br></details>

- <details><summary>Why is the domain size limited to 2³² grid points?</summary><br>The 32-bit unsigned integer grid index will overflow above this number. Using 64-bit index calculation would slow the simulation down by ~20%, as 64-bit uint is calculated on special function units and not the regular GPU cores. 2³² grid points with FP32/FP16 mixed-precision is equivalent to 225GB memory and single GPUs currently are only at 128GB, so it should be fine for a while to come. For higher resolutions above the single-domain limit, use multiple domains (typically 1 per GPU, but multiple domains on the same GPU also work).<br><br></details>

- <details><summary>Comparted to the benchmark numbers stated <a href="https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats">here</a>, efficiency seems much lower but performance is slightly better for most devices. How can this be?</summary><br>In that paper, the One-Step-Pull swap algorithm is implemented, using only misaligned reads and coalesced writes. On almost all GPUs, the performance penalty for misaligned writes is much larger than for misaligned reads, and sometimes there is almost no penalty for misaligned reads at all. Because of this, One-Step-Pull runs at peak bandwidth and thus peak efficiency.<br>Here, a different swap algorithm termed <a href="https://doi.org/10.3390/computation10060092">Esoteric-Pull</a> is used, a type of in-place streaming. This makes the LBM require much less memory (93 vs. 169 (FP32/FP32) or 55 vs. 93 (FP32/FP16) Bytes/cell for D3Q19), and also less memory bandwidth (153 vs. 171 (FP32/FP32) or 77 vs. 95 (FP32/FP16) Bytes/cell per time step for D3Q19) due to so-called implicit bounce-back boundaries. However memory access now is half coalesced and half misaligned for both reads and writes, so memory access efficiency is lower. For overall performance, these two effects approximately cancel out. The benefit of Esoteric-Pull - being able to simulate domains twice as large with the same amount of memory - clearly outweights the cost of slightly lower memory access efficiency, especially since performance is not reduced overall.<br><br></details>

- <details><summary>Why don't you use CUDA? Wouldn't that be more efficient?</summary><br>No, that is a wrong myth. OpenCL is exactly as efficient as CUDA on Nvidia GPUs if optimized properly. <a href="https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats">Here</a> I did roofline model and analyzed OpenCL performance on various hardware. OpenCL efficiency on modern Nvidia GPUs can be 100% with the right memory access pattern, so CUDA can't possibly be any more efficient. Without any performance advantage, there is no reason to use proprietary CUDA over OpenCL, since OpenCL is compatible with a lot more hardware.<br><br></details>

- <details><summary>Why no multi-relaxation-time (MRT) collision operator?</summary><br>The idea of MRT is to linearly transform the DDFs into "moment space" by matrix multiplication and relax these moments individually, promising better stability and accuracy. In practice, in the vast majority of cases, it has zero or even negative effects on stability and accuracy, and simple SRT is much superior. Apart from the kinematic shear viscosity and conserved terms, the remaining moments are non-physical quantities and their tuning is a blackbox. Although MRT can be implemented in an efficient manner with only a single matrix-vector multiplication in registers, leading to identical performance compared to SRT by remaining bandwidth-bound, storing the matrices vastly elongates and over-complicates the code for no real benefit.<br><br></details>

### Hardware

- <details><summary>Can FluidX3D run on multiple GPUs at the same time?</summary><br>Yes. The simulation grid is then split in domains, one for each GPU (domain decomposition method). The GPUs essentially pool their memory, enabling much larger grid resolution and higher performance. Rendering is parallelized across multiple GPUs as well; each GPU renders its own domain with a 3D offset, then rendered frames from all GPUs are overlayed with their z-buffers. Communication between domains is done over PCIe, so no SLI/Crossfire/NVLink/InfinityFabric is required. All GPUs must however be installed in the same node (PC/laptop/server). Even unholy combinations of Nvidia/AMD/Intel GPUs will work, although it is recommended to only use GPUs with similar memory capacity and bandwidth together. Using a fast gaming GPU and slow integrated GPU together would only decrease performance due to communication overhead.<br><br></details>

- <details><summary>I'm on a budget and have only a cheap computer. Can I run FluidX3D on my toaster PC/laptop?</summary><br>Absolutely. Today even the most inexpensive hardware, like integrated GPUs or entry-level gaming GPUs, support OpenCL. You might be a bit more limited on memory capacity and grid resolution, but you should be good to go. I've tested FluidX3D on very old and inexpensive hardware and even on my Samsung S9+ smartphone, and it runs just fine, although admittedly a bit slower.<br><br></details>

- <details><summary>I don't have an expensive workstation GPU, but only a gaming GPU. Will performance suffer?</summary><br>No. Efficiency on gaming GPUs is exactly as good as on their "professional"/workstation counterparts. Performance often is even better as gaming GPUs have higher boost clocks.<br><br></details>

- <details><summary>Do I need a GPU with ECC memory?</summary><br>No. Gaming GPUs work just fine. Some Nvidia GPUs automatically reduce memory clocks for compute applications to almost entirely eliminate memory errors.<br><br></details>

- <details><summary>My GPU does not support CUDA. Can I still use FluidX3D?</summary><br>Yes. FluidX3D uses OpenCL 1.2 and not CUDA, so it runs on any GPU from any vendor since around 2012.<br><br></details>

- <details><summary>I don't have a dedicated graphics card at all. Can I still run FluidX3D on my PC/laptop?</summary><br>Yes. FluidX3D also runs on all integrated GPUs since around 2012, and also on CPUs.<br><br></details>

- <details><summary>I need more memory than my GPU can offer. Can I run FluidX3D on my CPU as well?</summary><br>Yes. You only need to install the <a href="https://github.com/intel/llvm/releases/tag/2022-09">OpenCL Runtime for Intel CPUs</a>.<br><br></details>

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

- Lehmann, M.: [Computational study of microplastic transport at the water-air interface with a memory-optimized lattice Boltzmann method](https://doi.org/10.15495/EPub_UBT_00006977). PhD thesis, (2023)
- Lehmann, M.: [Esoteric Pull and Esoteric Push: Two Simple In-Place Streaming Schemes for the Lattice Boltzmann Method on GPUs](https://doi.org/10.3390/computation10060092). Computation, 10, 92, (2022)
- Lehmann, M., Krause, M., Amati, G., Sega, M., Harting, J. and Gekle, S.: [Accuracy and performance of the lattice Boltzmann method with 64-bit, 32-bit, and customized 16-bit number formats](https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats). Phys. Rev. E 106, 015308, (2022)
- Lehmann, M.: [Combined scientific CFD simulation and interactive raytracing with OpenCL](https://www.researchgate.net/publication/360501260_Combined_scientific_CFD_simulation_and_interactive_raytracing_with_OpenCL). IWOCL'22: International Workshop on OpenCL, 3, 1-2, (2022)
- Lehmann, M., Oehlschlägel, L.M., Häusl, F., Held, A. and Gekle, S.: [Ejection of marine microplastics by raindrops: a computational and experimental study](https://doi.org/10.1186/s43591-021-00018-8). Micropl.&Nanopl. 1, 18, (2021)
- Lehmann, M.: [High Performance Free Surface LBM on GPUs](https://doi.org/10.15495/EPub_UBT_00005400). Master's thesis, (2019)
- Lehmann, M. and Gekle, S.: [Analytic Solution to the Piecewise Linear Interface Construction Problem and Its Application in Curvature Calculation for Volume-of-Fluid Simulation Codes](https://doi.org/10.3390/computation10020021). Computation, 10, 21, (2022)



## Contact

- FluidX3D is solo-developed and maintained by Dr. Moritz Lehmann.
- For any questions, feedback or other inquiries, contact me at [moritz.lehmann@uni-bayreuth.de](mailto:moritz.lehmann@uni-bayreuth.de?subject=FluidX3D).
- Updates will be posted on Twitter via [@FluidX3D](https://twitter.com/FluidX3D) and [@ProjectPhysX](https://twitter.com/ProjectPhysX), under the hashtag [#FluidX3D](https://twitter.com/hashtag/FluidX3D?src=hashtag_click&f=live) or on my [YouTube channel](https://www.youtube.com/c/ProjectPhysX).