# Rain Simulation

This project is a C++ tool to calculate the optimal velocity at which bodies should move through rain to minimize wetness. The program able to estimate how wet a body gets when it travels a straight path under rain at a fixed velocity. It finds a velocity near the minimum using Brent's method and performs a parabolic fit around this point to better estimate the minimum. Bodies are defined in XML files, and simulations are configured through an `Input.xml` file. 

## Installation

Requirements:
- `clang++` (C++17 or newer)
- `make`

From the `Code` directory build with:
```bash
make
```
Other useful commands:
```bash
make debug   # build with debugging flags
make format  # format source files with clang-format
make clean   # remove object files and executable
```

The compiled binary is `main.exe`.

## Usage
The program requires:
- An `Input.xml` file in the same directory as `main.exe` and the Makefile, which specifies the simulation settings (see [Input.xml](Code/Input.xml) for an example).  
- A body XML file which describes the body shape and dynamics (see the [Bodies](Bodies/) folder for examples).

Run with:
```bash
make run
```
Or directly:
```bash
./main.exe   # Linux/macOS
.\main.exe   # Windows
```

## Body Models
<table>
  <tr>
    <td align="center" style="padding-left: 40px;">
      <img src="Media/Walk.gif" height="250px"><br>
    </td>
    <td align="center"style="padding-left: 40px;">
      <img src="Media/Run.gif" height="250px"><br>
    </td>
  </tr>
</table>

Renderings of the bodies defined in [WalkingMan.xml](Bodies/WalkingMan.xml) and [RunningMan.xml](Bodies/RunningMan.xml) respectively.

## Output
Below is an example output:
```csv
############################################################################################################
# Minimums of wetness found for the following parameters:
# Body = RunningMan.xml, vmax = 2 vfall
# dx = 0.01 m, nstep = 50
# nfit = 9, dv = 0.006 vfall
# Columns:
# vcross (vfall)   vtail (vfall)   vopt (vfall)   vopt_std (vfall)   Rmin (m^2)   Rmin_std (m^2)
############################################################################################################
0.000000000000	0.000000000000	2.000000000000	0.000000000000	0.438482950029	0.000000000000
0.000000000000	0.500000000000	0.874593172486	0.000800133595	0.268469915778	0.000006118143
0.000000000000	1.000000000000	1.201341426966	0.000861805107	0.163199562816	0.000020078286
0.500000000000	0.000000000000	2.000000000000	0.000000000000	0.458102375515	0.000000000000
0.500000000000	0.500000000000	1.123105402195	0.001327244707	0.333991782935	0.000008355485
0.500000000000	1.000000000000	1.267568725214	0.000436180555	0.216661264164	0.000005192611
```
Running the program without changing `Input.xml` should yeld this same output. Note that all velocities are expressed in units of the vertical velocity of the rain `vfall`, and that `Rmin` represents the minimum wetness achieved at the optimal velocity `vopt`.