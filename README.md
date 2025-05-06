# Spherical Harmonics
Spherical Harmonics are special functions defined on the surface of a sphere that arise in solving partial differential equations in spherical coordinates, especially in physics and engineering. They are denoted as _Yₗᵐ(θ, φ)_ where _l_ is the degree and _m_ is the order. These functions form an orthonormal basis for square-integrable functions on the sphere, making them fundamental in applications such as quantum mechanics, geophysics, computer graphics, and the analysis of spherical data. Their angular structure captures both symmetry and variation on a spherical domain.


# Capsid-Modeling
This code solves Surface harmonics to model Virus Capsids. The equation for surface harmonics is given as:\
![Spherical Harmonics Equation](https://latex.codecogs.com/png.image?\dpi{150}Y_\ell^m(\theta,\phi)=\sqrt{\frac{(2\ell+1)}{4\pi}\cdot\frac{(\ell-m)!}{(\ell+m)!}}P_\ell^m(\cos\theta)e^{im\phi}) \
Inline code comments can use visual studio extension to render equations.

# Usage

The minimum dependencies for executing the code are:\
cmake VERSION 3.16 \
C++ Version 23

# How to compile

git clone https://github.com/abashsharma/Capsid-Modeling.git \
mkdir build && cd build \
cmake .. \
cmake --build . --config Release \

Execute the code inside the built 'Capsid' directory as: \
./Phase \
In order to run a specific minimization routine for the simulation, you can pass it as an argument:\
./Phase -MC \

A proper test case is sumarized as: \
./Phase -MC -N 5 -q 100 \
where N is the desired number of Harmonic modes, and q is the desired number of grid points on the geometry. \

Other parameters such as the harmonic coefficients can be changed inside the code, and should be built again for now, while we explore the necessary parameters of interest for the study. 

# Output

The simulation outputs several files including the following: \
Initial.xyz (Initial Configuration) \
radius.txt (Average radius of the geometry as the geometry changes) \
E.txt (Total Energy of the system, to check if the system has achieved equilibrium) \
mc_...txt (the xyz cooardinates of the evolution of geometry)

# Visualization

The geometry can be visualized using softwares that render xyz coardinates. One suggested example is using vmd. \
vmd mc_..txt

