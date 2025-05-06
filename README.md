# Capsid-Modeling
This code solves Surface harmonics to model Virus Capsids. The equation for surface harmonics is given as:\
![Spherical Harmonics Equation](https://latex.codecogs.com/png.image?\dpi{150}Y_\ell^m(\theta,\phi)=\sqrt{\frac{(2\ell+1)}{4\pi}\cdot\frac{(\ell-m)!}{(\ell+m)!}}P_\ell^m(\cos\theta)e^{im\phi}) \
Inline code comments can use visual studio extension to render equations.

# Building

git clone https://github.com/abashsharma/Capsid-Modeling.git \
mkdir build && cd build \
cmake .. \
cmake --build . --config Release

# Code
Written in C++23 for mdspan etc
(No C++ compiler on BioHPC, so needs to be run in a conda env)

To dos: Search for Todo in the project to find recommendations.
