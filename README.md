# Capsid-Modeling
This code solves Surface harmonics to model Virus Capsids. The equation for surface harmonics is given as:\
Yₗᵐ(θ, φ) = √[(2ℓ + 1)/(4π) · (ℓ - m)! / (ℓ + m)!] · Pₗᵐ(cos θ) · e^(i m φ) \
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
