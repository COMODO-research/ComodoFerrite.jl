<img src="https://github.com/COMODO-research/ComodoFerrite/blob/main/assets/img/ComodoFerrite_logo.jpg" alt="ComodoFerrite uniaxial loading of a cube" width="50%"/>

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache-blue.svg)](https://github.com/COMODO-research/ComodoFerrite.jl/blob/main/LICENSE)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://COMODO-research.github.io/ComodoFerrite.jl/dev/)
[![Build Status](https://github.com/COMODO-research/ComodoFerrite.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/COMODO-research/ComodoFerrite.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/COMODO-research/ComodoFerrite.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/COMODO-research/ComodoFerrite.jl)

**A Julia package to combine the powers of [Comodo.jl](https://github.com/COMODO-research/Comodo.jl) and [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl)**

# About ComodoFerrite.jl
This Julia package enables one to combine [Comodo.jl](https://github.com/COMODO-research/Comodo.jl) and [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl). The latter is a powerful Julia library for finite element analysis, while the former offers powerful geometry processing and meshing tools. In addition, Comodo links with Makie for advanced Julia based visualisation. This library helps brings these capabilities together. Functionality provided here includes the conversion between Comodo meshes to Ferrite meshes and vice verse. 

<img src="https://github.com/COMODO-research/ComodoFerrite/blob/main/assets/anim/ComodoFerrite_cube_uniaxial_linear_elasticity.gif" alt="ComodoFerrite uniaxial loading of a cube" width="50%"/>

# Installation
To install
```
julia> ]

(@v1.11) pkg> add https://github.com/COMODO-research/ComodoFerrite.jl

```
<!--[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://COMODO-research.github.io/ComodoFerrite.jl/stable/)-->

# Documentation
Under construction. For now please see the examples in the [`examples` folder](https://github.com/COMODO-research/ComodoFerrite.jl/tree/main/examples). 

# License 
ComodoFerrite.jl is released open source under the [Apache 2.0 license](https://github.com/COMODO-research/ComodoFerrite.jl/blob/main/LICENSE).
