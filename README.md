# ACEMD

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ACEsuit.github.io/ACEMD.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ACEsuit.github.io/ACEMD.jl/dev/)
[![Build Status](https://github.com/ACEsuit/ACEMD.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ACEsuit/ACEMD.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a test package to make API for [ACE](https://github.com/ACEsuit) molecular dynamics that can be used by other packages.

## Install

Install ACE repo, if you haven't done so yet

```julia
] registry add https://github.com/ACEsuit/ACEregistry.git
```

After that

```julia
using Pkg
Pkg.add("ACEMD")
```

Note that, if you want to do ACE calculations in Molly, you need to add both Molly and ACEMD to **the same environment**. Otherwise the package extension does not load.
