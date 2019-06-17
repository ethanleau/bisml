Adaptive BRDF-Oriented Multiple Importance Sampling of Many Lights
===============

## Overview

This project contains the source code for the EGSR 2019 paper:

**Adaptive BRDF-Oriented Multiple Importance Sampling of Many Lights**, by

Yifan Liu, Kun Xu and Ling-Qi Yan.

The implementation is based on [pbrt-v3](https://github.com/mmp/pbrt-v3) and we mainly modified/added:

* `src/integrators/path.h(cpp)`: integrator of path tracing
* `src/acceleratiors/lightree.h(cpp)`: light tree
* `src/accelerators/lighttreesampler.h(cpp)`: direct illumination sampling techniques, core idea of our paper
* `src/materials/ltc.h(cpp)`: evaluation of brdf integral



Usage:
--------------

Usage should be relatively clear for developers familiar with PBRT.

```bash
$ git clone --recursive https://github.com/lyfxyz/blsml/
$ cd blsml
$ mkdir build
$ cd build
$ make ..
```
Please refer to  [pbrt-v3](https://github.com/mmp/pbrt-v3) for more details on building.



## Test Scene

In `scene/` folder, we contain an example scene: Cornell Box with thousands of small lights, which is used in our paper. To render this scene, just:

```bash
$ ./pbrt PATH2SCENE/cornell-box/scene.pbrt
```

