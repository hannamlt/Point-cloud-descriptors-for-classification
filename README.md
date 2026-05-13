# 3D Point Cloud Descriptors for Semantic Segmentation

Research internship · Inria Sophia Antipolis 
Contribution to the [CGAL open-source library](https://www.cgal.org/)

## Overview

Novel topological descriptors for semantic segmentation of 3D point clouds, 
implemented in C++ within the CGAL framework. Two complementary approaches 
were developed and integrated as open-source contributions to CGAL.


## Methods

**3D Compactness Descriptor**  
Extension of 2D compactness theory to 3D volumes. Quantifies spherical 
similarity using volume/surface area ratio, computed via CGAL alpha-wrapping. 
Scale and rotation invariant. Multi-scale analysis across multiple neighborhood radii.

**Skeletonization Pipeline**  
Mean curvature flow (MCF) skeletonization for medial axis extraction. 
Distance-to-skeleton features capture topological hierarchy of branching structures 
(e.g. tree trunks vs. foliage). Implemented using CGAL surface mesh skeletonization.


## Results

| Metric | Baseline | + Compactness descriptor |
|--------|----------|--------------------------|
| Accuracy | 98.58% | 99.57% (+1.0%) |
| F1-score | 0.9888 | 0.9946 (+0.58%) |

**Performance optimizations: 30%+ speedup** (5.2h → 3.5h) via k-d trees, 
planimetric grids, and OpenMP parallelization.


## Structure


---

## Stack

`C++` `CGAL` `CMake` `OpenMP` `Random Forest`

## Related

- [CGAL Alpha Wrapping](https://doc.cgal.org/latest/Alpha_wrap_3/)
- [CGAL Point Set Classification](https://doc.cgal.org/latest/Classification/)
- [CGAL Surface Mesh Skeletonization](https://doc.cgal.org/latest/Surface_mesh_skeletonization/)
