### What’s new in v2.0.0
1. Manifold initialization is now handled directly within the GrowFundCurv1D function — it no longer requires an initial segment computed externally.
2. The algorithm now computes an initial fundamental domain using the linear approximation of the manifold. Users can specify the desired distance of this domain from the fixed point.
3. Users can define both the approximate arclength of the manifold to compute and a maximum number of fundamental domain iterations. The algorithm will grow the manifold up to the specified arclength unless it reaches the iteration limit first.
4. For orientation-reversing manifolds (with negative eigenvalues), the algorithm now computes both branches using the first iterate of the map. In v1.0.0, it used the second iterate and handled each branch separately.
5. A new option allows computation of a single orbit on a two-dimensional manifold using one of the eigenvectors from its linear approximation.
6. The algorithm now stores the indices corresponding to each fundamental domain.
7. The inter_planes function no longer adds intersection points directly to the manifold. Instead, it records each intersection point along with the index of the preceding point in the manifold that intersected the plane.
