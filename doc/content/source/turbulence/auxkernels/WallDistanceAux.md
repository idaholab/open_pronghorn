# WallDistanceAux

`WallDistanceAux` is an auxiliary finite volume (FV) kernel that **computes the distance from each
quadrature point to the nearest solid wall boundary**. The resulting distance field can be used
by near-wall turbulence models, wall functions, or any model that requires a wall-normal length
scale (e.g., for computing $y^+$ or two-layer blending functions).

More precisely, for each quadrature point located at position $\mathbf{x}_p$ inside an element,
`WallDistanceAux` computes

\begin{equation}
d(\mathbf{x}_p) = \min_{\mathbf{x}_w \in \Gamma_{\text{walls}}} \lVert \mathbf{x}_p - \mathbf{x}_w \rVert,
\end{equation}

where $\Gamma_{\text{walls}}$ is the union of all boundary faces designated as solid walls, and
$\mathbf{x}_w$ is taken as the centroid of a boundary face in the mesh.

The computed value is stored in the auxiliary variable associated with this kernel.

## Usage

### Wall boundaries

The set of boundaries that represent solid walls is specified with the `walls` parameter:
these are the boundaries to which distances are computed.

- `walls`: a list of boundary names that correspond to solid walls.

Internally, MOOSE converts the boundary names in `walls` to boundary IDs and loops over all active
elements on those boundaries. For each such boundary face, the centroid position is used as
$\mathbf{x}_w$ in the distance calculation.

If any of the supplied boundary names is invalid (i.e., not present in the mesh), an error is raised.

### Mesh and discretization requirements

`WallDistanceAux` has a few important limitations:

- Replicated meshes only: the object explicitly checks that the mesh is replicated and
  will error out if run on a distributed mesh.
- Finite volume variables: the auxiliary variable associated with this kernel must be a
  finite volume variable (`MooseVariableFV<Real>`). If a non-FV variable is supplied, the kernel
  emits a parameter error indicating that it is currently programmed to use FV machinery only.

Standard MOOSE `block` and `boundary` restrictions can be used, as with any other `AuxKernel`,
to limit where the distance field is computed.

!syntax parameters /AuxKernels/WallDistanceAux

!syntax inputs /AuxKernels/WallDistanceAux

!syntax children /AuxKernels/WallDistanceAux
