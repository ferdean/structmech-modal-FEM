# FEM modal solver (structural mechanics)

Linear FEM solver for structural mechanics applications. Solves a modal problem on a given externally generated mesh, allowing multimaterial single-body domains.

Mass and stiffness matrices correspond to the force balance equation containing information of the N degrees of freedom of the system, this is 

$\mathbf{M \ddot{U}}+\mathbf{K U}=\mathbf{F}$

The current version includes: 

* *lib*, including all the assembly, solve and visualization function and methods.
* *tests*, including meshes, geometries and material properties to serve as example of application. The example corresponds to a 3D model of a railway wheel manipulated with a viscoelastic layer constrained with an aluminum layer to reduce the surface mechanical vibration. 
* *plots*, including some plotting examples.

