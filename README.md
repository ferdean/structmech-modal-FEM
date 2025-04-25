# FEM modal solver (structural mechanics)

> Note: This project corresponds to the author's BSc thesis

Linear FEM solver for structural mechanics applications. Solves a modal problem on a given externally generated mesh [1], allowing multimaterial single-body domains.

[1] The mesh must be input as a MATLAB structure with two fields, an array containing the coordinates of the nodes and an array containing the topology of the mesh (i.e. which nodes form each of the elements). In the *test* folder you can find examples of the architecture of these structures. If you have any questions, do not hesitate to contact me at deandresvertferran(at)gmail(dot)com.

Mass and stiffness matrices correspond to the force balance equation containing information of the N degrees of freedom of the system, this is 

$\mathbf{M \ddot{U}}+\mathbf{K U}=\mathbf{F}$

The current version includes: 

* *lib*, including all the assembly, solve and visualization function and methods.
* *tests*, including meshes, geometries and material properties to serve as example of application. The example corresponds to a 3D model of a railway wheel manipulated with a viscoelastic layer constrained with an aluminum layer to reduce the surface mechanical vibration. 
* *plots*, including some plotting examples.

### Examples of mode representation: 

* **Railway wheel case:**

![example_mode_Wheel](https://user-images.githubusercontent.com/92535468/164887157-5175ffca-6b09-4fce-a824-41712d63ef33.png)

* **Simple plate case:**

![example_mode_Plate](https://user-images.githubusercontent.com/92535468/164887158-10cd4958-8b0e-49a7-bc92-c11262909028.png)

