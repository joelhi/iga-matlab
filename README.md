# iga-matlab

### What is this?

This repository contains a set of MATLAB functions to do geoemtric nonlinear Isogeometric analysis.

### Why was it done?

This was developed as part of my master thesis in structural engineering at Chalmers University of Technology in the fall of 2018. Had previously experimented with actively bent forms, and in this thesis I wanted to investigate the use of nonlinear IGA as a means to design them. The core aim was to find a geometrically consistent (NURBS) way of relating deformed and undeformed design geometries to each other, accounting for the stiffness distribution.

More information can be found in the final [thesis report](https://hdl.handle.net/20.500.12380/301616).

It builds upon the work done by S. Almstedt and P. Safari Hesari in the their [thesis](https://hdl.handle.net/20.500.12380/301616)

### What is IGA?

*Nurbs geometry vs a polygonal mesh.*

***

### Functionality

- **Kirchhoff-Love (KL) Shells**

The main target was the geometrically nonlinear analysis of single IGA patches.

This repository contains a single patch implementation.

As a design feature, allows the input of an independent thickness parameterization, which can be used to manipulate

![](https://github.com/joelhi/IGA_MATLAB/blob/master/KL%20Shell/Resources/Surface_e11.gif)


- **Euler-Bernoulli (EB) Beams**

As a first step, prior to the shell implementation, and simplified version was made, as a proving ground. This was derived in a similar fashion to the shell, but for a 1-dimensiona object in 2-dimensional space. Effectively making it an Euler-Bernoulli beam (E-B).

This implementation was made as a test, to try different solution schemes etc. and is probably not the best IGA implementation of an E-B beam.

***

Below are two examples; the first is a displacement controlled solution to a simple elastica curve, and the second is a dynamic relaxation solution.

![alt text](https://github.com/joelhi/IGA_MATLAB/blob/master/Beam/Gifs/Elastica1.gif)

![alt text](https://github.com/joelhi/IGA_MATLAB/blob/master/Beam/Gifs/DR_Faster.gif)

***

### General Disclaimer

The implementation does have some issues, so should be used with caution.

Some examples are:

- **Membrane locking for KL shells** *need a very high degree (>10) for the locking not to be substantial*
- **Beam implementation not verified** *Sometimes looks to behave in a too stiff fashion, however, buckling load corresponds well to analytical solution.*

But if you are interested in IGA, this repo can serve as a good point of entry.
