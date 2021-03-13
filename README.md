# iga-matlab

## What is this?

This repository contains a set of matlab functions to do geoemtric nonlinear Isogeometric analysis.

### Why was it done?

This was developed as part of my master thesis in structural engineering at Chalmers University of Technology in the fall of 2018. Had previously experimented with actively bent forms, and in this thesis I wanted to investigate the use of nonlinear IGA as a means to design them. The core aim was to find a geometrically consistent (NURBS) way of relating deformed and undeformed design geometries to each other accounting for the stiffness distribution.

More information can be found in the final [thesis report](https://hdl.handle.net/20.500.12380/301616).

It builds upon the work done by S. Almstedt and P. Safari Hesari in the their [thesis](https://hdl.handle.net/20.500.12380/301616)

### What is IGA?

*Nurbs geometry vs a polygonal mesh.*

***

### Functionality

- **Kirchhoff-Love Shells**

The main target was the geometrically nonlinear analysis of single IGA patches.

![](https://github.com/joelhi/IGA_MATLAB/blob/master/KL%20Shell/Resources/Surface_e11.gif)


- **Euler-Bernoulli Beams**

Matrix and Dynamic Relaxation.

![alt text](https://github.com/joelhi/IGA_MATLAB/blob/master/Beam/Gifs/Elastica1.gif)

![alt text](https://github.com/joelhi/IGA_MATLAB/blob/master/Beam/Gifs/DR_Faster.gif)

***

### General Disclaimer
