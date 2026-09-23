# Using FEM for solving the deformation on a plate with Kirchhoff method.
## Plate specifications:
    Thickness t = 0.01m 
    pressure applied q = -0.7 kPa on z direction
    Young's modulus E =  210 x 10^3 MPa
    Poisson ν = 0.3
Below are the dimension and geometry of the plate:

<img width="474" height="257" alt="Screenshot 2026-09-23 at 21 40 29" src="https://github.com/user-attachments/assets/711fd833-b7bf-409f-9889-74e17ded7e98" />

The governing equation that describes the phenomenon is:

$$\frac{\partial^4 w}{\partial x^4} + 2\frac{\partial^4 w}{\partial x^2 \partial y^2} + \frac{\partial^4 w}{\partial y^4} = \frac{q(x,y)}{D_k}$$

where $D_k = \dfrac{Et^3}{12(1-\nu^2)}$ is the plate bending stiffness.

## The approach given by the code.
We use four noded rectangular elements with each node given three degrees of freedom: displacement over z axis and rotation over y and y axis. On the problem solved
the BCs are on the perimeter of the plate and and fixing only the displasment on the z axis.

## Solution to the problem.
For the 12 dof of each element we use pascals polynomial [w] = [x][a] where $$w(x,y) = \alpha_1 + \alpha_2 x + \alpha_3 y + \alpha_4 x^2 + \alpha_5 xy + \alpha_6 y^2 + \alpha_7 x^3 + \alpha_8 x^2y + \alpha_9 xy^2 + \alpha_{10} y^3 + \alpha_{11} x^3y + \alpha_{12} xy^3$$

$$w(x,y) = [x]\{\alpha\}$$

where [x] is the Pascal polynomial vector and $\{\alpha\}$ are the 12 unknown coefficients.

As with most problems using finite element method we will construct a stiffness matrix K. But firstly we should calculate the deformation matrix k. 
### Curvature-Displacement Matrix

$$k_x = -\frac{\partial^2 w}{\partial x^2}, \qquad k_y = -\frac{\partial^2 w}{\partial y^2}, \qquad k_{xy} = -2\frac{\partial^2 w}{\partial x \partial y}$$

and $w = [x]\{\alpha\}$, $\{\alpha\} = [A]^{-1}\{d\}$, the curvature vector is:

$$\{k\} = \begin{Bmatrix} k_x \\ k_y \\ k_{xy} \end{Bmatrix} = [\beta]\{\alpha\}$$

Therefore:

$$\{k\} = [\beta][A]^{-1}\{d\} \quad \Rightarrow \quad \{k\} = [B_\tau] \ \{d\}$$

where $[B_\tau] = [\beta][A]^{-1}$ is the **strain-displacement matrix**.

### Stiffness Matrix $[k_e]$

**Stress-strain relation:** $\{\sigma\} = [E]\{\varepsilon\}$, where $\{\varepsilon\} = z\{k\}$

$$[E] = \frac{E}{1-\nu^2} \begin{bmatrix} 1 & \nu & 0; \\ \nu & 1 & 0; \\ 0 & 0 & \frac{1-\nu}{2}; \end{bmatrix}$$

**Principle of Virtual Work:**

$$W_{int} = \int_{V_e} \{\bar{\varepsilon}\}^T [E] \{\varepsilon\} \, dV_e$$

Integrating through the thickness:

$$W_{int} = \int_{A_e} \{\bar{k}\}^T [E_k] \{k\} \, dA_e$$

where:

$$[E_k] = \frac{t^3}{12}[E] = D_k \begin{bmatrix} 1 & \nu & 0; \\ \nu & 1 & 0 ;\\ 0 & 0 & \frac{1-\nu}{2} \end{bmatrix}, \qquad D_k = \frac{Et^3}{12(1-\nu^2)}$$

**Element stiffness matrix:** Given $\{k\} = [B_\tau]\{d\}$:

$$[k_e] = \int_{A_e} [B_\tau]^T [E_k] [B_\tau] \, dA_e = \int_{-b}^{+b}\int_{-a}^{+a} [B_\tau]^T [E_k] [B_\tau] \, dx\, dy$$

where $[B_\tau] = [\beta][A]^{-1}$. 

The integrals are evaluated using Gauss quadrature.

Finally for the Moment vector:

$$\{M\} = \begin{bmatrix} M_x & M_y & M_{xy} \end{bmatrix}^T = [E_k][B_\tau]\{d\}$$

And from that is easy to calculate the stress vector:

$$\{\sigma\} = [E]\,z\,[B_\tau]\{d\} = [S]\{d\}$$

where $[S] = z[E][B_\tau]$ is the stress matrix.

# Results

The program runs on mesh made by 4x4, 8x8, 16x16, 32x32, 54x54, 76x76, 96x96 elements.
As the mesh gets denser the results converge for maximum displacement and maximum principal stress at the bottom of the plate.
We calculate those numbers only on the node that diagonals of the shape meet.

<img width="993" height="705" alt="image" src="https://github.com/user-attachments/assets/4fc92630-97d4-470e-8641-8b25715fb87b" />


<img width="993" height="705" alt="image" src="https://github.com/user-attachments/assets/f48b610c-f591-4033-b642-82c989aab70b" />

We use as a reference for the maximum principal stress a value given by NAFEMS LSB1 as σ* = 0.802 MPa.
We can see that the program nears to the given value even though it never reaches it.

The total displacement graph for the whole structure is:

<img width="1145" height="993" alt="image" src="https://github.com/user-attachments/assets/f0666f78-7429-417d-9ad0-23734cfcf0aa" />

The displacements are following what we would expect according to the boundary conditions and pressure affecting the plate.

Below are given the maximum values given by the different meshes:
   element_dim   n_elem   w_max (mm)  s_max (MPa)   |σ|/σ*
--------------------------------------------------------
           4       16      0.00286       0.2682    0.334
           8       64      0.00772       0.4566    0.569
          16      256      0.01119       0.6238    0.778
          32     1024      0.01290       0.7116    0.887
          54     2916      0.01361       0.7466    0.931
          76     5776      0.01392       0.7612    0.949
          96     9216      0.01408       0.7686    0.958


















