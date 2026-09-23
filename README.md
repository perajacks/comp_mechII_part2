# Using FEM for solving the deformation on a plate with Kirchhoff method.
## Plate specifications:
    Thickness t = 0.01m 
    Presure applied q = -0.7 kPa on z direction
    Young modulues E =  210 x 10^3 MPa
    Poison ν = 0.3
Below are the dimension and geometry of the plate:

<img width="474" height="257" alt="Screenshot 2026-09-23 at 21 40 29" src="https://github.com/user-attachments/assets/711fd833-b7bf-409f-9889-74e17ded7e98" />

The governing equetion that discribes the fenomenon is:

$$\frac{\partial^4 w}{\partial x^4} + 2\frac{\partial^4 w}{\partial x^2 \partial y^2} + \frac{\partial^4 w}{\partial y^4} = \frac{q(x,y)}{D_k}$$

where $D_k = \dfrac{Et^3}{12(1-\nu^2)}$ is the plate bending stiffness.

## The aprouch given by the code.
We use four noded rectangulare elements with each node given three degris of fredom: displasment over z axis and rotasion over y and y axis. On the problem solved
the BCs are on the perimeter of the plate and and fixing only the displasment on the z axis.

## Solution to the problem.
For the 12 dof of each element we use pascals polyonym [w] = [x][a] where $$w(x,y) = \alpha_1 + \alpha_2 x + \alpha_3 y + \alpha_4 x^2 + \alpha_5 xy + \alpha_6 y^2 + \alpha_7 x^3 + \alpha_8 x^2y + \alpha_9 xy^2 + \alpha_{10} y^3 + \alpha_{11} x^3y + \alpha_{12} xy^3$$

$$w(x,y) = [x]\{\alpha\}$$

where [x] is the Pascal polynomial vector and $\{\alpha\}$ are the 12 unknown coefficients.

As with most problems using finite element method we will constract a stifnes matrix K. But firstly we should calculate the deformation matrix k. 
### Displacement Matrix

$$k_x = -\frac{\partial^2 w}{\partial x^2}, \qquad k_y = -\frac{\partial^2 w}{\partial y^2}, \qquad k_{xy} = -2\frac{\partial^2 w}{\partial x \partial y}$$

and $w = [x]\{\alpha\}$, $\{\alpha\} = [A]^{-1}\{d\}$, the curvature vector is:

$$\{k\} = \begin{Bmatrix} k_x; \\ k_y; \\ k_{xy}; \end{Bmatrix} \{\alpha\} \implies \{k\} = [\beta]\{\alpha\}$$

Therefore:

$$\{k\} = [\beta][A]^{-1}\{d\} \quad \Rightarrow \quad \{k\} = [B_\tau] \ \{d\}$$

where $[B_\tau] = [\beta][A]^{-1}$ is the **strain-displacement matrix**.























