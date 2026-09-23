## Using FEM for solving the deformation on a plate with Kirchhoff method.
# Plate specifications:
    Thickness t = 0.01m 
    Presure applied q = -0.7 kPa on z direction
    Young modulues E =  210 x 10^3 MPa
    Poison ν = 0.3
Below are the dimension and geometry of the plate:

<img width="474" height="257" alt="Screenshot 2026-09-23 at 21 40 29" src="https://github.com/user-attachments/assets/711fd833-b7bf-409f-9889-74e17ded7e98" />

The governing equetion that discribes the fenomenon is:

$$\frac{\partial^4 w}{\partial x^4} + 2\frac{\partial^4 w}{\partial x^2 \partial y^2} + \frac{\partial^4 w}{\partial y^4} = \frac{q(x,y)}{D_k}$$

where $D_k = \dfrac{Et^3}{12(1-\nu^2)}$ is the plate bending stiffness.

# The aprouch given by the code.
We use four noded rectangulare elements with each node given three degris of fredom: displasment over z axis and rotasion over y and y axis. On the problem solved
the BCs are on the perimeter of the plate and and fixing only the displasment on the z axis.





















