import numpy as np
import sympy as sm

class KirchhoffPlateElement:
    def __init__(self, mesh, material):
        self.mesh = mesh
        self.mat = material

    def ke(self, elem_nodes):
        ke = np.zeros((12, 12))
  
        x, y = sm.symbols('x y')
        D = float(self.mat.D)
        nu = float(self.mat.nu)
        
        Ek = float(D) * sm.Matrix([
          [1.0,  nu, 0.0],
          [nu,  1.0, 0.0],
          [0.0, 0.0, (1.0 - nu)/2.0]])
          
        a=1
        b=1
        
        A = np.array([
        [1, -a, -b,  a**2,  a*b,  b**2,  -a**3,  -a**2*b,  -a*b**2,  -b**3,   a**3*b,   a*b**3],
        [0,  0,  1,   0,   -a,   -2*b,     0,      a**2,     2*a*b,   3*b**2,  -a**3,   -3*a*b**2],
        [0,  1,  0,  -2*a, -b,    0,    3*a**2,   2*a*b,    b**2,      0,    -3*a**2*b,  -b**3],
        [1,  a, -b,  a**2, -a*b,  b**2,   a**3,   -a**2*b,   a*b**2,  -b**3,  -a**3*b,  -a*b**3],
        [0,  0,  1,   0,    a,   -2*b,     0,      a**2,    -2*a*b,   3*b**2,   a**3,    3*a*b**2],
        [0,  1,  0,   2*a, -b,    0,    3*a**2,  -2*a*b,    b**2,      0,    -3*a**2*b,  -b**3],
        [1,  a,  b,  a**2,  a*b,  b**2,   a**3,    a**2*b,   a*b**2,   b**3,   a**3*b,   a*b**3],
        [0,  0,  1,   0,    a,    2*b,     0,      a**2,     2*a*b,   3*b**2,   a**3,    3*a*b**2],
        [0,  1,  0,   2*a,  b,    0,    3*a**2,   2*a*b,    b**2,      0,     3*a**2*b,   b**3],
        [1, -a,  b,  a**2, -a*b,  b**2,  -a**3,    a**2*b,  -a*b**2,   b**3,  -a**3*b,  -a*b**3],
        [0,  0,  1,   0,   -a,    2*b,     0,      a**2,    -2*a*b,   3*b**2,  -a**3,   -3*a*b**2],
        [0,  1,  0,  -2*a,  b,    0,    3*a**2,  -2*a*b,    b**2,      0,     3*a**2*b,   b**3]

        ], dtype=float)


        invA  = sm.Matrix(np.linalg.inv(A))


        b_matix = sm.Matrix([
            [0, 0, 0, 2, 0, 0, 6*x, 2*y, 0, 0, 6*x*y, 0],
            [0, 0, 0, 0, 0, 2, 0, 0, 2*x, 6*y, 0, 6*x*y],
            [0, 0, 0, 0, 2, 0, 0, 4*x, 4*y, 0, 6*x**2, 6*y**2]

        ], dtype=float)

        BT = b_matix * invA
        Ke_smbl = (BT.T)*Ek*BT
        Ke = sm.integrate(Ke_smbl, (x, -a, a))
        Ke = sm.integrate(Ke, (y, -b, b))
        ke = Ke
        return ke

class Assembler:
    def __init__(self, mesh, element):
        self.mesh = mesh
        self.element = element
        

    def assemble_stiffness(self):
        K = np.zeros((self.mesh.ndof, self.mesh.ndof))
        for conn in self.mesh.elements:
            elem_coords = self.mesh.nodes[conn]
            ke = self.element.ke(elem_coords)
            edofs = []
            for nid in conn:
                base = nid * self.mesh.ndof_per_node
                edofs.extend([base, base+1, base+2])
            for i in range(12):
                for j in range(12):
                    print(ke[i, j])
                    K[edofs[i], edofs[j]] += ke[i, j]
        return K

class Solver:
    def solve(self, K, F):
        return np.linalg.solve(K, F)
