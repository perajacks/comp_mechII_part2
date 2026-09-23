import numpy as np

# weights in gauss intergration are 1 as a result we will not calculate them

class KirchhoffPlateElement:
    def __init__(self, mesh, material):
        self.mesh = mesh
        self.mat = material

    def ke(self, elem_nodes):
        ke = np.zeros((12, 12))

        D = float(self.mat.D)
        nu = float(self.mat.nu)
        
        Ek = float(D) * np.array([
          [1.0,  nu, 0.0],
          [nu,  1.0, 0.0],
          [0.0, 0.0, (1.0 - nu)/2.0]])
          
              

        x1, y1 = elem_nodes[0]
        x2, y2 = elem_nodes[1]
        x4, y4 = elem_nodes[3]
        
        xi  = np.array([1/np.sqrt(3),-1/np.sqrt(3)])
        hta = np.array([1/np.sqrt(3),-1/np.sqrt(3)])

        
        
        ax = x2 - x1
        ay = y2 - y1

        bx = x4 - x1
        by = y4 - y1
       
        
        a = 1 
        b = 1


        Delta = ax*by - ay*bx
        detJ = np.abs(Delta) / 4.0

        invJ = (2.0/Delta) * np.array([
           [by, -bx],
           [-ay, ax]
            ])
        
        
        JinvT = np.transpose(invJ)

        Tn = np.array([
           [1, 0, 0],
           [0, -JinvT[0,0], -JinvT[0,1]],
           [0, -JinvT[1,0], -JinvT[1,1]]
        ])
                                       


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


        invA  = np.linalg.inv(A)

       


        Te = np.block([
        [Tn, np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3))],
        [np.zeros((3,3)), Tn, np.zeros((3,3)), np.zeros((3,3))],
        [np.zeros((3,3)), np.zeros((3,3)), Tn, np.zeros((3,3))],
        [np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3)), Tn]
        ])
        
        Te_inv = np.linalg.inv(Te)
        
        TB = (4.0 / (Delta**2)) * np.array([
       [by**2,      ay**2,      -ay*by],
       [bx**2,      ax**2,      -ax*bx],
       [-2*bx*by,  -2*ax*ay,    ax*by + ay*bx]
         ], dtype=float)

  
        for i in  xi:
            for j in hta:
                
                x = a * i
                y = b * j
                
                b_matix = np.array([
                    [0, 0, 0, 2, 0, 0, 6*x, 2*y, 0, 0, 6*x*y, 0],
                    [0, 0, 0, 0, 0, 2, 0, 0, 2*x, 6*y, 0, 6*x*y],
                    [0, 0, 0, 0, 2, 0, 0, 4*x, 4*y, 0, 6*x**2, 6*y**2]
        
                ], dtype=float)
                
                BT = b_matix @ invA
                
                Bx = TB @ BT @ Te_inv
        

                ke = ke + (np.transpose(Bx) @ Ek @ Bx) * abs(detJ)
                
           
        
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
                    K[edofs[i], edofs[j]] += ke[i, j]
        return K
    
    
    def assemble_load(self, fe_obj, q):         
        F = np.zeros(self.mesh.ndof)
        for conn in self.mesh.elements:
            elem_coords = self.mesh.nodes[conn]
            fe_local = fe_obj.fe(elem_coords, q)
            edofs = []
            for nid in conn:
                base = nid * self.mesh.ndof_per_node
                edofs.extend([base, base+1, base+2])
            for i in range(12):
                F[edofs[i]] += fe_local[i]
        return F

class Solver:
    def solve(self, K, F):
        return np.linalg.solve(K, F)
    
class Fe:
    def __init__(self, mesh, material):
        self.mesh = mesh
        self.mat = material
      
    def fe(self,elem_nodes, q):
        
        x1, y1 = elem_nodes[0]
        x2, y2 = elem_nodes[1]
        x4, y4 = elem_nodes[3]
        
        xi  = np.array([1/np.sqrt(3),-1/np.sqrt(3)])
        hta = np.array([1/np.sqrt(3),-1/np.sqrt(3)])

        
        
        ax = x2 - x1
        ay = y2 - y1

        bx = x4 - x1
        by = y4 - y1
       
        
        a = 1 
        b = 1
     

        Delta = ax*by - ay*bx
        detJ = np.abs(Delta) / 4.0
       

        
        
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


        invA  = np.linalg.inv(A)
        
        trans_inv_A = np.transpose(invA)
        
        g = np.zeros(12, dtype=float)
        
    
        
        for r in xi:
           for s in hta:

               x = a * r
               y = b * s
    
               X = np.array([
                   1, x, y, x**2, x*y, y**2,
                   x**3, x**2*y, x*y**2, y**3,
                   x**3*y, x*y**3
               ], dtype=float)
    
               g += X * q * abs(detJ)   
    
        fe =  trans_inv_A  @ g
        return fe
    
class Stress:
    def moments_at_point(elem_coords, d_disp, mat, x_local, y_local):

        
        x1, y1 = elem_coords[0]; x2, y2 = elem_coords[1]; x4, y4 = elem_coords[3]
        ax = x2-x1; ay = y2-y1; bx = x4-x1; by = y4-y1
        
        a = 1 
        b = 1
        
        Delta = ax*by - ay*bx

        invJ = (2.0/Delta) * np.array([[by, -bx], [-ay, ax]])
        JinvT = np.transpose(invJ)
        Tn = np.array([[1,0,0],[0,-JinvT[0,0],-JinvT[0,1]],[0,-JinvT[1,0],-JinvT[1,1]]])
        Te = np.block([[Tn,np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3))],
                       [np.zeros((3,3)),Tn,np.zeros((3,3)),np.zeros((3,3))],
                       [np.zeros((3,3)),np.zeros((3,3)),Tn,np.zeros((3,3))],
                       [np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3)),Tn]])
        Te_inv = np.linalg.inv(Te)

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
        invA = np.linalg.inv(A)

        x, y = x_local, y_local
        b_matrix = np.array([
            [0,0,0,2,0,0,6*x,2*y,0,0,6*x*y,0],
            [0,0,0,0,0,2,0,0,2*x,6*y,0,6*x*y],
            [0,0,0,0,2,0,0,4*x,4*y,0,6*x**2,6*y**2]
        ], dtype=float)
        
        TB = (4.0/Delta**2) * np.array([
        [by**2,      ay**2,      -ay*by],
        [bx**2,      ax**2,      -ax*bx],
        [-2*bx*by,  -2*ax*ay,    ax*by+ay*bx]
            ])


        BT = b_matrix @ invA
        Bx = TB @ BT @ Te_inv
        kappa = Bx @ d_disp                       # [kxx, kyy, 2kxy]

        nu = mat.nu
        Ek = mat.D * np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])
        return Ek @ kappa                            # [Mxx, Myy, Mxy]
            



    def principal_stresses(sx, sy, txy):
        center = 0.5*(sx+sy)
        radius = np.sqrt((0.5*(sx-sy))**2 + txy**2)
        return center+radius, center-radius          # sigma_1 (max), sigma_2 (min)

