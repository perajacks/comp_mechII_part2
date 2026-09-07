from pre_processor import Mesh, Material, BoundaryConditions
from solver import KirchhoffPlateElement, Assembler, Solver, Fe
from post_processor import PostProcessor
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    mesh_dimension = 9  # in nodes
    
    nodes = []
    
    for i in range(mesh_dimension):
        for j in range(mesh_dimension):
            s = j/(mesh_dimension-1)
            t = i/(mesh_dimension-1)
            x = s + 0.866*t
            y = 0.5*t
            nodes.append((x,y))
            
    elements = []
    nxt_row = mesh_dimension          
    element_dim = mesh_dimension - 1 
    
    for i in range(element_dim):       
        for j in range(element_dim):   
            n0 = i*nxt_row + j
            n1 = n0 + 1
            n2 = n0 + nxt_row
            n3 = n2 + 1
            elements.append((n0, n1, n3, n2))  
    
    mesh = Mesh(nodes, elements)
    mat = Material(E=210e9, nu=0.3, t=0.01)
    element = KirchhoffPlateElement(mesh, mat)
    asm = Assembler(mesh, element)
   
    K = asm.assemble_stiffness()
    bc = BoundaryConditions()
    
    q = -700.0
    fe_obj = Fe(mesh, mat)
    F = asm.assemble_load(fe_obj, q)   # σωστό, πλήρες assembly του global force vector
    
    # Simply supported σε όλη την περίμετρο: μόνο το w (dof 0) μηδενίζεται,
    # οι στροφές (dofs 1,2) μένουν ελεύθερες -> παραμένουν σχολιασμένες.
    for i in range(mesh_dimension):
        bc.fixed.append((i,0))
        #bc.fixed.append((i,1))
        #bc.fixed.append((i,2))
        bc.fixed.append((mesh_dimension**2 - mesh_dimension + i,0))
        #bc.fixed.append((mesh_dimension**2 - mesh_dimension + i,1))
        #bc.fixed.append((mesh_dimension**2 - mesh_dimension + i,2))
        if i != 0 and i != mesh_dimension-1:
            bc.fixed.append((mesh_dimension*i,0))
            #bc.fixed.append((mesh_dimension*i,1))
            #bc.fixed.append((mesh_dimension*i,2))
            bc.fixed.append((mesh_dimension*i+mesh_dimension-1,0))
            #bc.fixed.append((mesh_dimension*i+mesh_dimension-1,1))
            #bc.fixed.append((mesh_dimension*i+mesh_dimension-1,2))
    
    K, F = bc.apply(K, F)
    solver = Solver()
    U = solver.solve(K, F)
    post = PostProcessor(mesh)
    d = post.export_displacements(U)
    print("Displacements:\n", d)
    
    # μόνο τα w DOFs
    w = U[0::3]
    
    x = np.array([p[0] for p in nodes])
    y = np.array([p[1] for p in nodes])
    
    plt.figure(figsize=(6,5))
    plt.tripcolor(x, y, w, shading='flat', cmap='viridis')
    plt.colorbar(label="w (m)")
    plt.scatter(x, y, c='k', s=8)
    plt.gca().set_aspect('equal')
    plt.title("Plate displacement (top view)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
  