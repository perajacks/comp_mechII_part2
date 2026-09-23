from pre_processor import Mesh, Material, BoundaryConditions
from solver import KirchhoffPlateElement, Assembler, Solver,Fe, Stress
from post_processor import PostProcessor
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

def run_mesh(mesh_dimension):
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
    F = asm.assemble_load(fe_obj, q)
    

    
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

    
    # μόνο τα w DOFs
    w_all = U[0::3]
    w_max = np.max(np.abs(w_all))
    
    x = np.array([p[0] for p in nodes])
    y = np.array([p[1] for p in nodes])
    
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    
    triang = mtri.Triangulation(x, y)
    
    ax.plot_trisurf(triang, w_all, cmap='viridis', edgecolor='none', alpha=0.9)
    
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('w (m)')
    ax.set_title(f'Plate displacement (3D) elemnt_mesh:, {mesh_dimension-1}×{mesh_dimension-1}')
    
    mappable = plt.cm.ScalarMappable(cmap='viridis')
    mappable.set_array(w_all)
    fig.colorbar(mappable, ax=ax, shrink=0.5, label='w (m)')
    
    plt.tight_layout()
    plt.show()
  

#--------------------
#Stress calcuating
#--------------------
    
    i_c = mesh_dimension // 2
    j_c = mesh_dimension // 2
    elem_index = i_c * element_dim + j_c
    
    conn = elements[elem_index]              
    elem_coords = mesh.nodes[list(conn)]
    
    d_local = np.zeros(12)
    for k, nid in enumerate(conn):
        d_local[3*k : 3*k+3] = d[nid, :]     # d[nid,:] = [w, θx, θy] 
    
    x1, y1 = elem_coords[0]
    x2, y2 = elem_coords[1]
    x4, y4 = elem_coords[3]
    
    
    M = Stress.moments_at_point(elem_coords, d_local, mat, -1, -1)   
    t = mat.t
    sigma_bottom = -6.0 * M / t**2      # [sxx, syy, txy] 
    s1, s2 = Stress.principal_stresses(*sigma_bottom)
    s_max = max(abs(s1), abs(s2))
    
    return w_max, s_max, w_all, nodes, elements



if __name__ == "__main__":

    mesh_dims = [5, 9, 17, 33, 55, 77, 97]
    results = []


    print("element_dime  n_elem     w_max         s_max")
    for md in mesh_dims:
        w_max, s_max, w_all, nodes, elements = run_mesh(md)
        n_elem = (md-1)**2
        results.append((md, n_elem, w_max, s_max))
        print(f"{md-1:>6} {n_elem:>8} {w_max*1e3:>12.5f} {s_max/1e6:>12.4f}")


    n_elems = [r[1] for r in results]
    s_maxs  = [r[3]/1e6 for r in results]

    plt.figure(figsize=(7, 5))
plt.plot(n_elems, s_maxs, 's-', color='coral', label='FEM')
plt.axhline(0.802, color='gray', linestyle='--', label='NAFEMS 0.802 MPa')
plt.xlabel('Αριθμός στοιχείων')
plt.ylabel('σ_max (MPa)')
plt.title('Convergence of maximum stress')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()