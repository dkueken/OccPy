import os
import numpy as np
import open3d as o3d
import pyvista as pv
import laspy


def aabb_to_mesh(aabb):
    min_x, min_y, min_z = aabb.get_min_bound()
    max_x, max_y, max_z = aabb.get_max_bound()

    # vertices
    V = np.array([
        [min_x, min_y, min_z],
        [max_x, min_y, min_z],
        [max_x, max_y, min_z],
        [min_x, max_y, min_z],
        [min_x, min_y, max_z],
        [max_x, min_y, max_z],
        [max_x, max_y, max_z],
        [min_x, max_y, max_z],
    ], dtype=np.float64)

    F = np.array([
        [0, 1, 2], [0, 2, 3],       # bottom (z=min)
        [4, 5, 6], [4, 6, 7],       # top (z=max)
        [0, 1, 5], [0, 5, 4],       # y=min face
        [3, 2, 6], [3, 6, 7],       # y=max face
        [0, 3, 7], [0, 7, 4],       # x=min face
        [1, 2, 6], [1, 6, 5],       # x=max face
    ], dtype=np.int32)

    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(V)
    mesh.triangles = o3d.utility.Vector3iVector(F)
    mesh.compute_vertex_normals()
    return mesh

def batch_aabbs_to_mesh(aabbs: list[o3d.geometry.AxisAlignedBoundingBox]):
    all_vertices = []
    all_triangles = []
    offset = 0

    F = np.array([
        [0, 1, 2], [0, 2, 3],       # bottom (z=min)
        [4, 5, 6], [4, 6, 7],       # top (z=max)
        [0, 1, 5], [0, 5, 4],       # y=min face
        [3, 2, 6], [3, 6, 7],       # y=max face
        [0, 3, 7], [0, 7, 4],       # x=min face
        [1, 2, 6], [1, 6, 5],       # x=max face
    ], dtype=np.int32)

    for aabb in aabbs:
        min_x, min_y, min_z = aabb.get_min_bound()
        max_x, max_y, max_z = aabb.get_max_bound()

        # vertices
        V = np.array([
            [min_x, min_y, min_z],
            [max_x, min_y, min_z],
            [max_x, max_y, min_z],
            [min_x, max_y, min_z],
            [min_x, min_y, max_z],
            [max_x, min_y, max_z],
            [max_x, max_y, max_z],
            [min_x, max_y, max_z],
        ], dtype=np.float64)

        all_vertices.append(V)
        all_triangles.append(F + offset)
        offset += V.shape[0]

    all_vertices = np.vstack(all_vertices)
    all_triangles = np.vstack(all_triangles)

    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(all_vertices)
    mesh.triangles = o3d.utility.Vector3iVector(all_triangles)

    return mesh

def aabb_to_lineset(aabb):
    return o3d.geometry.LineSet.create_from_axis_aligned_bounding_box(aabb)

def batch_aabbs_to_lineset(aabbs: list[o3d.geometry.AxisAlignedBoundingBox]):
    ls = o3d.geometry.LineSet()
    if not aabbs:
        return ls

    tmpl = o3d.geometry.LineSet.create_from_axis_aligned_bounding_box(aabbs[0])
    E = np.asarray(tmpl.lines, dtype=np.int32)

    all_pts, all_lines = [], []
    offset = 0
    for bb in aabbs:
        P = np.asarray(bb.get_box_points())
        all_pts.append(P)
        all_lines.append(E + offset)
        offset += P.shape[0]  # += 8

    ls.points = o3d.utility.Vector3dVector(np.vstack(all_pts))
    ls.lines  = o3d.utility.Vector2iVector(np.vstack(all_lines))
    return ls


def o3d_mesh_to_pyvista(o3d_mesh):
    vertices = np.asarray(o3d_mesh.vertices)
    triangles = np.asarray(o3d_mesh.triangles)
    faces = np.hstack([np.full((triangles.shape[0], 1), 3), triangles]).astype(np.int64).ravel()
    pv_mesh = pv.PolyData(vertices, faces)
    
    if o3d_mesh.has_vertex_colors():
        colors = np.asarray(o3d_mesh.vertex_colors)
        pv_mesh.point_data["Colors"] = (colors * 255).astype(np.uint8)
    
    return pv_mesh
    
def o3d_lineset_to_pyvista(o3d_lineset):
    vertices = np.asarray(o3d_lineset.points)
    lines = np.asarray(o3d_lineset.lines)

    # PyVista expects a flat array like: [2, v0, v1, 2, v0, v1, ...]
    n_lines = lines.shape[0]
    cells = np.hstack([np.full((n_lines, 1), 2), lines]).astype(np.int64).ravel()

    # Create PolyData with line cells
    pv_lines = pv.PolyData()
    pv_lines.points = vertices
    pv_lines.lines = cells
    
    return pv_lines


def test_pyvista():

    # Create a cube mesh
    cube = pv.Cube()
    
    # Initialize the plotter
    plotter = pv.Plotter()
    
    # Add the cube mesh with transparency
    plotter.add_mesh(cube, opacity=0.5, color='red')
    
    # Display the plot
    plotter.show()

def pyvista_bboxs():
    
    aabbs = []
    for i in range(100):
        min_bound = np.random.rand(3) * 10
        max_bound = min_bound + [1.0, 1.0, 1.0]
        aabbs.append(o3d.geometry.AxisAlignedBoundingBox(min_bound, max_bound))
        
    mesh = batch_aabbs_to_mesh(aabbs)
    
    pv_mesh = o3d_mesh_to_pyvista(mesh)
    
    plotter = pv.Plotter()
    
    plotter.add_mesh(pv_mesh, opacity=0.5, color='red')
    
    plotter.show()

def pyvista_lineset():
    aabbs = []
    for i in range(100):
        min_bound = np.random.rand(3) * 10
        max_bound = min_bound + [1.0, 1.0, 1.0]
        aabbs.append(o3d.geometry.AxisAlignedBoundingBox(min_bound, max_bound))
    
    o3d_lines = batch_aabbs_to_lineset(aabbs)
    
    pv_lines = o3d_lineset_to_pyvista(o3d_lines)

    plotter = pv.Plotter()
    plotter.add_mesh(pv_lines, color="black", line_width=3)
    plotter.show()
    
def pyvista_bboxs_with_lines():
    aabbs = []
    for i in range(100):
        min_bound = np.random.rand(3) * 10
        max_bound = min_bound + [1.0, 1.0, 1.0]
        aabbs.append(o3d.geometry.AxisAlignedBoundingBox(min_bound, max_bound))
    
    o3d_lines = batch_aabbs_to_lineset(aabbs)
    
    pv_lines = o3d_lineset_to_pyvista(o3d_lines)
    
    mesh = batch_aabbs_to_mesh(aabbs)
    
    pv_mesh = o3d_mesh_to_pyvista(mesh)

    plotter = pv.Plotter()
    plotter.add_mesh(pv_lines, color="black", line_width=3)
    plotter.add_mesh(pv_mesh, opacity=0.4, color='red')
    plotter.show()
    

HORIZONTAL_BUFFER = 10
VERTICAL_BUFFER_BOTTOM = 5
VERTICAL_BUFFER_TOP = 20
VOX_DIM = 0.1

def occmap_vis_pyvista(occmap_file, pointcloud_file):
    # read npy
    occmap = np.load(occmap_file)
    # switch x and y
    occmap = np.swapaxes(occmap, 0, 1)
    
    las = laspy.read(pointcloud_file)
    points = np.vstack((las.x, las.y, las.z)).transpose()
    min_pc = np.min(points, axis=0)
    max_pc = np.max(points, axis=0)
    
    min_pc -= [HORIZONTAL_BUFFER, HORIZONTAL_BUFFER, VERTICAL_BUFFER_BOTTOM]
    max_pc += [HORIZONTAL_BUFFER, HORIZONTAL_BUFFER, VERTICAL_BUFFER_TOP]
    
    print(min_pc)
    print(max_pc)
    
    vox_dim = 0.1
    
    print((max_pc-min_pc)/vox_dim)
    
    dims = occmap.shape
    
    print(dims)
    
    crop_min = np.floor(np.divide(dims, 2)) - 40
    crop_max = np.floor(np.divide(dims, 2)) + 30
    
    bboxs_occl = []
    bbox_unobserved = []
    bbox_hit = []
    
    # where the fuck do these values come from
    # TODO: redo the occlusion mapping and very clearly save the min_bound and max_bound of the voxelgrid in 3D coordinates
    offset = [+0.2,+0.1,-0.6]
    # offset = [0,0,0]
    
    
    for x in range(int(crop_min[0]), int(crop_max[0])):
        for y in range(int(crop_min[1]), int(crop_max[1])):
            for z in range(int(crop_min[2]), int(crop_max[2])):
                min_bound = min_pc + np.array([x,y,z])*vox_dim + offset
                max_bound = min_bound + vox_dim
                if occmap[x,y,z] == 3:
                    # add mesh cube
                    bboxs_occl.append(o3d.geometry.AxisAlignedBoundingBox(min_bound, max_bound))
                elif occmap[x,y,z] == 4:
                    bbox_unobserved.append(o3d.geometry.AxisAlignedBoundingBox(min_bound, max_bound))
                elif occmap[x,y,z] == 1:
                    bbox_hit.append(o3d.geometry.AxisAlignedBoundingBox(min_bound, max_bound))
                    
    if len(bboxs_occl) > 0:
        mesh_occl = batch_aabbs_to_mesh(bboxs_occl)
        pv_mesh_occl = o3d_mesh_to_pyvista(mesh_occl)
    
    mesh_hit = batch_aabbs_to_mesh(bbox_hit)
    pv_mesh_hit = o3d_mesh_to_pyvista(mesh_hit)
    
    if len(bbox_unobserved) > 0:
        mesh_unobserved = batch_aabbs_to_mesh(bbox_unobserved)
        pv_mesh_unobserved = o3d_mesh_to_pyvista(mesh_unobserved)
    
    plotter = pv.Plotter()
    
    if len(bboxs_occl) > 0:
        plotter.add_mesh(pv_mesh_occl, opacity=0.2, color='red')
    plotter.add_mesh(pv_mesh_hit, opacity=0.1, color='green')
    if len(bbox_unobserved) > 0:
        plotter.add_mesh(pv_mesh_unobserved, opacity=0.2, color='blue')
        
    # add point cloud
    min_pc_crop = min_pc + crop_min*vox_dim + offset
    max_pc_crop = min_pc + crop_max*vox_dim + offset
    print(min_pc_crop)
    print(max_pc_crop)
    mask = np.all((points >= min_pc_crop) & (points <= max_pc_crop), axis=1)
    points_in_crop = points[mask]
    point_cloud = pv.PolyData(points_in_crop)
    # plotter.add_points(point_cloud, style="points", point_size=2, opacity=0.5)
    plotter.add_points(point_cloud, style="points", point_size=2)
    
    plotter.show()
    

def main():
 
    occmap_file = "C:/Users/wcherlet/OneDrive - UGent/PHD/temp/Classification_all.npy"
    pointcloud_file = "C:/Users/wcherlet/OneDrive - UGent/PHD/temp/ABI_2t_1cm_SOR_6_10.las"
    occmap_vis_pyvista(occmap_file, pointcloud_file)

if __name__ == "__main__":
    main()