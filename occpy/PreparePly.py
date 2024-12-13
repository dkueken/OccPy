import numpy as np
import OSToolBox as ost


def prepare_ply(vox_dim, PlotDim, data):
    data = array3Dto2D(data, vox_dim, PlotDim)
    # Generate mesh data
    verts, faces = generate_mesh_data(data, vox_dim, PlotDim)

    return verts, faces

def array3Dto2D(data, vox_dim, PlotDim):
    x = np.arange(PlotDim['minX'], PlotDim['maxX'], vox_dim)
    y = np.arange(PlotDim['minY'], PlotDim['maxY'], vox_dim)
    z = np.arange(PlotDim['minZ'], PlotDim['maxZ'], vox_dim)

    x, y, z = np.meshgrid(x, y, z)
    x, y, z, data = x.flatten(), y.flatten(), z.flatten(), data.flatten()
    data = np.c_[x,y,z,data]
    mask = np.ma.masked_where(data[:,-1]>0, data[:,-1]).mask
    data = data[mask]
    return data


def calculate_voxel_corners(data, vox_dim, PlotDim):
    offset = vox_dim/2.
    xyz_offsets = np.array([
        [-offset, -offset, -offset],
        [offset, -offset, -offset],
        [offset, offset, -offset],
        [-offset, offset, -offset],
        [-offset, -offset, offset],
        [offset, -offset, offset],
        [offset, offset, offset],
        [-offset, offset, offset]
    ])

    # Fetch scalars array and repeat each line the number of lines in xyz_offsets
    scalars = data[:,3:]
    scalars = np.repeat(scalars, xyz_offsets.shape[0], axis=0)

    # Add offsets for each points
    xyz = data[:,:3][:, np.newaxis] + xyz_offsets
    xyz = np.vstack(xyz)

    # Concatenate the points with their respective scalars
    xyz_scalars = np.c_[xyz, scalars]

    return xyz_scalars


def generate_mesh_data(data, vox_dim, PlotDim):
    # Calculate voxel corners
    corners = calculate_voxel_corners(data, vox_dim, PlotDim)

    # Define the face indices for a cube
    face_indices_0 = np.array([
        [0, 1, 2], [0, 2, 3],  # Front face
        [4, 5, 6], [4, 6, 7],  # Back face
        [0, 1, 5], [0, 5, 4],  # Bottom face
        [2, 3, 7], [2, 7, 6],  # Top face
        [1, 2, 6], [1, 6, 5],  # Right face
        [3, 0, 4], [3, 4, 7]   # Left face
    ])

    # Total number of voxels
    num_voxels = int(corners.shape[0] / 8)

    # Repeat face indices for each voxel
    faces = np.tile(face_indices_0, (num_voxels, 1))

    # Offset face indices for each voxel
    offset = np.arange(0, num_voxels * 8, 8).repeat(12)
    faces += offset[:, np.newaxis]

    return corners, faces