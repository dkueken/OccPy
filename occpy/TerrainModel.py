import rasterio
import numpy as np
from rasterio.enums import Resampling

from rasterio.transform import from_origin
from scipy.ndimage import gaussian_filter


# Value written to pixels with no hit voxel in their column. Kept as a finite
# sentinel (rather than NaN) so that `rasterio.fill.fillnodata`, as used in
# occpy.util.normalize_occlusion_output, can identify and interpolate them.
NODATA_VALUE = -9999.0

def _elevation_from_nhit_grid(nhit_grid, plot_dim, vox_dim, use_min):
    """
    Collapse an (nx, ny, nz) Nhit grid into a 2D elevation raster by picking,
    for each (x, y) column, the voxel-center height of either the lowest
    (use_min=True) or highest (use_min=False) voxel with at least one hit.

    Returns
    -------
    numpy.ndarray
        2D array of shape (ny, nx) in standard north-up raster order
        (row 0 = maxY, column 0 = minX), with columns that have no hit
        voxel set to NODATA_VALUE.
    """
    min_z = plot_dim[2]
    hit = nhit_grid > 0
    has_hit = np.any(hit, axis=2)

    if use_min:
        z_idx = np.argmax(hit, axis=2)  # index of first True along z
    else:
        nz = hit.shape[2]
        z_idx = nz - 1 - np.argmax(np.flip(hit, axis=2), axis=2)  # index of last True along z

    elevation = min_z + (z_idx.astype(np.float64) + 0.5) * vox_dim
    elevation = np.where(has_hit, elevation, NODATA_VALUE)

    # nhit_grid axes are (x, y, z); rasters are row-major (row=y, col=x)
    # with row 0 = maxY, so transpose to (ny, nx) and flip the y-axis.
    elevation = np.flipud(elevation.T)
    return elevation

def _smooth_elevation(elevation, sigma):
    """
    Gaussian-smooth an elevation raster to remove per-column min/max spike noise.

    NODATA_VALUE pixels are excluded from the convolution (via normalized
    convolution: the data and a validity mask are both smoothed, and the result
    is their ratio) so that gaps neither drag down neighboring elevations nor
    get smoothed into real data. NODATA_VALUE pixels themselves are left
    untouched; gap filling remains the job of `rasterio.fill.fillnodata`
    downstream (e.g. in occpy.util.normalize_occlusion_output).

    Parameters
    ----------
    elevation : numpy.ndarray
        2D elevation raster, with missing columns set to NODATA_VALUE.
    sigma : float
        Standard deviation of the Gaussian kernel, in pixels (= voxels).

    Returns
    -------
    numpy.ndarray
        Smoothed elevation raster, same shape as `elevation`.
    """
    if not sigma:
        return elevation

    valid = elevation != NODATA_VALUE
    if not np.any(valid):
        return elevation

    data = np.where(valid, elevation, 0.0)
    weights = valid.astype(np.float64)

    data_smooth = gaussian_filter(data, sigma=sigma)
    weight_smooth = gaussian_filter(weights, sigma=sigma)

    smoothed = elevation.copy()
    smoothed[valid] = (data_smooth[valid] / weight_smooth[valid])
    return smoothed

def _write_elevation_geotiff(elevation, plot_dim, vox_dim, out_file, crs=None):
    """Write a 2D north-up elevation array to a single-band GeoTIFF."""
    min_x, _, _, _, max_y, _ = plot_dim
    transform = from_origin(min_x, max_y, vox_dim, vox_dim)

    with rasterio.open(
        out_file, "w",
        driver="GTiff",
        height=elevation.shape[0],
        width=elevation.shape[1],
        count=1,
        dtype="float32",
        crs=crs,
        transform=transform,
        nodata=NODATA_VALUE,
    ) as dst:
        dst.write(elevation.astype("float32"), 1)

def get_dtm_from_nhit_grid(nhit_grid, config, out_file, smoothing_sigma=None):
    """
    Generate a DTM GeoTIFF from an Nhit grid by taking, for each (x, y) column,
    the height of the lowest voxel with at least one hit (i.e. the ground return).

    The resulting file has the same extent and pixel size as the voxel grid, so it
    can be used directly as the `dtm_file` argument of
    `occpy.util.normalize_occlusion_output`.

    Parameters
    ----------
    nhit_grid : numpy.ndarray
        3D array (nx, ny, nz) with the number of hits per voxel, as produced by
        the occpy raytracer.
    config : dict
        Configuration dictionary containing 'plot_dim' ([minX, minY, minZ, maxX,
        maxY, maxZ]) and 'vox_dim' (voxel size), matching the values used to
        derive `nhit_grid`. An optional 'crs' key (anything accepted by
        rasterio, e.g. an EPSG string) sets the output CRS.
    out_file : str
        Path to the output GeoTIFF file for the generated DTM.
    smoothing_sigma : float, optional
        Standard deviation of the Gaussian kernel (in voxels/pixels) for a
        nodata-aware smoothing pass to remove per-column max-height spike noise
        before writing. If None, no smoothing is applied.

    Returns
    -------
    numpy.ndarray
        2D DTM array as written to `out_file`, with shape (ny, nx).
    """
    plot_dim = config["plot_dim"]
    vox_dim = config["vox_dim"]

    dtm = _elevation_from_nhit_grid(nhit_grid, plot_dim, vox_dim, use_min=True)
    dtm = _smooth_elevation(dtm, smoothing_sigma)
    _write_elevation_geotiff(dtm, plot_dim, vox_dim, out_file, crs=config.get("crs"))
    return dtm

def get_dsm_from_nhit_grid(nhit_grid, config, out_file, smoothing_sigma=None):
    """
    Generate a canopy-surface GeoTIFF from an Nhit grid by taking, for each (x, y)
    column, the height of the highest voxel with at least one hit.

    The resulting file has the same extent and pixel size as the voxel grid, so it
    can be used directly as the `dsm_file` argument of
    `occpy.util.normalize_occlusion_output` (which derives a canopy height model
    internally as dsm - dtm).

    Parameters
    ----------
    nhit_grid : numpy.ndarray
        3D array (nx, ny, nz) with the number of hits per voxel, as produced by
        the occpy raytracer.
    config : dict
        Configuration dictionary containing 'plot_dim' ([minX, minY, minZ, maxX,
        maxY, maxZ]) and 'vox_dim' (voxel size), matching the values used to
        derive `nhit_grid`. An optional 'crs' key (anything accepted by
        rasterio, e.g. an EPSG string) sets the output CRS. 
    out_file : str
        Path to the output GeoTIFF file for the generated DSM.
    smoothing_sigma : float, optional
        Standard deviation of the Gaussian kernel (in voxels/pixels) for a
        nodata-aware smoothing pass to remove per-column max-height spike noise
        before writing. If None, no smoothing is applied.

    Returns
    -------
    numpy.ndarray
        2D DSM array as written to `out_file`, with shape (ny, nx).
    """
    plot_dim = config["plot_dim"]
    vox_dim = config["vox_dim"]

    dsm = _elevation_from_nhit_grid(nhit_grid, plot_dim, vox_dim, use_min=False)
    dsm = _smooth_elevation(dsm, smoothing_sigma)
    _write_elevation_geotiff(dsm, plot_dim, vox_dim, out_file, crs=config.get("crs"))
    return dsm


class TerrainModel():
    """
    DTM helper built on rasterio.

    Wraps GeoTIFF-backed elevation model and provides convenience methods for
    querying pixel values, retrieving the spatial extent, and cropping/resampling.
    """

    def __init__(self, path2geotiff):
        """
        Open a DTM GeoTIFF for reading.

        Parameters
        ----------
        path2geotiff : str
            Path to a GeoTIFF file (typically single-band) in the dataset CRS.
        """
        self.path2tiff = path2geotiff
        self.dtm = rasterio.open(path2geotiff)


    def get_terrainmodel_path(self):
        """
        Return the path to the currently opened DTM file.

        Returns
        -------
        str
            Absolute path of DTM file.
        """
        return self.dtm.name

    def get_pixel_value(self, x, y):
        """
        Sample the DTM at map coordinates.

        Parameters
        ----------
        x : float
            X coordinate in the dataset CRS.
        y : float
            Y coordinate in the dataset CRS.

        Returns
        -------
        float
            Elevation value at (x, y).
        """
        pix_val = self.dtm.sample([(x,y)])
        pix_val = list(pix_val)

        return pix_val[0][0] #TODO: check if this solution is generic!

    def get_extent(self):
        """
        Return the dataset bounding box.

        Returns
        -------
        rasterio.coords.BoundingBox
            Bounding box (left, bottom, right, top) in the dataset CRS.
        """
        return self.dtm.bounds

    # crop to extent and change resolution if not None
    def crop2extent(self, extent, out_file, res=None):
        """
        Crop the DTM to an extent and optionally resample to a target pixel size.
        Writes the cropped data to `out_file`. 
        If `res` is provided, the raster is resampled using bilinear interpolation 
        to the given pixel size. The instance is updated to reference the new file.

        Parameters
        ----------
        extent : tuple of float
            Bounding box in the dataset CRS as (min_x, max_y, max_x, min_y).
        out_file : str
            Output GeoTIFF path for the cropped (and optionally resampled) raster.
        res : float, optional
            Target pixel size (units per pixel). If None, the original resolution is kept.

        Returns
        -------
        numpy.ndarray
            Cropped (and optionally resampled) data with shape (count, height, width).
        """

        # first crop and then write back
        min_x, max_y, max_x, min_y = extent
        window = rasterio.windows.from_bounds(min_x, min_y, max_x, max_y, self.dtm.transform)
        data = self.dtm.read(1, window=window)
        if data.ndim == 2:
            data = data[np.newaxis, ...]
        # Calculate the transform for the window
        out_transform = self.dtm.window_transform(window)
        profile = self.dtm.profile.copy()
        profile.update({
            "height": data.shape[-2],
            "width": data.shape[-1],
            "transform": out_transform,
            "count": data.shape[0]
        })
        with rasterio.open(out_file, "w", **profile) as dataset:
            dataset.write(data)

        # also do resampling if necessary
        if res!=None:
            with rasterio.open(out_file) as dataset:
                scale_factor_x = dataset.res[0]/res
                scale_factor_y = dataset.res[1]/res

                profile = dataset.profile.copy()
                # resample data to target shape
                data = dataset.read(
                    out_shape=(
                        dataset.count,
                        int(dataset.height * scale_factor_y),
                        int(dataset.width * scale_factor_x)
                    ),
                    resampling=Resampling.bilinear
                )

                # scale image transform
                transform = dataset.transform * dataset.transform.scale(
                    (1 / scale_factor_x),
                    (1 / scale_factor_y)
                )
                profile.update({"height": data.shape[-2],
                                "width": data.shape[-1],
                                "transform": transform})
            with rasterio.open(out_file, "w", **profile) as dataset:
                dataset.write(data)

        # overwrite self.path2tiff to cropped file
        self.path2tiff = out_file
        self.dtm = rasterio.open(out_file)
        return data


