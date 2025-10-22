import rasterio
import numpy as np
from rasterio.enums import Resampling

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


