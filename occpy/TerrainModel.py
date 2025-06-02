import rasterio
from rasterio.enums import Resampling

class TerrainModel():
    def __init__(self, path2geotiff):
        self.path2tiff = path2geotiff
        self.dtm = rasterio.open(path2geotiff)


    def get_terrainmodel_path(self):
        return self.dtm.name

    def get_pixel_value(self, x, y):
        pix_val = self.dtm.sample([(x,y)])
        pix_val = list(pix_val)

        return pix_val[0][0] #TODO: check if this solution is generic!

    def get_extent(self):
        return self.dtm.bounds

    # crop to extent and change resolution if not None
    def crop2extent(self, extent, out_file, res):
        min_x, max_y, max_x, min_y = extent
        window = rasterio.windows.from_bounds(min_x, min_y, max_x, max_y, self.dtm.transform)

        dataset = self.dtm.read(window)
        # Calculate the transform for the window
        out_transform = dataset.transform * self.dtm.window_transform(window)

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
                out_transform = out_transform * dataset.transform.scale(
                    (1 / scale_factor_x),
                    (1 / scale_factor_y)
                )
                profile.update({"height": data.shape[-2],
                                "width": data.shape[-1],
                                "transform": out_transform})
        else:
            profile = self.dtm.profile.copy()
            profile.update({
                "height": data.shape[-2],
                "width": data.shape[-1],
                "transform": out_transform
            })
        
        with rasterio.open(out_file, "w", **profile) as dataset:
            dataset.write(data)

        # overwrite self.path2tiff to cropped file
        self.path2tiff = out_file
        self.dtm = rasterio.open(out_file)
        return data






