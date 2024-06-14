from osgeo import gdal
import rasterio
from rasterio.enums import Resampling
from rasterio.fill import fillnodata
import lidario as lio

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
        gdal.Translate(out_file, self.path2tiff, projWin=extent)

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


    def convert2PC(self):
        translator = lio.Translator("geotiff", "dataframe")
        point_cloud = translator.translate(self.path2tiff)
        return point_cloud
    





