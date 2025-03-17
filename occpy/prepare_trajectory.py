import numpy as np
import pandas as pd
from scipy import interpolate


def interpolate_traj(traj_time, traj_x, traj_y, traj_z, pts_gpstime):


    f_x = interpolate.interp1d(traj_time, traj_x, kind='linear', fill_value="extrapolate")
    sensor_x = f_x(pts_gpstime)
    f_y = interpolate.interp1d(traj_time, traj_y, kind='linear', fill_value="extrapolate")
    sensor_y = f_y(pts_gpstime)
    f_z = interpolate.interp1d(traj_time, traj_z, kind='linear', fill_value="extrapolate")
    sensor_z = f_z(pts_gpstime)

    d = {'time': pts_gpstime, 'sensor_x': sensor_x, 'sensor_y': sensor_y, 'sensor_z': sensor_z}

    df = pd.DataFrame(data=d)

    return df


