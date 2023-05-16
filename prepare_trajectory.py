import numpy as np
import pandas as pd
from scipy import interpolate


def interpolate_traj(traj, pts_gpstime):


    f_x = interpolate.interp1d(traj['%time'], traj['x'], kind='linear', fill_value="extrapolate")
    sensor_x = f_x(pts_gpstime)
    f_y = interpolate.interp1d(traj['%time'], traj['y'], kind='linear', fill_value="extrapolate")
    sensor_y = f_y(pts_gpstime)
    f_z = interpolate.interp1d(traj['%time'], traj['z'], kind='linear', fill_value="extrapolate")
    sensor_z = f_z(pts_gpstime)

    d = {'time': pts_gpstime, 'sensor_x': sensor_x, 'sensor_y': sensor_y, 'sensor_z': sensor_z}

    df = pd.DataFrame(data=d)

    return df


