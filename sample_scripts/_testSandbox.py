import numpy as np
import pandas as pd
import laspy

is_sorted = lambda a: np.all(a[:-1] <= a[1:])

laz_in = r"D:\_tmp_wdir\OcclusionMappingTests\Rameren\VZ400i\ScanPos001 - SINGLESCANS - 230207_112234 - Not impacted by exit aperture_sort.laz"

test = laspy.read(laz_in)

gps_time = test.gps_time.copy()
return_number = test.return_number.copy()
number_of_returns = test.number_of_returns.copy()

df = pd.DataFrame(data={'gps': gps_time[0:100], 'return_number': return_number[0:100], 'number_of_returns': number_of_returns[0:100]})

is_sorted(gps_time)


