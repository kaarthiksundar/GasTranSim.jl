

import json
import numpy as np
import matplotlib.pyplot as plt
import math
from pathlib import Path



dir = "output/solution/"
file_name_pattern = "sol-yamal"

num_rows = len(list(Path(dir).glob("{}*.json".format(file_name_pattern))))
print(num_rows)

file_path = "{}{}-{}.json".format(dir, file_name_pattern, 0)
with open(file_path, 'r') as f:
        data = json.load(f)

num_cols =  len(data["initial_pipe_pressure"]["1"]["distance"])

xt_mat_pressure = np.zeros((num_rows, num_cols))
xt_mat_flux = np.zeros((num_rows, num_cols))
time_hrs = []
for i in range(0, num_rows):
        file_path = "{}{}-{}.json".format(dir, file_name_pattern, i)
        with open(file_path, 'r') as f:
            data = json.load(f)
        pressure = data["initial_pipe_pressure"]["1"]["value"]
        flux = data["initial_pipe_flow"]["1"]["value"]
        time_hrs.append(data["time"]/3600)
        xt_mat_pressure[i] = pressure
        xt_mat_flux[i] = flux



fig, (ax1, ax2) = plt.subplots(figsize=(10, 6), nrows=2)
pos1 = ax1.imshow(xt_mat_pressure*1e-5, origin="lower", cmap="jet")
ax1.set_xlabel("Pressure")
ax1.set_ylabel("Time (hrs)")
ax1.set_xticks(())
ax1.set_yticks(())
fig.colorbar(pos1, ax=ax1, location='right', anchor=(0, 0.3), shrink=1.0)

ax2.plot(time_hrs, 1e-5*xt_mat_pressure[:, -1], "or-")
ax2.set_ylabel("Pressure")
ax2.set_xlabel("Time (hrs)")
ax2.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
plt.show()


# np.savetxt("case-2.txt", xt_mat_pressure)

