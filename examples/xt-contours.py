

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

for i in range(0, num_rows):
        file_path = "{}{}-{}.json".format(dir, file_name_pattern, i)
        with open(file_path, 'r') as f:
            data = json.load(f)
        pressure = data["initial_pipe_pressure"]["1"]["value"]
        flux = data["initial_pipe_flow"]["1"]["value"]
        xt_mat_pressure[i][:] = pressure
        xt_mat_flux[i][:] = flux

# print(xt_mat_flux[0], xt_mat_pressure[0])

fig, (ax1, ax2) = plt.subplots(figsize=(10, 6), nrows=2)
pos1 = ax1.imshow(xt_mat_flux, origin="lower", cmap="RdBu")
pos2 = ax2.imshow(xt_mat_pressure * 1e-5, origin="lower", cmap="RdBu")
ax1.set_xlabel("Flux")
ax1.set_ylabel("Time (hrs)")
ax1.set_xticks(())
ax1.set_yticks(())
# ax1.set_ylim([0, 12])

ax2.set_xlabel("Pressure")
ax2.set_ylabel("Time (hrs)")
ax2.set_xticks(())
ax2.set_yticks(())
# ax1.set_ylim([0, 12])

fig.colorbar(pos1, ax=ax1, location='right', anchor=(0, 0.3), shrink=0.7)
fig.colorbar(pos2, ax=ax2, location='right', anchor=(0, 0.3), shrink=0.7)
plt.show()


