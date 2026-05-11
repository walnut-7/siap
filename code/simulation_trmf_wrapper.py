'''
Input: 
- pdt
- repl
# --------- TRMF model implementation ---------
- rank
- time_lags
- lambda_w
- lambda_x
- lambda_theta
'''




# --------- config ---------
import sys

pdt = float(sys.argv[1]) # missingness prob
repl = int(sys.argv[2]) # replicate number

# --------- read data ----------
import datetime as dt
from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np

def ymd(value: int) -> int:
  return datetime.strptime(str(value), "%Y%m%d")

nc_data = Dataset("./data/ssi_synthetic_interpolated.nc", mode="r")

# print("Dimensions:")
# for name, dim in nc_data.dimensions.items():
#   print(f"  {name}: size {len(dim)}")

# print("Variables:")
# for name, var in nc_data.variables.items():
#   print(f"  {name}: shape {var.shape}, dtype {var.dtype}")


synthetic_range = [(ymd(20180314) - ymd(18491230)).days, (ymd(20230129) - ymd(18491230)).days + 1]
start = synthetic_range[0]
end = synthetic_range[1]
xt = np.asarray(nc_data.variables['ssi_interpolated'][start:end, :]).T
time = np.asarray(nc_data.variables['time'][:])
wvl = np.asarray(nc_data.variables['wavelength'][:])
t = [timedelta(seconds = (float(minute) + 1.5)*24*60*60) + ymd(18491230) for minute in time]
# t_str = [dt.strftime("%Y-%m-%d") for dt in t]
nc_data.close()

# --------- add missingness ------------
import json
from pathlib import Path

def get_index_col_major(rows, cols, index):
  if not 0 <= index < rows * cols:
    raise IndexError("index out of range")
  col, row = divmod(index, rows)
  return (row, col)

json_root = Path("./code/repl")
mask_path = json_root / f"mask_pdt{int(pdt*10)}_repl{repl}.json"
with mask_path.open() as fh:
  mask0 = json.load(fh)
mask_path.unlink(missing_ok=True)

x = xt.copy()
d1, d2 = x.shape
for index in mask0[0]['miss']:
  row, col = get_index_col_major(d1, d2, index)
  x[row, col] = 0 # np.nan
inds_cal = tuple(np.array(mask0[0]['cp']['S.cal.o']).T)
x[inds_cal] = 0
xt[inds_cal] = 0
x[:,mask0[0]['cp']['S.cal.w']] = 0
xt[:,mask0[0]['cp']['S.cal.w']] = 0

# -------- model --------
import sys
from pathlib import Path
MODULE_ROOT = Path(__file__).resolve().parent
if str(MODULE_ROOT) not in sys.path:
    sys.path.insert(0, str(MODULE_ROOT))

## Reload module to reflect changes during testing
import importlib
import trmf
importlib.reload(trmf)

from trmf import TRMF
from simulation_results import save_result

## Structural hyperparameters
rank = 10
time_lags = np.array([1, 2])
d = time_lags.shape[0]

## Initialize parameters
np.random.seed(int(pdt*10)*100 + repl)
W = 0.1 * np.random.rand(d1, rank)
X = 0.1 * np.random.rand(d2, rank)
theta = 0.1 * np.random.rand(d, rank)
init_para = {"W": W, "X": X, "theta": theta}

## Set hyparameters
lambda_w = 3 # lambda_1 in siap
lambda_x = 3 # alpha in siap
lambda_theta = 5
eta = 1 # lambda_x * eta = lambda_2 in siap
init_hyper = {"lambda_w": lambda_w, "lambda_x": lambda_x, "lambda_theta": lambda_theta, "eta": eta}
maxiter = 200

## Fit model
fit = TRMF(xt, x, init_para, init_hyper, time_lags, maxiter, mute = True)

# --------- save results ----------
res = {
  "pdt": pdt,
  "repl": repl,
  "trmf": fit
}

# Persist metrics so code/simulation_results.py can collate simulation summaries.
save_path = save_result("trmf", pdt, repl, fit)
print(f"Saved TRMF results to {save_path}")
