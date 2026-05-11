'''
Input: 
- pdt
- repl
# --------- LATC model implementation ---------
- time_lags
- alpha
- rho
- lambda0
- theta
- epsilon
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

synthetic_range = [(ymd(20180314) - ymd(18491230)).days, (ymd(20230129) - ymd(18491230)).days + 1]
start = synthetic_range[0]
end = synthetic_range[1]
xt = np.asarray(nc_data.variables['ssi_interpolated'][start:end, :]).T
time = np.asarray(nc_data.variables['time'][:])
wvl = np.asarray(nc_data.variables['wavelength'][:])
t = [timedelta(seconds = (float(minute) + 1.5)*24*60*60) + ymd(18491230) for minute in time]
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

import importlib
import latc
importlib.reload(latc)

from latc import latc
from simulation_results import save_result

def matrix_to_tensor(m: str, period_len: int) -> np.ndarray:
  """
  matrix: variable x time
  tensor: variable x day of period x period number
  """
  mat = m.copy()
  remainder = mat.shape[1] % period_len   # Pad with zeros to make columns a multiple of period_len
  pad_cols = 0
  if remainder != 0:
    pad_cols = period_len - remainder
    mat = np.pad(mat, ((0, 0), (0, pad_cols)), mode='constant', constant_values=0)
  tensor = mat.reshape(mat.shape[0], -1, period_len).transpose(0,2,1).copy()
  return tensor, pad_cols

period = 365
xt_tensor, pad_cols = matrix_to_tensor(xt, period)
x_tensor, _ = matrix_to_tensor(x, period)

time_lags = np.arange(1, 3)
alpha = np.ones(3) * 1 # alpha_0, alpha_1, alpha_2: weights of nuclear norm penalties
rho = 1e-4 # ADMM penalty parameter (not part of the loss)
theta = 10 # number of leading singular values treated as free (unpenalized)
lambda0 = 6 # lambda0 / 2: weight of the autoregressive (AR) penalty
fit = latc(xt_tensor, x_tensor, time_lags, alpha, rho, lambda0, theta, maxiter = 100, pad_cols = pad_cols, mute = True, standardize = True)

# --------- save results ----------
res = {
  "pdt": pdt,
  "repl": repl,
  "latc": fit
}

# Persist metrics so code/simulation_results.py can collate simulation summaries.
save_path = save_result("latc", pdt, repl, fit)
print(f"Saved LATC results to {save_path}")