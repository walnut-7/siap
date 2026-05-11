import numpy as np
import time


def ten2mat(tensor, mode):
    return np.reshape(np.moveaxis(tensor, mode, 0), (tensor.shape[mode], -1), order = 'F')

def mat2ten(mat, tensor_size, mode):
    index = list()
    index.append(mode)
    for i in range(tensor_size.shape[0]):
        if i != mode:
            index.append(i)
    return np.moveaxis(np.reshape(mat, tensor_size[index].tolist(), order = 'F'), 0, mode)

def svt_tnn(mat, tau, theta):
    [m, n] = mat.shape
    if 2 * m < n:
        u, s, v = np.linalg.svd(mat @ mat.T, full_matrices = 0)
        s = np.sqrt(s)
        idx = np.sum(s > tau)
        mid = np.zeros(idx)
        mid[: theta] = 1
        mid[theta : idx] = (s[theta : idx] - tau) / s[theta : idx]
        return (u[:, : idx] @ np.diag(mid)) @ (u[:, : idx].T @ mat)
    elif m > 2 * n:
        return svt_tnn(mat.T, tau, theta).T
    u, s, v = np.linalg.svd(mat, full_matrices = 0)
    idx = np.sum(s > tau)
    vec = s[: idx].copy()
    vec[theta : idx] = s[theta : idx] - tau
    return u[:, : idx] @ np.diag(vec) @ v[: idx, :]

def compute_mape(var, var_hat):
    return np.sum(np.abs(var - var_hat) / var) / var.shape[0]

def compute_rmse(var, var_hat):
    return  np.sqrt(np.sum((var - var_hat) ** 2) / var.shape[0])

def print_result(it, tol, var, var_hat):
    print('Iter: {}'.format(it))
    print('Tolerance: {:.6}'.format(tol))
    print('Imputation MAPE: {:.6}'.format(compute_mape(var, var_hat)))
    print('Imputation RMSE: {:.6}'.format(compute_rmse(var, var_hat)))
    print()

from scipy import sparse
from scipy.sparse.linalg import spsolve as spsolve

def generate_Psi(dim_time, time_lags):
    Psis = []
    max_lag = np.max(time_lags)
    for i in range(len(time_lags) + 1):
        row = np.arange(0, dim_time - max_lag)
        if i == 0:
            col = np.arange(0, dim_time - max_lag) + max_lag
        else:
            col = np.arange(0, dim_time - max_lag) + max_lag - time_lags[i - 1]
        data = np.ones(dim_time - max_lag)
        Psi = sparse.coo_matrix((data, (row, col)), shape = (dim_time - max_lag, dim_time))
        Psis.append(Psi)
    return Psis

def latc(dense_tensor, sparse_tensor, time_lags, alpha, rho0, lambda0, theta, 
         pad_cols = 0, standardize = False,
         epsilon = 1e-4, maxiter = 100, K = 3, mute = False) -> dict:
    """Low-Rank Autoregressive Tensor Completion (LATC)"""

    start = time.time()

    if standardize:
        mu = sparse_tensor.mean(axis=(1,2), where=(sparse_tensor!=0), dtype=np.float64)
        sd = sparse_tensor.std(axis=(1,2), where=(sparse_tensor!=0), dtype=np.float64)
        sparse_miss = np.where(sparse_tensor[:,:,:] == 0)
        dense_miss = np.where(dense_tensor[:,:,:] == 0)
        sparse_tensor = (sparse_tensor[:,:,:] - mu[:, np.newaxis, np.newaxis]) / sd[:, np.newaxis, np.newaxis]
        dense_tensor = (dense_tensor[:,:,:] - mu[:, np.newaxis, np.newaxis]) / sd[:, np.newaxis, np.newaxis]
        sparse_tensor[sparse_miss] = 0
        dense_tensor[dense_miss] = 0

    dim = np.array(sparse_tensor.shape)
    dim_time = int(np.prod(dim) / dim[0])
    d = len(time_lags)
    max_lag = np.max(time_lags)
    sparse_mat = ten2mat(sparse_tensor, 0)
    pos_missing = np.where(sparse_mat == 0)
    pos_test = np.where((dense_tensor != 0) & (sparse_tensor == 0)) # padding positions are ignored
    dense_test = dense_tensor[pos_test]
    # del dense_tensor

    empty_cols = np.where(np.array([(col == 0).all() for col in sparse_mat.T]))[0]
    if pad_cols > 0:
        empty_cols = empty_cols[:-pad_cols]
    # empty_cols_tensor = [[],[]]
    # for col in empty_cols:
    #     day_idx = int(col // dim[1])
    #     time_idx = int(col % dim[1])
    #     empty_cols_tensor[0].append(time_idx)
    #     empty_cols_tensor[1].append(day_idx)
    nonempty_cols = np.where(np.array([(col != 0).any() for col in sparse_mat.T]))[0]
    # nonempty_cols_tensor = [[],[]]
    # for col in nonempty_cols:
    #     day_idx = int(col // dim[1])
    #     time_idx = int(col % dim[1])
    #     nonempty_cols_tensor[0].append(time_idx)
    #     nonempty_cols_tensor[1].append(day_idx)
    # empty_cols_tensor = np.array(empty_cols_tensor)
    # nonempty_cols_tensor = np.array(nonempty_cols_tensor)

    pos_keys = pos_test[1] + pos_test[2] * dim[1]

    mask_w = np.isin(pos_keys, empty_cols)
    pos_test_w = (pos_test[0][mask_w], 
                  pos_test[1][mask_w],
                  pos_test[2][mask_w])
    mask_o = np.isin(pos_keys, nonempty_cols)
    pos_test_o = (pos_test[0][mask_o], 
                  pos_test[1][mask_o],
                  pos_test[2][mask_o])

    
    T = np.zeros(dim)
    Z_tensor = sparse_tensor.copy()
    Z = sparse_mat.copy()
    A = 0.001 * np.random.rand(dim[0], d)
    Psis = generate_Psi(dim_time, time_lags)
    iden = sparse.coo_matrix((np.ones(dim_time), (np.arange(0, dim_time), np.arange(0, dim_time))), 
                             shape = (dim_time, dim_time))
    it = 0
    ind = np.zeros((d, dim_time - max_lag), dtype = np.int_)
    for i in range(d):
        ind[i, :] = np.arange(max_lag - time_lags[i], dim_time - time_lags[i])
    last_mat = sparse_mat.copy()
    snorm = np.linalg.norm(sparse_mat, 'fro')
    rho = rho0
    while True:
        temp = []
        for m in range(dim[0]):
            Psis0 = Psis.copy()
            for i in range(d):
                Psis0[i + 1] = A[m, i] * Psis[i + 1]
            B = Psis0[0] - sum(Psis0[1 :])
            temp.append(B.T @ B)
        for k in range(K):
            rho = min(rho * 1.05, 1e5)
            tensor_hat = np.zeros(dim)
            for p in range(len(dim)):
                tensor_hat += alpha[p] * mat2ten(svt_tnn(ten2mat(Z_tensor - T / rho, p), 
                                                         alpha[p] / rho, theta), dim, p)
            temp0 = rho / lambda0 * ten2mat(tensor_hat + T / rho, 0)
            mat = np.zeros((dim[0], dim_time))
            for m in range(dim[0]):
                mat[m, :] = spsolve(temp[m] + rho * iden / lambda0, temp0[m, :])
            Z[pos_missing] = mat[pos_missing]
            Z_tensor = mat2ten(Z, dim, 0)
            T = T + rho * (tensor_hat - Z_tensor)
        for m in range(dim[0]):
            A[m, :] = np.linalg.lstsq(Z[m, ind].T, Z[m, max_lag :], rcond = None)[0]
        mat_hat = ten2mat(tensor_hat, 0)
        tol = np.linalg.norm((mat_hat - last_mat), 'fro') / snorm
        last_mat = mat_hat.copy()
        it += 1
        if it % 200 == 0:
            print_result(it, tol, dense_test, tensor_hat[pos_test])
        if (tol < epsilon) or (it >= maxiter):
            break
    if not mute:
        print_result(it, tol, dense_test, tensor_hat[pos_test])

    if standardize:
        tensor_hat = tensor_hat * sd[:, np.newaxis, np.newaxis] + mu[:, np.newaxis, np.newaxis]
        dense_tensor = dense_tensor * sd[:, np.newaxis, np.newaxis] + mu[:, np.newaxis, np.newaxis]
        sparse_tensor = sparse_tensor * sd[:, np.newaxis, np.newaxis] + mu[:, np.newaxis, np.newaxis]
        dense_tensor[dense_miss] = 0
        sparse_tensor[sparse_miss] = 0
        del sparse_miss, dense_miss

    mape = compute_mape(dense_tensor[pos_test], tensor_hat[pos_test])
    rmse = compute_rmse(dense_tensor[pos_test], tensor_hat[pos_test])
    mape_w = compute_mape(dense_tensor[pos_test_w], tensor_hat[pos_test_w])
    mape_o = compute_mape(dense_tensor[pos_test_o], tensor_hat[pos_test_o])

    end = time.time()
    res = {
        "tensor_hat": tensor_hat,
        "time_lags": time_lags.tolist(),
        "iter": it,
        "hyperparam": {
            "alpha": alpha,
            "rho": rho,
            "lambda0": lambda0,
            "theta": theta,
            "epsilon": epsilon,
        },
        "mrae": mape,
        "mrae_o": mape_o,
        "mrae_w": mape_w,
        "rmse": rmse,
        "runtime": end - start
    }
    return res
