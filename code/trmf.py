import numpy as np
from numpy.linalg import inv as inv
import time


def TRMF(dense_mat, sparse_mat, init_para, init_hyper, time_lags, maxiter,
        standardize = False, mute = False) -> dict:
    """Temporal Regularized Matrix Factorization, TRMF."""

    start = time.time()

    if standardize:
        mu = sparse_mat.mean(axis=1, where=(sparse_mat != 0), dtype=np.float64)
        sd = sparse_mat.std(axis=1, where=(sparse_mat != 0), dtype=np.float64)
        sparse_miss = np.where(sparse_mat == 0)
        dense_miss = np.where(dense_mat == 0)
        sparse_mat = (sparse_mat - mu[:, np.newaxis]) / sd[:, np.newaxis]
        dense_mat = (dense_mat - mu[:, np.newaxis]) / sd[:, np.newaxis]
        sparse_mat[sparse_miss] = 0
        dense_mat[dense_miss] = 0

    ## Initialize parameters
    W = init_para["W"]
    X = init_para["X"]
    theta = init_para["theta"]
    
    ## Set hyperparameters
    lambda_w = init_hyper["lambda_w"]
    lambda_x = init_hyper["lambda_x"]
    lambda_theta = init_hyper["lambda_theta"]
    eta = init_hyper["eta"]
    
    dim1, dim2 = sparse_mat.shape
    pos_train = np.where(sparse_mat != 0)
    pos_test = np.where((dense_mat != 0) & (sparse_mat == 0))
    empty_cols = np.where(np.array([(col == 0).all() for col in sparse_mat.T]))[0]
    nonempty_cols = np.where(np.array([(col != 0).any() for col in sparse_mat.T]))[0]
    pos_test_w = (pos_test[0][np.isin(pos_test[1], empty_cols)], pos_test[1][np.isin(pos_test[1], empty_cols)])
    pos_test_o = (pos_test[0][np.isin(pos_test[1], nonempty_cols)], pos_test[1][np.isin(pos_test[1], nonempty_cols)])


    binary_mat = sparse_mat.copy()
    binary_mat[pos_train] = 1
    d, rank = theta.shape
    
    for it in range(maxiter):
        ## Update spatial matrix W
        for i in range(dim1):
            pos0 = np.where(sparse_mat[i, :] != 0)
            Xt = X[pos0[0], :]
            vec0 = Xt.T @ sparse_mat[i, pos0[0]]
            mat0 = inv(Xt.T @ Xt + lambda_w * np.eye(rank))
            W[i, :] = mat0 @ vec0
        ## Update temporal matrix X
        for t in range(dim2):
            pos0 = np.where(sparse_mat[:, t] != 0)
            Wt = W[pos0[0], :]
            Mt = np.zeros((rank, rank))
            Nt = np.zeros(rank)
            if t < np.max(time_lags):
                Pt = np.zeros((rank, rank))
                Qt = np.zeros(rank)
            else:
                Pt = np.eye(rank)
                Qt = np.einsum('ij, ij -> j', theta, X[t - time_lags, :])
            if t < dim2 - np.min(time_lags):
                if t >= np.max(time_lags) and t < dim2 - np.max(time_lags):
                    index = list(range(0, d))
                else:
                    index = list(np.where((t + time_lags >= np.max(time_lags)) & (t + time_lags < dim2)))[0]
                for k in index:
                    Ak = theta[k, :]
                    Mt += np.diag(Ak ** 2)
                    theta0 = theta.copy()
                    theta0[k, :] = 0
                    Nt += np.multiply(Ak, X[t + time_lags[k], :]
                                      - np.einsum('ij, ij -> j', theta0, X[t + time_lags[k] - time_lags, :]))
            vec0 = Wt.T @ sparse_mat[pos0[0], t] + lambda_x * Nt + lambda_x * Qt
            mat0 = inv(Wt.T @ Wt + lambda_x * Mt + lambda_x * Pt + lambda_x * eta * np.eye(rank))
            X[t, :] = mat0 @ vec0
        ## Update AR coefficients theta
        for k in range(d):
            theta0 = theta.copy()
            theta0[k, :] = 0
            mat0 = np.zeros((dim2 - np.max(time_lags), rank))
            for L in range(d):
                mat0 += X[np.max(time_lags) - time_lags[L] : dim2 - time_lags[L] , :] @ np.diag(theta0[L, :])
            VarPi = X[np.max(time_lags) : dim2, :] - mat0
            var1 = np.zeros((rank, rank))
            var2 = np.zeros(rank)
            for t in range(np.max(time_lags), dim2):
                B = X[t - time_lags[k], :]
                var1 += np.diag(np.multiply(B, B))
                var2 += np.diag(B) @ VarPi[t - np.max(time_lags), :]
            theta[k, :] = inv(var1 + lambda_theta * np.eye(rank) / lambda_x) @ var2

        mat_hat = W @ X.T
        mape = np.sum(np.abs(dense_mat[pos_test] - mat_hat[pos_test]) 
                      / dense_mat[pos_test]) / dense_mat[pos_test].shape[0]
        mape_o = np.sum(np.abs(dense_mat[pos_test_o] - mat_hat[pos_test_o]) 
                      / dense_mat[pos_test_o]) / dense_mat[pos_test_o].shape[0]
        mape_w = np.sum(np.abs(dense_mat[pos_test_w] - mat_hat[pos_test_w]) 
                      / dense_mat[pos_test_w]) / dense_mat[pos_test_w].shape[0]        
        rmse = np.sqrt(np.sum((dense_mat[pos_test] - mat_hat[pos_test]) ** 2)/dense_mat[pos_test].shape[0])
        
        if (it + 1) % 100 == 0 and not mute:
            print('Iter: {}'.format(it + 1))
            print('Imputation MAPE: {:.6}'.format(mape))
            print('Imputation RMSE: {:.6}'.format(rmse))
            print()

    if standardize:
        mat_hat = mat_hat * sd[:, np.newaxis] + mu[:, np.newaxis]
        dense_mat = dense_mat * sd[:, np.newaxis] + mu[:, np.newaxis]
        sparse_mat = sparse_mat * sd[:, np.newaxis] + mu[:, np.newaxis]
        dense_mat[dense_miss] = 0
        sparse_mat[sparse_miss] = 0
        del sparse_miss, dense_miss
        mape = np.sum(np.abs(dense_mat[pos_test] - mat_hat[pos_test])
                      / dense_mat[pos_test]) / dense_mat[pos_test].shape[0]
        mape_o = np.sum(np.abs(dense_mat[pos_test_o] - mat_hat[pos_test_o])
                        / dense_mat[pos_test_o]) / dense_mat[pos_test_o].shape[0]
        mape_w = np.sum(np.abs(dense_mat[pos_test_w] - mat_hat[pos_test_w])
                        / dense_mat[pos_test_w]) / dense_mat[pos_test_w].shape[0]
        rmse = np.sqrt(np.sum((dense_mat[pos_test] - mat_hat[pos_test]) ** 2)
                       / dense_mat[pos_test].shape[0])

    end = time.time()
    res = {
        "W": W,
        "X": X,
        "theta": theta,
        "rank": rank,
        "time_lags": time_lags.tolist(),
        "iter": it,
        "hyperparam": {
            "lambda_w": lambda_w,
            "lambda_x": lambda_x,
            "lambda_theta": lambda_theta,
            "eta": eta,
        },
        "mrae": mape,
        "mrae_o": mape_o,
        "mrae_w": mape_w,
        "rmse": rmse,
        "runtime": end - start
    }
    return res