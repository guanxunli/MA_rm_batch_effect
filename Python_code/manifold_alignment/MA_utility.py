import numpy as np

def manifold_alignment(Wx, Wy, Wxy = None, mu = 0.9, d = 30, normlized = True, eps = 1e-5):
    Wx = np.array(Wx)
    Wy = np.array(Wy)
    n = Wx.shape[0]

    if normlized == True:
        Wx = Wx/np.max(np.abs(Wx))
        Wy = Wy/np.max(np.abs(Wy))
        # Wx = Wx + 1
        # Wy = Wy + 1

    if Wxy == None:
        Wxy = np.diag(np.repeat(1, n))

    Wxy = mu * (np.sum(Wx) + np.sum(Wy)) / (2 * np.sum(Wxy)) * Wxy
    W = np.asarray(np.bmat(((Wx, Wxy), (Wxy.T, Wy))))
    W = (W + W.T)/2
    D = np.diag(np.sum(W, axis = 0))
    L = D - W

    vals, vecs = np.linalg.eig(L)
    idx = np.argsort(vals)
    for i in range(len(idx)):
        if vals[idx[i]] >= eps:
            break

    alignedNet = vecs.real[:,idx[i:(i + d)]]
    for i in range(alignedNet.shape[1]):
        alignedNet[:,i] /= np.linalg.norm(alignedNet[:,i])

    return alignedNet