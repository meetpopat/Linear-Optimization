import numpy as np


def karmarkar_algorithm(c, A, epsilon=1e-6, max_iter=100):
    m, n = A.shape
    X = np.ones(n) / n
    X = X.T
    Y = X
    Y0 = X
    r = np.sqrt(1 / (n * (n - 1)))
    alpha = (n - 1) / (3 * n)
    k = 0
    z = 1
    while k < max_iter and z > epsilon:
        D = np.diag(X)
        # print("D")
        # print(D)

        P = A.dot(D)
        # print("P")
        # print(P)

        temp = np.ones(n)
        P = np.vstack([P, temp])
        cd = c.dot(D)
        # print("CD")
        # print(cd)

        pp = np.linalg.inv(np.dot(P, P.T))
        q = P.T @ pp @ P
        # print("q")
        # print(q)

        Cp = np.dot((np.eye(n) - q), cd.T)
        # print("Cp")
        # print(Cp)

        cpnorm = np.linalg.norm(Cp, ord=2)
        # print("cpnorm")
        # print(cpnorm)

        Ynew = Y0 - alpha * r * (Cp / cpnorm)
        # print("Y")
        # print(Ynew)

        Xnew = np.dot(D, Ynew) / np.sum(np.dot(D, Ynew))
        # print("X")
        # print(Xnew)

        X = Xnew
        Y = Ynew
        z = c.dot(X)
        # print("z")
        # print(z)
        k += 1
        # print()
        # print()
    return X


c = np.array([1, -3, 3])  # minimize
A = np.array([[1, -3, 2]])

x = karmarkar_algorithm(c, A)
print("X:", x)
print("optimal Value:", c.dot(x))


"""
GIVEN INPUT:
min z = x1 - 3x2 + 3x3
s.t.
    x1 - 3x2 + 2x3 = 0
    x1 + x2 + x3 = 1
    xi >= 0 for  i = 1,2,3
    
Optimal Ans:
Z = 0
X = (3/4,1/4,0)
"""

"""
SOLVES IN CANONICAL FORM:
min cT.x
s.t. Ax = b
     eT.x = 1
     x ≥ 0
Any General LPP can be converted in this form
"""