
import numpy as np

from capytaine.matrices.block_matrices import BlockMatrix
from capytaine.matrices.block_toeplitz_matrices import BlockSymmetricToeplitzMatrix, BlockSymmetricCirculantMatrix


def random_block_matrix(list_of_x_shapes, list_of_y_shapes):
    A = []
    for x in list_of_x_shapes:
        line = []
        for y in list_of_y_shapes:
            line.append(np.random.rand(x, y))
        A.append(line)
    return BlockMatrix(A)


def full_like(A, value):
    if isinstance(A, BlockSymmetricCirculantMatrix):
        new_matrix = []
        for i in range(len(A._c_blocks)):
            new_matrix.append(full_like(A._c_blocks[i], value))
        return BlockSymmetricCirculantMatrix(I)
    elif isinstance(A, BlockSymmetricToeplitzMatrix):
        new_matrix = []
        for i in range(A.nb_blocks[0]):
            new_matrix.append(full_like(A._t_blocks[i], value))
        return BlockSymmetricToeplitzMatrix(I)
    elif isinstance(A, BlockMatrix):
        new_matrix = []
        for i in range(A.nb_blocks[0]):
            line = []
            for j in range(A.nb_blocks[1]):
                line.append(full_like(A.blocks[i][j], value))
            new_matrix.append(line)
        return BlockMatrix(new_matrix)
    elif isinstance(A, np.ndarray):
        return np.full_like(A, value)


def zeros_like(A):
    return full_like(A, 0.0)


def ones_like(A):
    return full_like(A, 1.0)


def identity_like(A):
    if isinstance(A, BlockSymmetricCirculantMatrix):
        I = [identity_like(A._c_blocks[0])]
        for i in range(1, len(A._c_blocks)):
            I.append(zeros_like(A._c_blocks[i]))
        return BlockSymmetricCirculantMatrix(I)
    elif isinstance(A, BlockSymmetricToeplitzMatrix):
        I = [identity_like(A._t_blocks[0])]
        for i in range(1, A.nb_blocks[0]):
            I.append(zeros_like(A._t_blocks[i]))
        return BlockSymmetricToeplitzMatrix(I)
    elif isinstance(A, BlockMatrix):
        I = []
        for i in range(A.nb_blocks[0]):
            line = []
            for j in range(A.nb_blocks[1]):
                if i == j:
                    line.append(identity_like(A.blocks[i][j]))
                else:
                    line.append(zeros_like(A.blocks[i][j]))
            I.append(line)
        return BlockMatrix(I)
    elif isinstance(A, np.ndarray):
        return np.eye(A.shape[0], A.shape[1])


if __name__ == '__main__':
    A = random_block_matrix([2, 3, 1], [2, 3, 1])
    A.plot_shape()

    O = zeros_like(A)
    print(repr(O))
    print(O.full_matrix())

    A2 = random_block_matrix([2, 3, 1], [2, 3, 1])
    B = BlockSymmetricToeplitzMatrix([A, A2])
    B.plot_shape()

    I = identity_like(B)
    print(repr(I))
    print(I.full_matrix())

