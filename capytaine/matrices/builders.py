
import numpy as np

from capytaine.matrices.block_matrices import BlockMatrix
from capytaine.matrices.block_toeplitz_matrices import BlockSymmetricToeplitzMatrix


def random_block_matrix(list_of_x_shapes, list_of_y_shapes):
    A = []
    for x in list_of_x_shapes:
        line = []
        for y in list_of_y_shapes:
            line.append(np.random.rand(x, y))
        A.append(line)
    return BlockMatrix(A)


def zeros_like(A):
    if isinstance(A, BlockSymmetricToeplitzMatrix):
        I = []
        for i in range(A.nb_blocks[0]):
            I.append(zeros_like(A.t_blocks[i]))
        return BlockSymmetricToeplitzMatrix(I)
    elif isinstance(A, BlockMatrix):
        I = []
        for i in range(A.nb_blocks[0]):
            line = []
            for j in range(A.nb_blocks[1]):
                line.append(zeros_like(A.blocks[i][j]))
            I.append(line)
        return BlockMatrix(I)
    elif isinstance(A, np.ndarray):
        return np.zeros_like(A)


def identity_like(A):
    if isinstance(A, BlockSymmetricToeplitzMatrix):
        I = [identity_like(A.t_blocks[0])]
        for i in range(1, A.nb_blocks[0]):
            I.append(zeros_like(A.t_blocks[i]))
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

