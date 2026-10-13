import numpy as np
import pytest

from capytaine.tools.perturbation_expansion import PerturbationExpansion


def test_addition_scalar():
    a = PerturbationExpansion([1, 2, 3])
    b = PerturbationExpansion([4, 5, 6])
    result = a + b
    assert result.order == 2
    assert result.coefficients == (5, 7, 9)


def test_subtraction_and_negation():
    a = PerturbationExpansion([1, 2, 3])
    b = PerturbationExpansion([4, 5, 6])
    assert (a - b).coefficients == (-3, -3, -3)
    assert (-a).coefficients == (-1, -2, -3)
    assert (10 - a).coefficients == (9, -2, -3)


def test_elementwise_multiplication_truncated_cauchy_product():
    a = PerturbationExpansion([1, 2, 3])
    b = PerturbationExpansion([4, 5, 6])
    result = a * b
    assert result.order == 2
    assert result.coefficients == (4, 13, 28)


def test_scalar_broadcasting():
    a = PerturbationExpansion([np.array([1.0, 2.0]), np.array([3.0, 4.0])])
    result = a + 1.0
    np.testing.assert_allclose(result.coefficients[0], [2.0, 3.0])
    np.testing.assert_allclose(result.coefficients[1], [3.0, 4.0])

    result2 = 2.0 * a
    np.testing.assert_allclose(result2.coefficients[0], [2.0, 4.0])
    np.testing.assert_allclose(result2.coefficients[1], [6.0, 8.0])


def test_order_mismatch_raises():
    a = PerturbationExpansion([1, 2])
    b = PerturbationExpansion([1, 2, 3])
    with pytest.raises(ValueError):
        a + b
    with pytest.raises(ValueError):
        a * b


def test_getitem_indexes_into_contained_arrays():
    a = PerturbationExpansion([np.array([1.0, 2.0, 3.0]), np.array([4.0, 5.0, 6.0])])

    element = a[1]
    assert isinstance(element, PerturbationExpansion)
    assert element.order == a.order
    assert element.coefficients[0] == 2.0
    assert element.coefficients[1] == 5.0

    sliced = a[1:]
    np.testing.assert_allclose(sliced.coefficients[0], [2.0, 3.0])
    np.testing.assert_allclose(sliced.coefficients[1], [5.0, 6.0])

    with pytest.raises(IndexError):
        a[10]


def test_len_and_iter_act_on_contained_arrays():
    a = PerturbationExpansion([np.array([1.0, 2.0, 3.0]), np.array([4.0, 5.0, 6.0])])
    assert len(a) == 3

    elements = list(a)
    assert len(elements) == 3
    for i, element in enumerate(elements):
        assert isinstance(element, PerturbationExpansion)
        assert element.coefficients[0] == a.coefficients[0][i]
        assert element.coefficients[1] == a.coefficients[1][i]


def test_matmul_cauchy_product():
    I = np.eye(2)
    A1 = np.array([[0.0, 1.0], [2.0, 0.0]])
    A2 = np.array([[1.0, 0.0], [0.0, 1.0]])
    x = PerturbationExpansion([I, A1, A2])

    B0 = np.array([[2.0, 0.0], [0.0, 2.0]])
    B1 = np.array([[0.0, 0.0], [1.0, 0.0]])
    B2 = np.array([[0.0, 1.0], [0.0, 0.0]])
    y = PerturbationExpansion([B0, B1, B2])

    result = x @ y
    expected0 = I @ B0
    expected1 = I @ B1 + A1 @ B0
    expected2 = I @ B2 + A1 @ B1 + A2 @ B0
    np.testing.assert_allclose(result.coefficients[0], expected0)
    np.testing.assert_allclose(result.coefficients[1], expected1)
    np.testing.assert_allclose(result.coefficients[2], expected2)


def test_evaluate():
    a = PerturbationExpansion([1.0, 2.0, 3.0])
    assert a.evaluate(0.1) == pytest.approx(1.0 + 2.0 * 0.1 + 3.0 * 0.1**2)


def inverse_deformation_cofactor():
    I = np.eye(3)
    nabla_u = np.random.rand(3, 3)
    J = PerturbationExpansion([I, nabla_u])
    A = J.det() * J.inv().T
    expected1 = np.trace(nabla_u) * np.eye(3) - nabla_u.T
    assert np.allclose(A.coefficients[1], expected1)
