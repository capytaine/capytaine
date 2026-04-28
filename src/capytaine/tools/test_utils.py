
import numpy as np

from capytaine.bem.engines import MatrixEngine


def random_complex_array(rng, shape):
    return rng.random(shape) + 1j * rng.random(shape)

class TestEngine(MatrixEngine):
    """Engine that does not compute any physically meaningful matrix, but is very fast and thus can be used for testing the rest of the code (solvers, postprocessing, etc.) without having to wait for the computation of the matrices."""

    def __init__(self, random_seed=9876):
        self.RNG = np.random.default_rng(random_seed)
        self.exportable_settings = {'engine': 'TestEngine'}

    def build_matrices(self, mesh1, mesh2, **kwargs):
        return (
                random_complex_array(self.RNG, (mesh1.nb_faces, mesh2.nb_faces)),
                np.eye(mesh1.nb_faces, mesh2.nb_faces) + random_complex_array(self.RNG, (mesh1.nb_faces, mesh2.nb_faces)),
                )

    def build_S_matrix(self, mesh1, mesh2, **kwargs):
        return random_complex_array(self.RNG, (mesh1.nb_faces, mesh2.nb_faces)),

    def build_fullK_matrix(self, mesh1, mesh2, **kwargs):
        return random_complex_array(self.RNG, (3, mesh1.nb_faces, mesh2.nb_faces)),

    def linear_solver(self, A: np.ndarray, b: np.ndarray) -> np.ndarray:
        return random_complex_array(self.RNG, (A.shape[1],))

    def compute_ram_estimation(self, problem) -> float:
        return 1e6
