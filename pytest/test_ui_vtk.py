import pytest
import capytaine as cpt

try:
    import vtk
except ImportError:
    vtk = None

try:
    import IPython
except ImportError:
    IPython = None

@pytest.mark.skipif(vtk is None or IPython is None,
                    reason='vtk and/or Ipython are not installed')
def test_animation_of_dofs():
    body = cpt.Sphere()
    body.add_translation_dof(name="Heave")
    animation = body.animate({"Heave": 0.2}, loop_duration=1.0)
    animation.embed_in_notebook()

@pytest.mark.skipif(vtk is None, reason='vtk is not installed')
def test_animate_missing_dof():
    body = cpt.Sphere()
    with pytest.raises(ValueError, match=".*no dof of this name.*"):
        body.animate({"Heave": 0.1}, loop_duration=1.0)
