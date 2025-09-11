import pyg4ometry
import pymcnp


def test_oneBlockVis(vis=False):
    """
    test visualising one block mesh with VTK viewer
    :param vis: visualisation
    :type vis: boolean
    :return: none
    """
    reg = pyg4ometry.mcnp.Registry()
    print(f"block1...")
    b1 = pymcnp.blockphantom.Block("full", rotationSteps=[0, 0, 0], reg=reg)

    if vis:
        v = pyg4ometry.visualisation.VtkViewer()
        v.addAxes()
        v.addMeshSimple(b1.mesh(), name="b1.mesh")

        v.view()


# remove comment when debugging
test_oneBlockVis(True)

if __name__ == "__main__":
    test_oneBlockVis(True)
