import pyg4ometry
import pymcnp


def test_twoBlocksVis(vis=False):
    """
    test visualising two blocks mesh with VTK viewer
    :param vis: visualisation
    :type vis: boolean
    :return: none
    """
    reg = pyg4ometry.mcnp.Registry()
    # CELLS
    # --- PHANTOM ---
    block1 = pymcnp.blockphantom.Block("full", rotationSteps=[0, 1, 1], reg=reg)
    [block2, con2] = block1.makeNewConnectedBlock("full", 2, 22, makeConnector=True, reg=reg)
    block2 = block2.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0], reg=reg)

    if vis:
        v = pyg4ometry.visualisation.VtkViewer()
        v.addAxes()
        v.addMeshSimple(block1.mesh(), name="b1.mesh")
        v.addMeshSimple(block2.mesh(), name="b2.mesh")
        v.addMeshSimple(con2.mesh(), name="b2c.mesh")

        v.view()


# remove comment when debugging
test_twoBlocksVis(True)

if __name__ == "__main__":
    test_twoBlocksVis(True)
