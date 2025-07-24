import pyg4ometry
import pymcnp

visulising = 1
def test_twoBlocksVis(vis=False):
    """
    test visualising two blocks mesh with VTK viewer
    :param vis: visualisation
    :type vis: boolean
    :return: none
    """
    print(f"block1...")
    b1 = pymcnp.blockphantom.Block("full", rotationSteps=[0, 1, 1])
    print(f"block2...")
    [b2, c2] = b1.makeNewConnectedBlock("full", 2, 17, makeConnector=True)
    print(f"block3...")
    [b3, c3] = b2.makeNewConnectedBlock("full", 3, 1, makeConnector=True)
    print(f"block4...")
    [b4, c4] = b3.makeNewConnectedBlock("full", 3, 1, makeConnector=True)
    print(f"block5...")
    [b5, c5] = b4.makeNewConnectedBlock("full", 3, 1, makeConnector=True)
    print(f"block6...")
    [b6, c6] = b5.makeNewConnectedBlock("full", 3, 1, makeConnector=True)
    #print(f"rotating...")
    #b2 = b2.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0])

    if vis:
        v = pyg4ometry.visualisation.VtkViewer()
        v.addAxes()
        v.addMeshSimple(b1.mesh(), name="b1.mesh")
        v.addMeshSimple(b2.mesh(), name="b2.mesh")
        v.addMeshSimple(c2.mesh(), name="b2c.mesh")
        v.addMeshSimple(b3.mesh(), name="b3.mesh")
        v.addMeshSimple(c3.mesh(), name="b3c.mesh")
        v.addMeshSimple(b4.mesh(), name="b4.mesh")
        v.addMeshSimple(c4.mesh(), name="b4c.mesh")
        v.addMeshSimple(b5.mesh(), name="b5.mesh")
        v.addMeshSimple(c5.mesh(), name="b5c.mesh")
        v.addMeshSimple(b6.mesh(), name="b6.mesh")
        v.addMeshSimple(c6.mesh(), name="b6c.mesh")

        v.view()


# remove comment when debugging
#test_twoBlocksVis(True)

if __name__ == "__main__":
    test_twoBlocksVis(True)
