import pyg4ometry
import pymcnp

def test_fourBlocksRotVis(vis=False):
    """
    test rotating blocks about connection points
    :param vis: visualisation
    :type vis: boolean
    :return: none
    """
    print(f"block1...")
    b1 = pymcnp.blockphantom.Block("full", rotationSteps=[0, 1, 1])
    print(f"block2...")
    [b2, c2] = b1.makeNewConnectedBlock("full", 2, 22, makeConnector=True)
    print(f"rotating...")
    b2 = b2.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0])
    print(f"block3...")
    [b3, c3] = b2.makeNewConnectedBlock("full", 3, 13, makeConnector=True)
    print(f"rotating...")
    b3 = b3.rotateAboutConnection(hole=3, rotationSteps=[1, 0, 0])
    print(f"block4...")
    [b4, c4] = b3.makeNewConnectedBlock("full", 8, 5, makeConnector=True)

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

        v.view()


# remove comment when debugging
#test_fourBlocksRotVis(True)

if __name__ == "__main__":
    test_fourBlocksRotVis(True)
