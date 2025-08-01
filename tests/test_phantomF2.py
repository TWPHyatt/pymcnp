import pyg4ometry
import pymcnp


def test_phantomF2(vis=False, write=False):
    """
    test universal whole body phantom F2
    :param vis: visualisation
    :type vis: boolean
    :param write: write to file
    :type write: boolean
    :return: none
    """
    reg = pyg4ometry.mcnp.Registry()

    # crotch
    b1 = pymcnp.blockphantom.Block("full", rotationSteps=[0, 1, 1])

    # left leg
    [b2, b2c1] = b1.makeNewConnectedBlock("full", 2, 22, makeConnector=True)
    b2 = b2.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0])
    [b3, b3c1] = b2.makeNewConnectedBlock("full", 2, 0, makeConnector=True)
    b3 = b3.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0])
    [b4, b4c1] = b3.makeNewConnectedBlock("full", 2, 0, makeConnector=True)
    b4 = b4.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0])
    [b5, b5c1] = b2.makeNewConnectedBlock("half", 13, 13, makeConnector=True)

    # right leg
    [b6, b6c1] = b1.makeNewConnectedBlock("full", 2, 17, makeConnector=True)
    b6 = b6.rotateAboutConnection(hole=2, rotationSteps=[0, 1,   0])
    [b7, b7c1] = b6.makeNewConnectedBlock("full", 2, 0, makeConnector=True)
    b7 = b7.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0])
    [b8, b8c1] = b7.makeNewConnectedBlock("full", 2, 0, makeConnector=True)
    b8 = b8.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0])
    #[b9, b9c1] = b6.makeNewConnectedBlock("half", 13, 20, makeConnector=True)

    # lower abdomen
    #[b10, b10c1] = b1.makeNewConnectedBlock("full", 0, 12, makeConnector=True)
    #b10 = b10.rotateAboutConnection(hole=0, rotationSteps=[0, -1, 0])
    #[b11, b11c1] = b10.makeNewConnectedBlock("full", 20, 13, makeConnector=True)
    #[b12, b12c1] = b10.makeNewConnectedBlock("full", 13, 20, makeConnector=True)
    #[b13, b13c1] = b11.makeNewConnectedBlock("full", 21, 11, makeConnector=True)
    #[b14, b14c1] = b12.makeNewConnectedBlock("full", 12, 11, makeConnector=True)

    # torso
    #[b15, b15c1] = b10.makeNewConnectedBlock("full", 14, 3, makeConnector=True)
    #b15 = b15.rotateAboutConnection(hole=14, rotationSteps=[0, -1, 0])
    #[b16, b16c1] = b15.makeNewConnectedBlock("full", 14, 21, makeConnector=True)
    #b16 = b16.rotateAboutConnection(hole=14, rotationSteps=[0, -1, 0])
    #[b17, b17c1] = b10.makeNewConnectedBlock("full", 0, 14, makeConnector=True)

    # TODO add materials, importance, void & world cells and geometries, etc... for writer

    if write:
        f = pyg4ometry.mcnp.Writer(columnMax=75)
        f.setTitle("PHANTOM F2")
        f.addGeometry(reg=reg)
        f.write("i-phantomF2.txt")

    if vis:
        v = pyg4ometry.visualisation.VtkViewer()
        v.addAxes()
        v.addMeshSimple(b1.mesh(), name="mesh-b1")
        print(f"1 done")
        v.addMeshSimple(b2.mesh(), name="mesh-b2")
        v.addMeshSimple(b2c1.mesh(), name="mesh-b2c1")
        print(f"2 done")
        v.addMeshSimple(b3.mesh(), name="mesh-b3")
        v.addMeshSimple(b3c1.mesh(), name="mesh-b3c1")
        print(f"3 done")
        v.addMeshSimple(b4.mesh(), name="mesh-b4")
        v.addMeshSimple(b4c1.mesh(), name="mesh-b4c1")
        print(f"4 done")
        v.addMeshSimple(b5.mesh(), name="mesh-b5")
        v.addMeshSimple(b5c1.mesh(), name="mesh-b5c1")
        print(f"5 done")
        v.addMeshSimple(b6.mesh(), name="mesh-b6")
        v.addMeshSimple(b6c1.mesh(), name="mesh-b6c1")
        print(f"6 done")
        v.addMeshSimple(b7.mesh(), name="mesh-b7")
        v.addMeshSimple(b7c1.mesh(), name="mesh-b7c1")
        print(f"7 done")
        v.addMeshSimple(b8.mesh(), name="mesh-b8")
        v.addMeshSimple(b8c1.mesh(), name="mesh-b8c1")
        print(f"8 done")
        #v.addMeshSimple(b9.mesh(), name="meshb9")
        #v.addMeshSimple(b9c1.mesh(), name="meshb9c")
        #print(f"9 done")
        #v.addMeshSimple(b10.mesh(), name="meshb10")
        #v.addMeshSimple(b10c1.mesh(), name="meshb10c")
        #print(f"10 done")

        v.view()


# remove comment when debugging
test_phantomF2(vis=True, write=False)

if __name__ == "__main__":
    test_phantomF2(True)
