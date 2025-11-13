import pyg4ometry
import pymcnp


def test_phantomF2Vis(vis=False, write=False):
    """
    test universal whole body phantom F2
    :param vis: visualisation
    :type vis: boolean
    :param write: write to file
    :type write: boolean
    :return: none
    """
    reg = pyg4ometry.mcnp.Registry()

    # CROTCH
    b1 = pymcnp.blockphantom.Block("full", rotationSteps=[0, 1, 1])

    # LEFT LEG
    [b2, b2c1] = b1.makeNewConnectedBlock("full", 2, 22, makeConnector=True)
    b2 = b2.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0])
    [b3, b3c1] = b2.makeNewConnectedBlock("full", 2, 0, makeConnector=True)
    b3 = b3.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0])
    [b4, b4c1] = b3.makeNewConnectedBlock("full", 2, 0, makeConnector=True)
    b4 = b4.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0])
    [b5, b5c1] = b2.makeNewConnectedBlock("half", 13, 13, makeConnector=True)

    # RIGHT LEG
    [b6, b6c1] = b1.makeNewConnectedBlock("full", 2, 17, makeConnector=True)
    b6 = b6.rotateAboutConnection(hole=2, rotationSteps=[0, 1,   0])
    [b7, b7c1] = b6.makeNewConnectedBlock("full", 2, 0, makeConnector=True)
    b7 = b7.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0])
    [b8, b8c1] = b7.makeNewConnectedBlock("full", 2, 0, makeConnector=True)
    b8 = b8.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0])
    [b9, b9c1] = b6.makeNewConnectedBlock("half", 13, 20, makeConnector=True)

    # LOWER ABDOMEN
    [b10, b10c1] = b1.makeNewConnectedBlock("full", 0, 12, makeConnector=True)
    b10 = b10.rotateAboutConnection(hole=0, rotationSteps=[0, 1, 0])
    [b11, b11c1] = b10.makeNewConnectedBlock("full", 20, 13, makeConnector=True)
    [b12, b12c1] = b10.makeNewConnectedBlock("full", 13, 20, makeConnector=True)
    [b13, b13c1] = b11.makeNewConnectedBlock("half", 11, 12, makeConnector=True)
    [b14, b14c1] = b12.makeNewConnectedBlock("half", 11, 21, makeConnector=True)

    # TORSO
    [b15, b15c1] = b10.makeNewConnectedBlock("full", 14, 3, makeConnector=True)
    b15 = b15.rotateAboutConnection(hole=14, rotationSteps=[0, 1, 0])
    [b16, b16c1] = b15.makeNewConnectedBlock("full", 13, 20, makeConnector=True)
    b16 = b16.rotateAboutConnection(hole=13, rotationSteps=[0, 1, 0])
    [b17, b17c1] = b16.makeNewConnectedBlock("full", 0, 19, makeConnector=True)
    b17 = b17.rotateAboutConnection(hole=0, rotationSteps=[0, 1, 0])
    [b18, b18c1] = b17.makeNewConnectedBlock("full", 20, 13, makeConnector=True)
    [b19, b19c1] = b17.makeNewConnectedBlock("full", 13, 20, makeConnector=True)

    # LEFT ARM DOWN
    #[b20, b20c1] = b18.makeNewConnectedBlock("full", 20, 13, makeConnector=True)
    #[b21, b21c1] = b20.makeNewConnectedBlock("full", 2, 0, makeConnector=True)
    #b21 = b21.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0])

    # LEFT ARM UP
    #[b20, b20c1] = b18.makeNewConnectedBlock("full", 20, 13, makeConnector=True)
    #b20 = b20.rotateAboutConnection(hole=20, rotationSteps=[-1, 0, 0])
    #[b21, b21c1] = b20.makeNewConnectedBlock("full", 2, 0, makeConnector=True)
    #b21 = b21.rotateAboutConnection(hole=2, rotationSteps=[0, 0, -1])

    # LEFT ARM 50% UP
    [b20, b20c1] = b18.makeNewConnectedBlock("full", 20, 13, makeConnector=True)
    b20 = b20.rotateAboutConnectionPartial(hole=20, rotationSteps=[-0.5, 0, 0])
    [b21, b21c1] = b20.makeNewConnectedBlock("full", 2, 0, makeConnector=True)
    b21 = b21.rotateAboutConnectionPartial(hole=2, rotationSteps=[0, -1, 0])

    # RIGHT ARM
    [b22, b22c1] = b19.makeNewConnectedBlock("full", 13, 20, makeConnector=True)
    [b23, b23c1] = b22.makeNewConnectedBlock("full", 2, 0, makeConnector=True)
    b23 = b23.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0])

    # NECK AND HEAD
    [b24, b24c1] = b17.makeNewConnectedBlock("full", 0, 2, makeConnector=True)
    b24 = b24.rotateAboutConnection(hole=0, rotationSteps=[0, 1, 0])
    [b25, b25c1] = b24.makeNewConnectedBlock("full", 20, 4, makeConnector=True)
    [b26, b26c1] = b24.makeNewConnectedBlock("half", 16, 12, makeConnector=True)
    [b27, b27c1] = b24.makeNewConnectedBlock("half", 15, 19, makeConnector=True)

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
        #v.addMeshSimple(b2c1.mesh(), name="mesh-b2c1")
        print(f"2 done")
        v.addMeshSimple(b3.mesh(), name="mesh-b3")
        #v.addMeshSimple(b3c1.mesh(), name="mesh-b3c1")
        print(f"3 done")
        v.addMeshSimple(b4.mesh(), name="mesh-b4")
        #v.addMeshSimple(b4c1.mesh(), name="mesh-b4c1")
        print(f"4 done")
        v.addMeshSimple(b5.mesh(), name="mesh-b5")
        #v.addMeshSimple(b5c1.mesh(), name="mesh-b5c1")
        print(f"5 done")
        v.addMeshSimple(b6.mesh(), name="mesh-b6")
        #v.addMeshSimple(b6c1.mesh(), name="mesh-b6c1")
        print(f"6 done")
        v.addMeshSimple(b7.mesh(), name="mesh-b7")
        #v.addMeshSimple(b7c1.mesh(), name="mesh-b7c1")
        print(f"7 done")
        v.addMeshSimple(b8.mesh(), name="mesh-b8")
        #v.addMeshSimple(b8c1.mesh(), name="mesh-b8c1")
        print(f"8 done")
        v.addMeshSimple(b9.mesh(), name="mesh-b9")
        #v.addMeshSimple(b9c1.mesh(), name="mesh-b9c1")
        print(f"9 done")
        v.addMeshSimple(b10.mesh(), name="mesh-b10")
        #v.addMeshSimple(b10c1.mesh(), name="mesh-b10c1")
        print(f"10 done")
        v.addMeshSimple(b11.mesh(), name="mesh-b11")
        #v.addMeshSimple(b11c1.mesh(), name="mesh-b11c1")
        print(f"11 done")
        v.addMeshSimple(b12.mesh(), name="mesh-b12")
        #v.addMeshSimple(b12c1.mesh(), name="mesh-b12c1")
        print(f"12 done")
        v.addMeshSimple(b13.mesh(), name="mesh-b13")
        #v.addMeshSimple(b13c1.mesh(), name="mesh-b13c1")
        print(f"13 done")
        v.addMeshSimple(b14.mesh(), name="mesh-b14")
        #v.addMeshSimple(b14c1.mesh(), name="mesh-b14c1")
        print(f"14 done")
        v.addMeshSimple(b15.mesh(), name="mesh-b15")
        #v.addMeshSimple(b15c1.mesh(), name="mesh-b15c1")
        print(f"15 done")
        v.addMeshSimple(b16.mesh(), name="mesh-b16")
        #v.addMeshSimple(b16c1.mesh(), name="mesh-b16c1")
        print(f"16 done")
        v.addMeshSimple(b17.mesh(), name="mesh-b17")
        #v.addMeshSimple(b17c1.mesh(), name="mesh-b17c1")
        print(f"17 done")
        v.addMeshSimple(b18.mesh(), name="mesh-b18")
        #v.addMeshSimple(b18c1.mesh(), name="mesh-b18c1")
        print(f"18 done")
        v.addMeshSimple(b19.mesh(), name="mesh-b19")
        #v.addMeshSimple(b19c1.mesh(), name="mesh-b19c1")
        print(f"19 done")
        v.addMeshSimple(b20.mesh(), name="mesh-b20")
        #v.addMeshSimple(b20c1.mesh(), name="mesh-b20c1")
        print(f"20 done")
        v.addMeshSimple(b21.mesh(), name="mesh-b21")
        #v.addMeshSimple(b21c1.mesh(), name="mesh-b21c1")
        print(f"21 done")
        v.addMeshSimple(b22.mesh(), name="mesh-b22")
        #v.addMeshSimple(b22c1.mesh(), name="mesh-b22c1")
        print(f"22 done")
        v.addMeshSimple(b23.mesh(), name="mesh-b23")
        #v.addMeshSimple(b23c1.mesh(), name="mesh-b23c1")
        print(f"23 done")
        v.addMeshSimple(b24.mesh(), name="mesh-b24")
        #v.addMeshSimple(b24c1.mesh(), name="mesh-b24c1")
        print(f"24 done")
        v.addMeshSimple(b25.mesh(), name="mesh-b25")
        #v.addMeshSimple(b25c1.mesh(), name="mesh-b25c1")
        print(f"25 done")
        v.addMeshSimple(b26.mesh(), name="mesh-b26")
        #v.addMeshSimple(b26c1.mesh(), name="mesh-b27c1")
        print(f"26 done")
        v.addMeshSimple(b27.mesh(), name="mesh-b27")
        #v.addMeshSimple(b27c1.mesh(), name="mesh-b27c1")
        print(f"27 done")

        v.view()


# remove comment when debugging
test_phantomF2Vis(vis=True, write=False)

if __name__ == "__main__":
    test_phantomF2Vis(True)
