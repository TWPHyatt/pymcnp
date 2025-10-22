import pyg4ometry
import pymcnp


def test_phantomF2Write(write=False):
    """
    test writing whole body phantom F2
    :param write: write to file
    :type write: boolean
    :return: none
    """
    reg = pyg4ometry.mcnp.Registry()

    # CELLS
    # --- PHANTOM ---

    # CROTCH
    b1 = pymcnp.blockphantom.Block("full", rotationSteps=[0, 1, 1], reg=reg)

    # LEFT LEG
    [b2, b2c1] = b1.makeNewConnectedBlock("full", 2, 22, makeConnector=True, reg=reg)
    b2 = b2.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0], reg=reg)
    [b3, b3c1] = b2.makeNewConnectedBlock("full", 2, 0, makeConnector=True, reg=reg)
    b3 = b3.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0], reg=reg)
    [b4, b4c1] = b3.makeNewConnectedBlock("full", 2, 0, makeConnector=True, reg=reg)
    b4 = b4.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0], reg=reg)
    [b5, b5c1] = b2.makeNewConnectedBlock("half", 13, 13, makeConnector=True, reg=reg)

    # RIGHT LEG
    [b6, b6c1] = b1.makeNewConnectedBlock("full", 2, 17, makeConnector=True, reg=reg)
    b6 = b6.rotateAboutConnection(hole=2, rotationSteps=[0, 1,   0], reg=reg)
    [b7, b7c1] = b6.makeNewConnectedBlock("full", 2, 0, makeConnector=True, reg=reg)
    b7 = b7.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0], reg=reg)
    [b8, b8c1] = b7.makeNewConnectedBlock("full", 2, 0, makeConnector=True, reg=reg)
    b8 = b8.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0], reg=reg)
    [b9, b9c1] = b6.makeNewConnectedBlock("half", 13, 20, makeConnector=True, reg=reg)

    # LOWER ABDOMEN
    [b10, b10c1] = b1.makeNewConnectedBlock("full", 0, 12, makeConnector=True, reg=reg)
    b10 = b10.rotateAboutConnection(hole=0, rotationSteps=[0, 1, 0], reg=reg)
    [b11, b11c1] = b10.makeNewConnectedBlock("full", 20, 13, makeConnector=True, reg=reg)
    [b12, b12c1] = b10.makeNewConnectedBlock("full", 13, 20, makeConnector=True, reg=reg)
    [b13, b13c1] = b11.makeNewConnectedBlock("half", 11, 12, makeConnector=True, reg=reg)
    [b14, b14c1] = b12.makeNewConnectedBlock("half", 11, 21, makeConnector=True, reg=reg)

    # TORSO
    [b15, b15c1] = b10.makeNewConnectedBlock("full", 14, 3, makeConnector=True, reg=reg)
    b15 = b15.rotateAboutConnection(hole=14, rotationSteps=[0, 1, 0], reg=reg)
    [b16, b16c1] = b15.makeNewConnectedBlock("full", 13, 20, makeConnector=True, reg=reg)
    b16 = b16.rotateAboutConnection(hole=13, rotationSteps=[0, 1, 0], reg=reg)
    [b17, b17c1] = b16.makeNewConnectedBlock("full", 0, 19, makeConnector=True, reg=reg)
    b17 = b17.rotateAboutConnection(hole=0, rotationSteps=[0, 1, 0], reg=reg)
    [b18, b18c1] = b17.makeNewConnectedBlock("full", 20, 13, makeConnector=True, reg=reg)
    [b19, b19c1] = b17.makeNewConnectedBlock("full", 13, 20, makeConnector=True, reg=reg)

    # --- WORLD ---
    cWorld = pyg4ometry.mcnp.Cell(reg=reg)
    cVoid = pyg4ometry.mcnp.Cell(reg=reg)

    # SURFACES
    blocks = [b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19]
    connectors = [b2c1, b3c1, b4c1, b5c1, b6c1, b7c1, b8c1, b9c1, b10c1, b11c1,
                  b12c1, b13c1, b14c1, b15c1, b16c1, b17c1, b18c1, b19c1]
    sSO1 = pyg4ometry.mcnp.SO(25, reg=reg)
    for b in blocks:
        cWorld.addSurface(b.geometry)
    for c in connectors:
        cWorld.addSurface(c.geometry)
    cWorld.addSurface(sSO1)
    cVoid.addSurface(sSO1)

    # GEOMETRY
    geoOut = sSO1
    geoIn = pyg4ometry.mcnp.Complement(sSO1)

    geoB1 = geoIn
    for b in blocks:
        geoB2 = pyg4ometry.mcnp.Intersection(geoB1, pyg4ometry.mcnp.Complement(b))
        geoB1 = geoB2
    geoC1 = geoB1
    for c in connectors:
        geoC2 = pyg4ometry.mcnp.Intersection(geoC1, pyg4ometry.mcnp.Complement(c))
        geoC1 = geoC2

    cWorld.addGeometry(geoC1)
    cVoid.addGeometry(geoOut)

    # MATERIAL
    m0 = pyg4ometry.mcnp.Material(0, reg=reg)
    # material numbers 1 & 2 are used for the block and connector
    m3 = pyg4ometry.mcnp.Material(3, -0.001225, reg=reg)
    cWorld.addMaterial(m3)
    cVoid.addMaterial(m0)

    # IMPORTANCE
    i0 = pyg4ometry.mcnp.IMP("p", 0)
    i1 = pyg4ometry.mcnp.IMP("p", 1)
    for b in blocks:
        b.addImportance(i1)
    for c in connectors:
        c.addImportance(i1)
    cWorld.addImportance(i1)
    cVoid.addImportance(i0)

    if write:
        f = pyg4ometry.mcnp.Writer(columnMax=60)
        f.setTitle("F2 PHANTOM BLOCKS")
        f.addGeometry(reg=reg)
        f.write("i-F2Phantom.txt")


# remove comment when debugging
test_phantomF2Write(True)

if __name__ == "__main__":
    test_phantomF2Write(True)