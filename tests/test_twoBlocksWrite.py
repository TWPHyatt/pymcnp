import pyg4ometry
import pymcnp


def test_twoBlocksWrite(write=False):
    """
    test writing two blocks to MCNP input file
    :param write: write to file
    :type write: boolean
    :return: none
    """
    reg = pyg4ometry.mcnp.Registry()

    # CELLS
    # --- PHANTOM ---
    block1 = pymcnp.blockphantom.Block("full", rotationSteps=[0, 1, 1], reg=reg)
    [block2, con2] = block1.makeNewConnectedBlock("full", 2, 22, makeConnector=True, reg=reg)
    block2 = block2.rotateAboutConnection(hole=2, rotationSteps=[0, 1, 0], reg=reg)
    # --- WORLD ---
    cWorld = pyg4ometry.mcnp.Cell(reg=reg)
    cVoid = pyg4ometry.mcnp.Cell(reg=reg)

    # SURFACES
    sSO1 = pyg4ometry.mcnp.SO(25, reg=reg)
    cWorld.addSurface(block1.geometry)  # cell2 is between block1 /
    cWorld.addSurface(block2.geometry)  # and block2 /
    cWorld.addSurface(con2.geometry)  # and connector1 /
    cWorld.addSurface(sSO1)  # / and s1
    cVoid.addSurface(sSO1)  # cell3 is outside s1

    # GEOMETRY
    geoOut = sSO1
    geoIn = pyg4ometry.mcnp.Complement(sSO1)
    geo1 = pyg4ometry.mcnp.Intersection(geoIn, pyg4ometry.mcnp.Complement(block1))
    geo2 = pyg4ometry.mcnp.Intersection(geo1, pyg4ometry.mcnp.Complement(block2))
    geoWorld = pyg4ometry.mcnp.Intersection(geo2, pyg4ometry.mcnp.Complement(con2))
    cWorld.addGeometry(geoWorld)
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
    block1.addImportance(i1)
    block2.addImportance(i1)
    con2.addImportance(i1)
    cWorld.addImportance(i1)
    cVoid.addImportance(i0)

    if write:
        f = pyg4ometry.mcnp.Writer(columnMax=60)
        f.setTitle("TWO CONNECTED BLOCKS")
        f.addGeometry(reg=reg)
        f.write("i-twoBlocks.txt")


# remove comment when debugging
test_twoBlocksWrite(True)

if __name__ == "__main__":
    test_twoBlocksWrite(True)
