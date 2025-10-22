import pyg4ometry
import pymcnp


def test_oneBlockWrite(write=False):
    """
    test writing one block to MCNP input file
    :param write: write to file
    :type write: boolean
    :return: none
    """
    reg = pyg4ometry.mcnp.Registry()

    # CELLS
    # --- PHANTOM ---
    block1 = pymcnp.blockphantom.Block("full", rotationSteps=[0, 0, 1], reg=reg)

    # --- WORLD ---
    cWorld = pyg4ometry.mcnp.Cell(reg=reg)
    cVoid = pyg4ometry.mcnp.Cell(reg=reg)

    # SURFACES
    sSO1 = pyg4ometry.mcnp.SO(25, reg=reg)
    cWorld.addSurface(block1.geometry)  # cell2 is between block1 /
    cWorld.addSurface(sSO1)  # / and s1
    cVoid.addSurface(sSO1)  # cell3 is outside s1

    # MATERIAL
    m2 = pyg4ometry.mcnp.Material(2, -0.001225, reg=reg)
    m0 = pyg4ometry.mcnp.Material(0, reg=reg)  # vacuum
    cWorld.addMaterial(m2)
    cVoid.addMaterial(m0)

    # GEOMETRY
    geoOut = sSO1
    geoIn = pyg4ometry.mcnp.Complement(sSO1)
    geoWorld = pyg4ometry.mcnp.Intersection(geoIn, pyg4ometry.mcnp.Complement(block1))
    cWorld.addGeometry(geoWorld)
    cVoid.addGeometry(geoOut)

    # IMPORTANCE
    i0 = pyg4ometry.mcnp.IMP("p", 0)
    i1 = pyg4ometry.mcnp.IMP("p", 1)
    block1.addImportance(i1)
    cWorld.addImportance(i1)
    cVoid.addImportance(i0)

    if write:
        f = pyg4ometry.mcnp.Writer(columnMax=60)
        f.setTitle("SINGLE ROTATED ISOLATED BLOCK")
        f.addGeometry(reg=reg)
        f.write("i-oneBlockRot.txt")


# remove comment when debugging
#test_oneBlockRotWrite(True)

if __name__ == "__main__":
    test_oneBlockWrite(True)
