import pyg4ometry
import pymcnp


def test_smallBoxWithHole(vis=True, write=False):
    reg = pyg4ometry.mcnp.Registry()

    # CELLS
    # --- PHANTOM ---
    cBox = pyg4ometry.mcnp.Cell(reg=reg)

    # --- WORLD ---
    cWorld = pyg4ometry.mcnp.Cell(reg=reg)
    cVoid = pyg4ometry.mcnp.Cell(reg=reg)

    # SURFACES
    # box surfaces
    dim = [6.0, 10, 5.5]  # dimensions of a full block
    sPX1 = pyg4ometry.mcnp.PX((-dim[0] / 2), surfaceNumber=1, reg=reg)  # px1 (left)
    sPX2 = pyg4ometry.mcnp.PX((dim[0] / 2), surfaceNumber=2, reg=reg)  # px2 (right)
    sPY1 = pyg4ometry.mcnp.PY((-dim[1] / 2), surfaceNumber=3, reg=reg)  # py1 (bottom)
    sPY2 = pyg4ometry.mcnp.PY((dim[1] / 2), surfaceNumber=4, reg=reg)  # py2 (top)
    sPZ1 = pyg4ometry.mcnp.PZ((-dim[2] / 2), surfaceNumber=5, reg=reg)  # pz1 (back)
    sPZ2 = pyg4ometry.mcnp.PZ((dim[2] / 2), surfaceNumber=6, reg=reg)  # pz2 (front)
    # hole surface
    small = 0.02
    holePos = [0, 0, (dim[2]+small)/2]
    holeDir = [0, 0, -2-(small/2)]
    sRCC1 = pyg4ometry.mcnp.RCC(*holePos, *holeDir, 1.5, surfaceNumber=7, reg=reg)
    cBox.addSurfaces([sPX1, sPX2, sPY1, sPY2, sPZ1, sPZ2, sRCC1])
    # world surfaces
    sSO1 = pyg4ometry.mcnp.SO(15, reg=reg)
    cWorld.addSurface(cBox.geometry)  # cell2 is between box /
    cWorld.addSurface(sSO1)  # / and s1
    cVoid.addSurface(sSO1)  # cell3 is outside s1

    # MATERIAL
    m0 = pyg4ometry.mcnp.Material(0, reg=reg)
    m1 = pyg4ometry.mcnp.Material(1, -0.92, reg=reg)
    m2 = pyg4ometry.mcnp.Material(2, -0.001225, reg=reg)
    cVoid.addMaterial(m0)
    cBox.addMaterial(m1)
    cWorld.addMaterial(m2)

    # GEOMETRY
    geoOut = sSO1
    geoIn = pyg4ometry.mcnp.Complement(sSO1)
    # box geometry
    geomBoxX = pyg4ometry.mcnp.Intersection(sPX1, pyg4ometry.mcnp.Complement(sPX2))
    geomBoxY = pyg4ometry.mcnp.Intersection(sPY1, pyg4ometry.mcnp.Complement(sPY2))
    geomBoxZ = pyg4ometry.mcnp.Intersection(sPZ1, pyg4ometry.mcnp.Complement(sPZ2))
    geom1 = pyg4ometry.mcnp.Intersection(geomBoxZ, pyg4ometry.mcnp.Intersection(geomBoxY, geomBoxX))
    geomBox = pyg4ometry.mcnp.Intersection(geom1, sRCC1)
    cBox.addGeometry(geomBox)
    # world geometry
    geoWorld = pyg4ometry.mcnp.Intersection(sSO1, pyg4ometry.mcnp.Complement(cBox))
    cWorld.addGeometry(geoWorld)
    cVoid.addGeometry(geoOut)

    # IMPORTANCE
    i0 = pyg4ometry.mcnp.IMP("p", 0)
    i1 = pyg4ometry.mcnp.IMP("p", 1)
    cWorld.addImportance(i1)
    cBox.addImportance(i1)
    cVoid.addImportance(i0)

    if write:
        f = pyg4ometry.mcnp.Writer(columnMax=75)
        f.setTitle("SMALL BOX WITH HOLE")
        f.addGeometry(reg=reg)
        f.write("i-smallBoxWithHole.txt")

    if vis:
        v = pyg4ometry.visualisation.VtkViewer()
        v.addAxes()
        v.addMeshSimple(cBox.mesh())
        v.view()


# remove comment when debugging
test_smallBoxWithHole(True, True)

if __name__ == "__main__":
    test_smallBoxWithHole(True, True)
