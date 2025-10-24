import pyg4ometry
import pymcnp
import numpy as np

def test_gridBlocksWrite(write=False):
    """
    test writing two blocks to MCNP input file
    :param write: write to file
    :type write: boolean
    :return: none
    """

    reg = pyg4ometry.mcnp.Registry()

    # CELLS
    # --- 4X4 GRID OF BLOCKS ---
    gridRowNum = 5
    gridColNum = 5
    removeNum = 7
    connectorNum = 6
    holeOrder = [0, 1, 2, 3, 7, 9, 4, 6]
    startX, startY, zDisp = -40, -40, 0
    deltaX, deltaY = 20, 20
    grid = [[[] for _ in range(gridColNum)] for _ in range(gridRowNum)]

    totalNumGridElements = (gridRowNum * gridColNum) - removeNum
    NumFullRows, remainder = divmod(totalNumGridElements, gridColNum)

    grid = []

    for r in range(NumFullRows):
        grid.append([0 for _ in range(gridColNum)])

    if remainder:
        grid.append([0 for _ in range(remainder)])

    for r, row in enumerate(grid):
        for c, el in enumerate(row):
            x = startX + c * deltaX
            y = startY + r * deltaY
            block = pymcnp.blockphantom.Block("full", rotationSteps=[0, 0, 0], reg=reg)
            block_p = block.transform(translation=[x, y, zDisp], rotation=[0, 0, 0], isRotationMatrix=False)
            grid[r][c] = [block_p]
            reg.addCell(block_p, replace=True)
            reg.addSurfaces(block_p.surfaceList, replace=True)
            for i in range(connectorNum):
                connector = block_p.addConnector(localHole=holeOrder[i], reg=reg)
                grid[r][c].append(connector)
                reg.addCell(connector, replace=True)
                reg.addSurfaces(connector.surfaceList, replace=True)

    # --- WORLD ---
    cWorld = pyg4ometry.mcnp.Cell(reg=reg)
    cVoid = pyg4ometry.mcnp.Cell(reg=reg)

    # SURFACES
    sSO1 = pyg4ometry.mcnp.SO(100, reg=reg)
    for r, row in enumerate(grid):
        for c, gridElement in enumerate(row):
            for cell in gridElement:
                cWorld.addSurface(cell.geometry)
    cWorld.addSurface(sSO1)
    cVoid.addSurface(sSO1)

    # GEOMETRY
    geoOut = sSO1
    geo = pyg4ometry.mcnp.Complement(sSO1)
    for r, row in enumerate(grid):
        for c, gridElement in enumerate(row):
            for cell in gridElement:
                geo = pyg4ometry.mcnp.Intersection(geo, pyg4ometry.mcnp.Complement(cell))
    cWorld.addGeometry(geo)
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
    for r, row in enumerate(grid):
        for c, gridElement in enumerate(row):
            for cell in gridElement:
                cell.addImportance(i1)
    cWorld.addImportance(i1)
    cVoid.addImportance(i0)

    if write:
        f = pyg4ometry.mcnp.Writer(columnMax=60)
        f.setTitle(f"{(gridRowNum*gridColNum)-removeNum} BLOCKS IN A {gridRowNum}X{gridColNum} GRID")
        f.addGeometry(reg=reg)
        f.write(f"i-{(gridRowNum*gridColNum)-removeNum}Blocks-{gridRowNum}X{gridColNum}Grid-Cons{connectorNum}.txt")

# remove comment when debugging
test_gridBlocksWrite(True)

if __name__ == "__main__":
    test_gridBlocksWrite(True)
