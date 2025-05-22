import pyg4ometry


class Block:
    def __init__(self, hole):
        self.hole = hole
        self.surface = None

    def makeBlock(self, hole, block=None, transform=[]):
        gometry = self._makeBlockGeometry(11.0, 16.5, 5.5)
        self.surface = gometry
        # transform...
        return self.surface

    def makeHalfBlock(self, hole, block=None, transform=[]):
        gometry = self._makeBlockGeometry(11.0, 16.5, 2.5)
        self.surface = gometry
        # transform...
        return self.surface

    def _makeBlockGeometry(self, x, y, z):
        small = 0.02  # to extend past the planes of the box by 0.01 cm so no inf small surface mesh covering hole
        unit = 5.5/2

        # polythene box
        px1 = pyg4ometry.mcnp.PX(-x/2)
        px2 = pyg4ometry.mcnp.PX(x/2)
        py1 = pyg4ometry.mcnp.PY(-y/2)
        py2 = pyg4ometry.mcnp.PY(y/2)
        pz1 = pyg4ometry.mcnp.PZ(-z/2)
        pz2 = pyg4ometry.mcnp.PZ(z/2)
        geomBoxY = pyg4ometry.mcnp.Intersection(py1, pyg4ometry.mcnp.Complement(py2))
        geomBoxX = pyg4ometry.mcnp.Intersection(px1, pyg4ometry.mcnp.Complement(px2))
        geomBoxZ = pyg4ometry.mcnp.Intersection(pz1, pyg4ometry.mcnp.Complement(pz2))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBoxZ, pyg4ometry.mcnp.Intersection(geomBoxY, geomBoxX))

        # source tubes bottom to top of block
        tube1 = pyg4ometry.mcnp.RCC(-unit, -(y+small)/2, 0, 0, y+small, 0, 0.4)  # left
        tube2 = pyg4ometry.mcnp.RCC(unit, -(y+small)/2, 0, 0, y+small, 0, 0.4)  # right
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(tube1))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(tube2))

        # connector holes front
        holeF1 = pyg4ometry.mcnp.RCC(-unit, +unit*2, (z+small)/2, 0, 0, -1-(small/2), 0.4)  # top left
        holeF2 = pyg4ometry.mcnp.RCC(+unit, +unit*2, (z+small)/2, 0, 0, -1-(small/2), 0.4)  # top right
        holeF3 = pyg4ometry.mcnp.RCC(-unit, 0, (z+small)/2, 0, 0, -1-(small/2), 0.4)  # middle left
        holeF4 = pyg4ometry.mcnp.RCC(0, 0, (z+small)/2, 0, 0, -1-(small/2), 0.4)  # middle center
        holeF5 = pyg4ometry.mcnp.RCC(+unit, 0, (z+small)/2, 0, 0, -1-(small/2), 0.4)  # middle right
        holeF6 = pyg4ometry.mcnp.RCC(-unit, -unit*2, (z+small)/2, 0, 0, -1-(small/2), 0.4)  # bottom left
        holeF7 = pyg4ometry.mcnp.RCC(+unit, -unit*2, (z+small)/2, 0, 0, -1-(small/2), 0.4)  # bottom right
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeF1))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeF2))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeF3))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeF4))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeF5))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeF6))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeF7))

        # connector holes back
        holeB1 = pyg4ometry.mcnp.RCC(-unit, +unit*2, -(z+small)/2, 0, 0, 1+(small/2), 0.4)  # top left
        holeB2 = pyg4ometry.mcnp.RCC(+unit, +unit*2, -(z+small)/2, 0, 0, 1+(small/2), 0.4)  # top right
        holeB3 = pyg4ometry.mcnp.RCC(-unit, 0, -(z+small)/2, 0, 0, 1+(small/2), 0.4)  # middle left
        holeB4 = pyg4ometry.mcnp.RCC(0, 0, -(z+small)/2, 0, 0, 1+(small/2), 0.4)  # middle center
        holeB5 = pyg4ometry.mcnp.RCC(+unit, 0, -(z+small)/2, 0, 0, 1+(small/2), 0.4)  # middle right
        holeB6 = pyg4ometry.mcnp.RCC(-unit, -unit*2, -(z+small)/2, 0, 0, 1+(small/2), 0.4)  # bottom left
        holeB7 = pyg4ometry.mcnp.RCC(+unit, -unit*2, -(z+small)/2, 0, 0, 1+(small/2), 0.4)  # bottom right
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeB1))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeB2))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeB3))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeB4))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeB5))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeB6))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeB7))

        # connector holes left
        holeL1 = pyg4ometry.mcnp.RCC(-x/2, +unit*2, 0, 1+(small/2), 0, 0, 0.4)  # top
        holeL2 = pyg4ometry.mcnp.RCC(-x/2, 0, 0, 1+(small/2), 0, 0, 0.4)  # middle
        holeL3 = pyg4ometry.mcnp.RCC(-x/2, -unit*2, 0, 1+(small/2), 0, 0, 0.4)  # bottom
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeL1))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeL2))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeL3))

        # connector holes right
        holeR1 = pyg4ometry.mcnp.RCC(x/2, +unit*2, 0, -1-(small/2), 0, 0, 0.4)  # top
        holeR2 = pyg4ometry.mcnp.RCC(x/2, 0, 0, -1-(small/2), 0, 0, 0.4)  # middle
        holeR3 = pyg4ometry.mcnp.RCC(x/2, -unit*2, 0, -1-(small/2), 0, 0, 0.4)  # bottom
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeR1))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeR2))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeR3))

        return geomBlock
