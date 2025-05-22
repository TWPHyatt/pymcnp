import pyg4ometry


class Block:
    def __init__(self, hole, x, y, z):
        small = 0.02  # to extend past the planes of the box by 0.01 cm so no inf small surface mesh covering hole
        unit = 5.5/2
        # x,y,z co-ordinates of the holes, and which axis they are aligned with
        self.holeInfo = [[[-unit, -(y+small)/2, 0], [0, 1, 0]],  # tube1, and y-aligned
                         [[unit, -(y+small)/2, 0], [0, 1, 0]]  # tube2, and y-aligned
                         [[], [0, 0, 1]],  # holeF1, and z-aligned
                         [[], [0, 0, 1]],  # holeF2, and z-aligned
                         [[], [0, 0, 1]],  # holeF3, and z-aligned
                         [[], [0, 0, 1]],  # holeF4, and z-aligned
                         [[], [0, 0, 1]],  # holeF5, and z-aligned
                         [[], [0, 0, 1]],  # holeF6, and z-aligned
                         [[], [0, 0, 1]],  # holeF7, and z-aligned
                         [[], [0, 0, 1]],  # holeB1, and z-aligned
                         [[], [0, 0, 1]],  # holeB2, and z-aligned
                         [[], [0, 0, 1]],  # holeB3, and z-aligned
                         [[], [0, 0, 1]],  # holeB4, and z-aligned
                         [[], [0, 0, 1]],  # holeB5, and z-aligned
                         [[], [0, 0, 1]],  # holeB6, and z-aligned
                         [[], [0, 0, 1]],  # holeB7, and z-aligned
                         [[], [1, 0, 0]],  # holeL1, and x-aligned
                         [[], [1, 0, 0]],  # holeL2, and x-aligned
                         [[], [1, 0, 0]],  # holeL3, and x-aligned
                         [[], [1, 0, 0]],  # holeR1, and x-aligned
                         [[], [1, 0, 0]],  # holeR2, and x-aligned
                         [[], [1, 0, 0]],  # holeR3, and x-aligned
                         ]
        self.hole = hole
        self.surface = None


    def makeBlock(self, hole, block=None, transform=[]):
        self.surface = self._makeBlockGeometry(11.0, 16.5, 5.5)

        # transform...
        return self.surface

    def makeHalfBlock(self, hole, block=None, transform=[]):
        self.surface = self._makeBlockGeometry(11.0, 16.5, 2.5)

        # transform...
        return self.surface

    def _makeBlockGeometry(self, x, y, z):


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

    def _makeBlockConnection(self, block1, block2, hole1, hole2):
        connector = pyg4ometry.mcnp.RCC(0, 0, 0, 0, 0, 2, 0.3)  #

"""
class Phantom
    init(registry)

    def makeBlock()
        block1 = Block()
        block2 = Block()

    def makeConnection()

Class Block
    init
        self.geometry = makeBlockGeometry()
        self.connections = {'top'    : { 1 : None, 2 : None}
                            'bottom' : { 1 : None, 2 : None}
                            'front'  : { 1 : None, 2 : None, 3 : None, 4 : None  ... }
                            'back'   : { 1 : None, 2 : None, 3 : None, 4 : None  ...}
                            'left    : { 1 : None, 2 : None, 3 : None }
                            'right   : { 1 : None, 2 : None, 3 : None }
                            }  # points to the other block objects?
        

    

"""