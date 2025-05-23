import pyg4ometry
import numpy as _np

class Block:
    def __init__(self, blockType, pos1, pos2=None, rot=None):
        pos1 = None  # first block has pos1 = None
        pos2 = None
        if blockType == "full":
            x, y, z = 11.0, 16.5, 5.5  # cm
        elif blockType == "half":
            x, y, z = 11.0, 16.5, 2.5  # cm
        else:
            msg = "Block type can only be `full` of `half`"
            raise TypeError(msg)

        identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        theta = 0

        if rot is None:
            pass  # no rotation
        else:
            if not all(isinstance(i, int) for i in rot):
                msg = "rot elements can only be integers : [x,y,z]" \
                      "\n rot[0] = 1 is a 90 degree rotation around the x axis, " \
                      "\n rot[1] = 1 is a 90 degree rotation around the y axis," \
                      "\n rot[2] = 1 is a 90 degree rotation around the z axis,"
                raise TypeError(msg)

        self.small = 0.02  # to extend past the planes of the box by 0.01 cm so no inf small surface mesh covering hole
        self.unit = y/(3*2)  # hole separation unit on the surface of the block
        # HoleInfo is the [xyz co-ordinates, xyz axis-vector] for each hole of a block
        self.holeInfo = [[[-self.unit, -(y+self.small)/2, 0], [0, y+self.small, 0]],                 # tube1
                         [[+self.unit, -(y+self.small)/2, 0], [0, y+self.small, 0]],                 # tube2
                         [[-self.unit, +self.unit*2, (z+self.small)/2], [0, 0, -1-(self.small/2)]],  # holeF1
                         [[+self.unit, +self.unit*2, (z+self.small)/2], [0, 0, -1-(self.small/2)]],  # holeF2
                         [[-self.unit, 0, (z+self.small)/2], [0, 0, -1-(self.small/2)]],             # holeF3
                         [[0, 0, (z+self.small)/2], [0, 0, -1-(self.small/2)]],                      # holeF4
                         [[+self.unit, 0, (z+self.small)/2], [0, 0, -1-(self.small/2)]],             # holeF5
                         [[-self.unit, -self.unit*2, (z+self.small)/2], [0, 0, -1-(self.small/2)]],  # holeF6
                         [[+self.unit, -self.unit*2, (z+self.small)/2], [0, 0, -1-(self.small/2)]],  # holeF7
                         [[-self.unit, +self.unit*2, -(z+self.small)/2], [0, 0, 1+(self.small/2)]],  # holeB1
                         [[+self.unit, +self.unit*2, -(z+self.small)/2], [0, 0, 1+(self.small/2)]],  # holeB2
                         [[-self.unit, 0, -(z+self.small)/2], [0, 0, 1+(self.small/2)]],             # holeB3
                         [[0, 0, -(z+self.small)/2], [0, 0, 1+(self.small/2)]],                      # holeB4
                         [[+self.unit, 0, -(z+self.small)/2], [0, 0, 1+(self.small/2)]],             # holeB5
                         [[-self.unit, -self.unit*2, -(z+self.small)/2], [0, 0, 1+(self.small/2)]],  # holeB6
                         [[+self.unit, -self.unit*2, -(z+self.small)/2], [0, 0, 1+(self.small/2)]],  # holeB7
                         [[-x/2, +self.unit*2, 0], [1+(self.small/2), 0, 0]],                        # holeL1
                         [[-x/2, 0, 0], [1+(self.small/2), 0, 0]],                                   # holeL2
                         [[-x/2, -self.unit*2, 0], [1+(self.small/2), 0, 0]],                        # holeL3
                         [[x/2, +self.unit*2, 0], [-1-(self.small/2), 0, 0]],                        # holeR1
                         [[x/2, 0, 0], [-1-(self.small/2), 0, 0]],                                   # holeR2
                         [[x/2, -self.unit*2, 0], [-1-(self.small/2), 0, 0]],                        # holeR3
                         ]

        if rot is not None:
            for h in range(0, len(self.holeInfo)):
                if rot[0] > 0:
                    # rotate in x
                    theta = rot[0] * _np.pi / 2
                    rotation = [
                        [1.0, 0.0, 0.0],
                        [0.0, _np.cos(theta), -_np.sin(theta)],
                        [0.0, _np.sin(theta), _np.cos(theta)]
                    ]
                    self.holeInfo[h][0] = _np.array(rotation) @ _np.array(self.holeInfo[h][0])  # xyz coordinates
                    self.holeInfo[h][1] = _np.array(rotation) @ _np.array(self.holeInfo[h][1])  # xyz axis-vector

                if rot[1] > 0:
                    # rotate in y
                    theta = rot[1] * _np.pi / 2
                    rotation = [
                        [_np.cos(theta), 0.0, _np.sin(theta)],
                        [0.0, 1.0, 0.0],
                        [-_np.sin(theta), 0.0, _np.cos(theta)]
                    ]
                    self.holeInfo[h][0] = _np.array(rotation) @ _np.array(self.holeInfo[h][0])  # xyz coordinates
                    self.holeInfo[h][1] = _np.array(rotation) @ _np.array(self.holeInfo[h][1])  # xyz axis-vector

                if rot[2] > 0:
                    # rotate in z
                    theta = rot[2] * _np.pi / 2
                    rotation = [
                        [_np.cos(theta), -_np.sin(theta), 0.0],
                        [_np.sin(theta), _np.cos(theta), 0.0],
                        [0.0, 0.0, 1.0]
                    ]
                    self.holeInfo[h][0] = _np.array(rotation) @ _np.array(self.holeInfo[h][0])  # xyz coordinates
                    self.holeInfo[h][1] = _np.array(rotation) @ _np.array(self.holeInfo[h][1])  # xyz axis-vector

        self.surface = self._makeBlockGeometry(x, y, z)

    def getHolePos(self, holeNum):
        if holeNum < 0 or holeNum > 21:
            msg = "There are 22 holes per block. Please chose a hole number between 0 and 21."
            raise TypeError(msg)
        return self.holeInfo[holeNum][0]

    def _makeBlockGeometry(self, x, y, z):

        # polythene box
        px1 = pyg4ometry.mcnp.P(1, 0, 0, -x/2)
        px2 = pyg4ometry.mcnp.P(1, 0, 0, x/2)
        py1 = pyg4ometry.mcnp.P(0, 1, 0, -y/2)
        py2 = pyg4ometry.mcnp.P(0, 1, 0, y/2)
        pz1 = pyg4ometry.mcnp.P(0, 0, 1, -z/2)
        pz2 = pyg4ometry.mcnp.P(0, 0, 1, z/2)
        geomBoxY = pyg4ometry.mcnp.Intersection(py1, pyg4ometry.mcnp.Complement(py2))
        geomBoxX = pyg4ometry.mcnp.Intersection(px1, pyg4ometry.mcnp.Complement(px2))
        geomBoxZ = pyg4ometry.mcnp.Intersection(pz1, pyg4ometry.mcnp.Complement(pz2))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBoxZ, pyg4ometry.mcnp.Intersection(geomBoxY, geomBoxX))

        # source tubes bottom to top of block
        tube1 = pyg4ometry.mcnp.RCC(*self.holeInfo[0][0], *self.holeInfo[0][1], 0.4)  # left
        tube2 = pyg4ometry.mcnp.RCC(*self.holeInfo[1][0], *self.holeInfo[1][1], 0.4)  # right
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(tube1))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(tube2))

        # connector holes front
        holeF1 = pyg4ometry.mcnp.RCC(*self.holeInfo[2][0], *self.holeInfo[2][1], 0.4)  # top left
        holeF2 = pyg4ometry.mcnp.RCC(*self.holeInfo[3][0], *self.holeInfo[3][1], 0.4)  # top right
        holeF3 = pyg4ometry.mcnp.RCC(*self.holeInfo[4][0], *self.holeInfo[4][1], 0.4)  # middle left
        holeF4 = pyg4ometry.mcnp.RCC(*self.holeInfo[5][0], *self.holeInfo[5][1], 0.4)  # middle center
        holeF5 = pyg4ometry.mcnp.RCC(*self.holeInfo[6][0], *self.holeInfo[6][1], 0.4)  # middle right
        holeF6 = pyg4ometry.mcnp.RCC(*self.holeInfo[7][0], *self.holeInfo[7][1], 0.4)  # bottom left
        holeF7 = pyg4ometry.mcnp.RCC(*self.holeInfo[8][0], *self.holeInfo[8][1], 0.4)  # bottom right
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeF1))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeF2))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeF3))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeF4))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeF5))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeF6))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeF7))

        # connector holes back
        holeB1 = pyg4ometry.mcnp.RCC(*self.holeInfo[9][0], *self.holeInfo[9][1], 0.4)   # top left
        holeB2 = pyg4ometry.mcnp.RCC(*self.holeInfo[10][0], *self.holeInfo[10][1], 0.4)  # top right
        holeB3 = pyg4ometry.mcnp.RCC(*self.holeInfo[11][0], *self.holeInfo[11][1], 0.4)  # middle left
        holeB4 = pyg4ometry.mcnp.RCC(*self.holeInfo[12][0], *self.holeInfo[12][1], 0.4)  # middle center
        holeB5 = pyg4ometry.mcnp.RCC(*self.holeInfo[13][0], *self.holeInfo[13][1], 0.4)  # middle right
        holeB6 = pyg4ometry.mcnp.RCC(*self.holeInfo[14][0], *self.holeInfo[14][1], 0.4)  # bottom left
        holeB7 = pyg4ometry.mcnp.RCC(*self.holeInfo[15][0], *self.holeInfo[15][1], 0.4)  # bottom right
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeB1))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeB2))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeB3))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeB4))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeB5))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeB6))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeB7))

        # connector holes left
        holeL1 = pyg4ometry.mcnp.RCC(*self.holeInfo[16][0], *self.holeInfo[16][1], 0.4)  # top
        holeL2 = pyg4ometry.mcnp.RCC(*self.holeInfo[17][0], *self.holeInfo[17][1], 0.4)  # middle
        holeL3 = pyg4ometry.mcnp.RCC(*self.holeInfo[18][0], *self.holeInfo[18][1], 0.4)  # bottom
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeL1))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeL2))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeL3))

        # connector holes right
        holeR1 = pyg4ometry.mcnp.RCC(*self.holeInfo[19][0], *self.holeInfo[19][1], 0.4)  # top
        holeR2 = pyg4ometry.mcnp.RCC(*self.holeInfo[20][0], *self.holeInfo[20][1], 0.4)  # middle
        holeR3 = pyg4ometry.mcnp.RCC(*self.holeInfo[21][0], *self.holeInfo[21][1], 0.4)  # bottom
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeR1))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeR2))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(holeR3))

        return geomBlock

    def _makeBlockConnection(self, block1, block2, hole1, hole2):
        connector = pyg4ometry.mcnp.RCC(0, 0, 0, 0, 0, 2, 0.3)  #
