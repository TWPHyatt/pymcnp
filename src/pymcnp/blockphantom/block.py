import pyg4ometry
import numpy as _np

class Block:
    def __init__(self, blockType, blockNumber=None, position=None):
        self.blockType = blockType
        self.dim = [0, 0, 0]  # dimensions [width in x, height in y, depth in z]
        if position is None:
            self.position = [0, 0, 0]  # first block is centered on the origin
        else:
            self.position = position
        if blockType == "full":
            self.dim = [11.0, 16.5, 5.5]  # cm
        elif blockType == "half":
            self.dim = [11.0, 16.5, 2.5]  # cm
        else:
            msg = "Block type can only be `full` of `half`"
            raise TypeError(msg)

        self.holePattern = [[False*3], [False*3], [False*3], [False*3], [False*3], [False*3]]
        self.small = 0.02  # to extend past the planes of the box by 0.01 cm so no inf small surface mesh covering hole
        self.unit = self.dim[1]/(3*2)  # hole separation unit on the surface of the block
        # HoleInfo is the [xyz co-ordinates, xyz axis-vector] for each hole of a block
        self.holeInfo = [[[-self.unit, -(self.dim[1]+self.small)/2, 0], [0, self.dim[1]+self.small, 0]],       # tube1
                         [[+self.unit, -(self.dim[1]+self.small)/2, 0], [0, self.dim[1]+self.small, 0]],       # tube2
                         [[-self.unit, +self.unit*2, (self.dim[2]+self.small)/2], [0, 0, -1-(self.small/2)]],  # holeF1
                         [[+self.unit, +self.unit*2, (self.dim[2]+self.small)/2], [0, 0, -1-(self.small/2)]],  # holeF2
                         [[-self.unit, 0, (self.dim[2]+self.small)/2], [0, 0, -1-(self.small/2)]],             # holeF3
                         [[0, 0, (self.dim[2]+self.small)/2], [0, 0, -1-(self.small/2)]],                      # holeF4
                         [[+self.unit, 0, (self.dim[2]+self.small)/2], [0, 0, -1-(self.small/2)]],             # holeF5
                         [[-self.unit, -self.unit*2, (self.dim[2]+self.small)/2], [0, 0, -1-(self.small/2)]],  # holeF6
                         [[+self.unit, -self.unit*2, (self.dim[2]+self.small)/2], [0, 0, -1-(self.small/2)]],  # holeF7
                         [[-self.unit, +self.unit*2, -(self.dim[2]+self.small)/2], [0, 0, 1+(self.small/2)]],  # holeB1
                         [[+self.unit, +self.unit*2, -(self.dim[2]+self.small)/2], [0, 0, 1+(self.small/2)]],  # holeB2
                         [[-self.unit, 0, -(self.dim[2]+self.small)/2], [0, 0, 1+(self.small/2)]],             # holeB3
                         [[0, 0, -(self.dim[2]+self.small)/2], [0, 0, 1+(self.small/2)]],                      # holeB4
                         [[+self.unit, 0, -(self.dim[2]+self.small)/2], [0, 0, 1+(self.small/2)]],             # holeB5
                         [[-self.unit, -self.unit*2, -(self.dim[2]+self.small)/2], [0, 0, 1+(self.small/2)]],  # holeB6
                         [[+self.unit, -self.unit*2, -(self.dim[2]+self.small)/2], [0, 0, 1+(self.small/2)]],  # holeB7
                         [[-self.dim[0]/2, +self.unit*2, 0], [1+(self.small/2), 0, 0]],                        # holeL1
                         [[-self.dim[0]/2, 0, 0], [1+(self.small/2), 0, 0]],                                   # holeL2
                         [[-self.dim[0]/2, -self.unit*2, 0], [1+(self.small/2), 0, 0]],                        # holeL3
                         [[self.dim[0]/2, +self.unit*2, 0], [-1-(self.small/2), 0, 0]],                        # holeR1
                         [[self.dim[0]/2, 0, 0], [-1-(self.small/2), 0, 0]],                                   # holeR2
                         [[self.dim[0]/2, -self.unit*2, 0], [-1-(self.small/2), 0, 0]],                        # holeR3
                         ]

        self.surfaces = self._makeSurfaces(*self.dim)
        self.geometry = self._makeGeometry(self.surfaces)
        self.blockNumber = blockNumber  # cell number

    def getHolePos(self, holeNum):
        if holeNum < 0 or holeNum > 21:
            msg = "There are 22 holes per block. Please chose a hole number between 0 and 21."
            raise TypeError(msg)
        return self.holeInfo[holeNum][0]

    def _makeSurfaces(self, x, y, z):
        surfaces = []
        # polythene box
        surfaces.append(pyg4ometry.mcnp.PX(-x/2))  # px1
        surfaces.append(pyg4ometry.mcnp.PX(x/2))  # px2
        surfaces.append(pyg4ometry.mcnp.PY(-y/2))  # py1
        surfaces.append(pyg4ometry.mcnp.PY(y/2))  # py2
        surfaces.append(pyg4ometry.mcnp.PZ(-z/2))  # pz1
        surfaces.append(pyg4ometry.mcnp.PZ(z/2))  # pz2

        # source tubes bottom to top of block
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[0][0], *self.holeInfo[0][1], 0.4))  # left tube
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[1][0], *self.holeInfo[1][1], 0.4))  # right tube

        # connector holes front
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[2][0], *self.holeInfo[2][1], 0.4))  # front top left
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[3][0], *self.holeInfo[3][1], 0.4))  # front top right
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[4][0], *self.holeInfo[4][1], 0.4))  # front middle left
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[5][0], *self.holeInfo[5][1], 0.4))  # front middle center
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[6][0], *self.holeInfo[6][1], 0.4))  # front middle right
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[7][0], *self.holeInfo[7][1], 0.4))  # front bottom left
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[8][0], *self.holeInfo[8][1], 0.4))  # front bottom right

        # connector holes back
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[9][0], *self.holeInfo[9][1], 0.4))    # back top left
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[10][0], *self.holeInfo[10][1], 0.4))  # back top right
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[11][0], *self.holeInfo[11][1], 0.4))  # back middle left
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[12][0], *self.holeInfo[12][1], 0.4))  # back middle center
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[13][0], *self.holeInfo[13][1], 0.4))  # back middle right
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[14][0], *self.holeInfo[14][1], 0.4))  # back bottom left
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[15][0], *self.holeInfo[15][1], 0.4))  # back bottom right

        # connector holes left
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[16][0], *self.holeInfo[16][1], 0.4))  # left-side top
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[17][0], *self.holeInfo[17][1], 0.4))  # left-side middle
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[18][0], *self.holeInfo[18][1], 0.4))  # left-side bottom

        # connector holes right
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[19][0], *self.holeInfo[19][1], 0.4))  # right-side top
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[20][0], *self.holeInfo[20][1], 0.4))  # right-side middle
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[21][0], *self.holeInfo[21][1], 0.4))  # right-side bottom

        return surfaces

    def _makeGeometry(self, surfaces):
        # polythene box
        geomBoxX = pyg4ometry.mcnp.Intersection(surfaces[0], pyg4ometry.mcnp.Complement(surfaces[1]))
        geomBoxY = pyg4ometry.mcnp.Intersection(surfaces[2], pyg4ometry.mcnp.Complement(surfaces[3]))
        geomBoxZ = pyg4ometry.mcnp.Intersection(surfaces[4], pyg4ometry.mcnp.Complement(surfaces[5]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBoxZ, pyg4ometry.mcnp.Intersection(geomBoxY, geomBoxX))

        # source tubes bottom to top of block
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[6]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[7]))

        # connector holes front
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[8]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[9]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[10]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[11]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[12]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[13]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[14]))

        # connector holes back
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[15]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[16]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[17]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[18]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[19]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[20]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[21]))

        # connector holes left
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[22]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[23]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[24]))

        # connector holes right
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[25]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[26]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[27]))

        return geomBlock

    def transform(self, rotation=[0, 0, 0], translation=[0, 0, 0]):
        rot = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

        if rotation == [0, 0, 0] and translation == [0, 0, 0]:
            # no rotation and no translation
            return self
        else:
            # rotation

            # check for correct integer input
            if len(rotation) != 3 or not isinstance(rotation, list):
                msg = "rotation must be a list of length 3 for [x,y,z] rotations"
                raise TypeError(msg)
            if not all(isinstance(i, int) for i in rotation):
                msg = "rotation elements can only be integers : [x,y,z]" \
                      "\n rot[0] = 1 is a 90 degree positive rotation around the x axis " \
                      "\n rot[1] = 1 is a 90 degree positive rotation around the y axis" \
                      "\n rot[2] = 1 is a 90 degree positive rotation around the z axis"
                raise TypeError(msg)

            for h in range(0, len(self.holeInfo)):
                if rotation[0] > 0:
                    # rotate in x
                    theta = rotation[0] * _np.pi / 2
                    rot = [
                        [1.0, 0.0, 0.0],
                        [0.0, _np.cos(theta), -_np.sin(theta)],
                        [0.0, _np.sin(theta), _np.cos(theta)]
                    ]
                if rotation[1] > 0:
                    # rotate in y
                    theta = rotation[1] * _np.pi / 2
                    rot = [
                        [_np.cos(theta), 0.0, _np.sin(theta)],
                        [0.0, 1.0, 0.0],
                        [-_np.sin(theta), 0.0, _np.cos(theta)]
                    ]
                if rotation[2] > 0:
                    # rotate in z
                    theta = rotation[2] * _np.pi / 2
                    rot = [
                        [_np.cos(theta), -_np.sin(theta), 0.0],
                        [_np.sin(theta), _np.cos(theta), 0.0],
                        [0.0, 0.0, 1.0]
                    ]

            # translation
            if len(translation) != 3 or not isinstance(translation, list):
                msg = "translation must be a list of length 3 for [x,y,z] translation"
                raise TypeError(msg)
            trans = translation

        for surface in self.surfaces:  # do the transformation
            surface.transform(rotation=rot, translation=trans)

        return Block(self.blockType, self.position)

    def addToRegistry(self, registry, replace=False):
        for surface in self.surfaces:
            registry.addSurface(surface, replace=replace)
        return Block(self.blockType, self.position, )

    def _makeBlockConnection(self, block1, block2, hole1, hole2):
        connector = pyg4ometry.mcnp.RCC(0, 0, 0, 0, 0, 2, 0.3)  #
