import pyg4ometry
import numpy as _np


class Block:
    def __init__(self, blockType, position=None, cellNumber=None):
        self.blockType = blockType
        self.dim = [0, 0, 0]  # dimensions [width in x, height in y, depth in z]
        if blockType == "full":
            self.dim = [11.0, 16.5, 5.5]  # cm
        elif blockType == "half":
            self.dim = [11.0, 16.5, 2.5]  # cm
        else:
            msg = "Block type can only be `full` of `half`"
            raise TypeError(msg)

        self.small = 0.02  # to extend past the planes of the box by 0.01 cm so no inf small surface mesh covering hole
        self.unit = self.dim[1] / (3 * 2)  # hole separation unit on the surface of the block
        # HoleInfo is the [xyz co-ordinates, xyz axis-vector] for each hole of a block
        self.holeInfo = [[[-self.unit, -(self.dim[1] + self.small) / 2, 0], [0, self.dim[1] + self.small, 0]],  # tube1
                         [[+self.unit, -(self.dim[1] + self.small) / 2, 0], [0, self.dim[1] + self.small, 0]],  # tube2
                         [[-self.dim[0] / 2, +self.unit * 2, 0], [1 + (self.small / 2), 0, 0]],  # holeL1
                         [[-self.dim[0] / 2, 0, 0], [1 + (self.small / 2), 0, 0]],  # holeL2
                         [[-self.dim[0] / 2, -self.unit * 2, 0], [1 + (self.small / 2), 0, 0]],  # holeL3
                         [[self.dim[0] / 2, +self.unit * 2, 0], [-1 - (self.small / 2), 0, 0]],  # holeR1
                         [[self.dim[0] / 2, 0, 0], [-1 - (self.small / 2), 0, 0]],  # holeR2
                         [[self.dim[0] / 2, -self.unit * 2, 0], [-1 - (self.small / 2), 0, 0]],  # holeR3
                         [[-self.unit, +self.unit * 2, (self.dim[2] + self.small) / 2], [0, 0, -1 - (self.small / 2)]],
                         # holeF1
                         [[+self.unit, +self.unit * 2, (self.dim[2] + self.small) / 2], [0, 0, -1 - (self.small / 2)]],
                         # holeF2
                         [[-self.unit, 0, (self.dim[2] + self.small) / 2], [0, 0, -1 - (self.small / 2)]],  # holeF3
                         [[0, 0, (self.dim[2] + self.small) / 2], [0, 0, -1 - (self.small / 2)]],  # holeF4
                         [[+self.unit, 0, (self.dim[2] + self.small) / 2], [0, 0, -1 - (self.small / 2)]],  # holeF5
                         [[-self.unit, -self.unit * 2, (self.dim[2] + self.small) / 2], [0, 0, -1 - (self.small / 2)]],
                         # holeF6
                         [[+self.unit, -self.unit * 2, (self.dim[2] + self.small) / 2], [0, 0, -1 - (self.small / 2)]],
                         # holeF7
                         ]
        if self.blockType == "full":
            self.holeInfo.append(
                [[-self.unit, +self.unit * 2, -(self.dim[2] + self.small) / 2], [0, 0, 1 + (self.small / 2)]])  # holeB1
            self.holeInfo.append(
                [[+self.unit, +self.unit * 2, -(self.dim[2] + self.small) / 2], [0, 0, 1 + (self.small / 2)]])  # holeB2
            self.holeInfo.append(
                [[-self.unit, 0, -(self.dim[2] + self.small) / 2], [0, 0, 1 + (self.small / 2)]])  # holeB3
            self.holeInfo.append([[0, 0, -(self.dim[2] + self.small) / 2], [0, 0, 1 + (self.small / 2)]])  # holeB4
            self.holeInfo.append(
                [[+self.unit, 0, -(self.dim[2] + self.small) / 2], [0, 0, 1 + (self.small / 2)]])  # holeB5
            self.holeInfo.append(
                [[-self.unit, -self.unit * 2, -(self.dim[2] + self.small) / 2], [0, 0, 1 + (self.small / 2)]])  # holeB6
            self.holeInfo.append(
                [[+self.unit, -self.unit * 2, -(self.dim[2] + self.small) / 2], [0, 0, 1 + (self.small / 2)]])  # holeB7

        if position is None:
            self.position = [0, 0, 0]  # first block is centered on the origin
        else:
            self.position = position
            for i, hInfo in enumerate(self.holeInfo):
                for j, hCoords in enumerate(hInfo[0]):
                    self.holeInfo[i][0][j] = hCoords + self.position[j]

        self.surfaces = self._makeSurfaces()
        self.geometry = self._makeGeometry(self.surfaces)
        self.cellNumber = cellNumber  # cell number

    def getHole(self, holeNumber):
        if holeNumber < 0 or holeNumber > 21:
            msg = "There are 22 holes per block. Please chose a hole number between 0 and 21."
            raise TypeError(msg)
        return self.holeInfo[holeNumber]  # [0] position, [1] vector (orientation and height)

    def _makeSurfaces(self):
        surfaces = []
        # polythene box
        surfaces.append(pyg4ometry.mcnp.PX((-self.dim[0] / 2) + self.position[0]))  # px1 (left)
        surfaces.append(pyg4ometry.mcnp.PX((self.dim[0] / 2) + self.position[0]))  # px2 (right)
        surfaces.append(pyg4ometry.mcnp.PY((-self.dim[1] / 2) + self.position[1]))  # py1 (bottom)
        surfaces.append(pyg4ometry.mcnp.PY((self.dim[1] / 2) + self.position[1]))  # py2 (top)
        surfaces.append(pyg4ometry.mcnp.PZ((-self.dim[2] / 2) + self.position[2]))  # pz1 (back)
        surfaces.append(pyg4ometry.mcnp.PZ((self.dim[2] / 2) + self.position[2]))  # pz2 (front)

        # source tubes bottom to top of block
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[0][0], *self.holeInfo[0][1], 0.4))  # left tube
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[1][0], *self.holeInfo[1][1], 0.4))  # right tube

        # connector holes left
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[2][0], *self.holeInfo[2][1], 0.4))  # left-side top
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[3][0], *self.holeInfo[3][1], 0.4))  # left-side middle
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[4][0], *self.holeInfo[4][1], 0.4))  # left-side bottom

        # connector holes right
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[5][0], *self.holeInfo[5][1], 0.4))  # right-side top
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[6][0], *self.holeInfo[6][1], 0.4))  # right-side middle
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[7][0], *self.holeInfo[7][1], 0.4))  # right-side bottom

        # connector holes front
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[8][0], *self.holeInfo[8][1], 0.4))  # front top left
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[9][0], *self.holeInfo[9][1], 0.4))  # front top right
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[10][0], *self.holeInfo[10][1], 0.4))  # front middle left
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[11][0], *self.holeInfo[11][1], 0.4))  # front middle cent
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[12][0], *self.holeInfo[12][1], 0.4))  # front middle right
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[13][0], *self.holeInfo[13][1], 0.4))  # front bottom left
        surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[14][0], *self.holeInfo[14][1], 0.4))  # front bottom right

        if self.blockType == "full":
            # connector holes back
            surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[15][0], *self.holeInfo[15][1], 0.4))  # back top left
            surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[16][0], *self.holeInfo[16][1], 0.4))  # back top right
            surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[17][0], *self.holeInfo[17][1], 0.4))  # back middle left
            surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[18][0], *self.holeInfo[18][1], 0.4))  # back middle cent
            surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[19][0], *self.holeInfo[19][1], 0.4))  # back middle right
            surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[20][0], *self.holeInfo[20][1], 0.4))  # back bottom left
            surfaces.append(pyg4ometry.mcnp.RCC(*self.holeInfo[21][0], *self.holeInfo[21][1], 0.4))  # back bottom right

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

        # connector holes left
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[8]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[9]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[10]))

        # connector holes right
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[11]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[12]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[13]))

        # connector holes front
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[14]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[15]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[16]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[17]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[18]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[19]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[20]))

        if self.blockType == "full":
            # connector holes back
            geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[21]))
            geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[22]))
            geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[23]))
            geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[24]))
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

        # new block
        block_p = Block(blockType=self.blockType, position=self.position, cellNumber=self.cellNumber)

        # override position
        block_p.position = [self.position[0] + trans[0], self.position[1] + trans[1], self.position[2] + trans[2]]

        # override holeInfo
        for i, hInfo in enumerate(self.holeInfo):
            for j, hCoords in enumerate(hInfo[0]):
                block_p.holeInfo[i][0][j] = hCoords + block_p.position[j]

        # override surfaces
        surfaces_p = []
        for surface in self.surfaces:  # do the transformation
            surfaces_p.append(surface.transform(rotation=rot, translation=trans))
        block_p.surfaces = surfaces_p

        # override geometry
        block_p.geometry = block_p._makeGeometry(block_p.surfaces)

        return block_p

    def addToRegistry(self, registry, replace=False):
        registry.addCell(self, replace=replace)

    def connectBlock(self, hole1, hole2):
        """
        hole1: where on ANOTHER block to connect this block to EG. blockOther.getHolePosition(hole=3)
        hole2: where on THIS block to connect this block to EG. blockThis.getHolePosition(hole=3)
        """
        print("h1: ", hole1)
        print("h2: ", hole2)
        # check if two holes are oriented opposite directions
        h1 = _np.array(hole1[1])
        h2 = _np.array(hole2[1])
        h1Norm = h1 / _np.linalg.norm(hole1[1])
        h2Norm = h2 / _np.linalg.norm(hole2[1])

        if _np.abs(_np.dot(h1Norm, h2Norm)) < 1e-6:
            # hole1 and hole2 are opposite
            rot = [0, 0, 0]
        else:
            # hole1 and hole2 are not opposite
            # the target is normalised hole1 vector scaled to match hole2 vector magnitude
            target = -h1Norm * _np.linalg.norm(h2)
            h2_p = [0, 0, 0]
            rot = [0, 0, 0]

            for ix in range(4):
                for iy in range(4):
                    for iz in range(4):
                        for a in range(ix):
                            h2_p = [h2[0], -h2[2], h2[1]]
                        for b in range(iy):
                            h2_p = [h2_p[0], -h2_p[2], h2_p[1]]
                        for c in range(iz):
                            h2_p = [h2_p[0], -h2_p[2], h2_p[1]]
                        if _np.allclose(h2_p, target):
                            rot = [ix, iy, iz]
                        else:
                            rot = [0, 0, 0]

        block_p = self.transform(rotation=rot, translation=hole1[0]).transform(rotation=[0, 0, 0], translation=[-el for el in hole2[0]])
        return block_p


class Connector:
    def __init__(self, hole, cellNumber=None):
        self.hole = hole

        connector = pyg4ometry.mcnp.RCC(*block_p.position, *self.hole[1], 0.3)
        block_p.geometry = pyg4ometry.mcnp.Intersection(block_p.geometry, connector)