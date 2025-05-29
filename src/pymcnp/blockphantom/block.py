import pyg4ometry
import numpy as _np


class Block:
    def __init__(self, blockType, translation=[0, 0, 0], rotation=[0, 0, 0], cellNumber=None):
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
        # HoleInfo is the [xyz co-ordinates, xyz axis-vector] for each hole FROM THE CENTER OF EACH BLOCK
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

        self.surfaces = self._makeSurfaces()
        self.geometry = self._makeGeometry()
        self.transform(rotation=rotation, translation=translation)
        self.cellNumber = cellNumber  # cell number

    def getHole(self, holeNumber):
        if holeNumber < 0 or holeNumber > 21:
            msg = "There are 22 holes per block. Please chose a hole number between 0 and 21."
            raise TypeError(msg)
        return self.holeInfo[holeNumber]  # [0] position, [1] vector (orientation and height)

    def _makeSurfaces(self):
        surfaces = []
        # polythene box
        surfaces.append(pyg4ometry.mcnp.PX((-self.dim[0] / 2)))  # px1 (left)
        surfaces.append(pyg4ometry.mcnp.PX((self.dim[0] / 2)))   # px2 (right)
        surfaces.append(pyg4ometry.mcnp.PY((-self.dim[1] / 2)))  # py1 (bottom)
        surfaces.append(pyg4ometry.mcnp.PY((self.dim[1] / 2)))   # py2 (top)
        surfaces.append(pyg4ometry.mcnp.PZ((-self.dim[2] / 2)))  # pz1 (back)
        surfaces.append(pyg4ometry.mcnp.PZ((self.dim[2] / 2)))   # pz2 (front)

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

    def _makeGeometry(self):
        # polythene box
        geomBoxX = pyg4ometry.mcnp.Intersection(self.surfaces[0], pyg4ometry.mcnp.Complement(self.surfaces[1]))
        geomBoxY = pyg4ometry.mcnp.Intersection(self.surfaces[2], pyg4ometry.mcnp.Complement(self.surfaces[3]))
        geomBoxZ = pyg4ometry.mcnp.Intersection(self.surfaces[4], pyg4ometry.mcnp.Complement(self.surfaces[5]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBoxZ, pyg4ometry.mcnp.Intersection(geomBoxY, geomBoxX))

        for i in range(6, 20+1):
            # connector holes tube, left, right, front
            geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(self.surfaces[i]))

        if self.blockType == "full":
            for i in range(21, 27+1):
                # connector holes back
                geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(self.surfaces[i]))

        return geomBlock

    def transform(self, rotation=[0, 0, 0], translation=[0, 0, 0]):
        rot = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

        if rotation == [0, 0, 0] and translation == [0, 0, 0]:
            # no rotation and no translation
            return self
        else:
            # rotation
            print(f"rotation")
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

        rotMat = _np.array(rot)
        transMat = _np.array(trans)

        # new block
        block_p = Block(blockType=self.blockType, cellNumber=self.cellNumber)  # copy block

        surfaces_new = []
        holeInfo_new = []

        for s in block_p.surfaces:
            s_new = s.transform(rotation=rot, translation=trans)
            surfaces_new.append(s_new)

        for hi in _np.array(self.holeInfo):
            pos_new = rotMat @ hi[0] + transMat
            dir_new = rotMat @ hi[1]
            holeInfo_new.append([pos_new, dir_new])

        block_p.surfaces = surfaces_new
        block_p.geometry = block_p._makeGeometry()
        block_p.holeInfo = holeInfo_new

        return block_p

    def holeToHoleTransform(self):
        pass

    def addToRegistry(self, registry, replace=False):
        registry.addCell(self, replace=replace)

    def makeNewConnectBlock(self, localHole, foreignHole, angle, makeConnector):
        """
        localHole: where on THIS block to connect this block to EG. blockThis.getHolePosition(hole=3)
        foreignHole: where on ANOTHER block to connect this block to EG. blockOther.getHolePosition(hole=3)
        angle: rotate around the axis of connection
        makeConnector : True
            returns [block, connector]
        makeConnector : False
            returns [block]
        """

        # work out transform from localHole to foreignHole

        # build block

        # transform block

        # build connector

        # transform connector

        #return [block_p, connector]


class Connector:
    def __init__(self, hole, cellNumber=None):
        self.hole = hole

        connector = pyg4ometry.mcnp.RCC(*block_p.position, *self.hole[1], 0.3)
        block_p.geometry = pyg4ometry.mcnp.Intersection(block_p.geometry, connector)