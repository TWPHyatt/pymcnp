import pyg4ometry
import numpy as _np

class Block(pyg4ometry.mcnp.Cell):
    def __init__(self, blockType, translation=[0, 0, 0], rotation=[0, 0, 0], cellNumber=None, reg=None):
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

        surfaces = self._makeSurfaces()

        # apply transformation
        rot = _np.eye(3)
        if rotation != [0, 0, 0]:
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

            if rotation[0] > 0:
                # rotate in x
                theta = rotation[0] * _np.pi / 2
                rotX = _np.array([
                    [1.0, 0.0, 0.0],
                    [0.0, _np.cos(theta), -_np.sin(theta)],
                    [0.0, _np.sin(theta), _np.cos(theta)]
                ])
                rot = rotX @ rot
            if rotation[1] > 0:
                # rotate in y
                theta = rotation[1] * _np.pi / 2
                rotY = _np.array([
                    [_np.cos(theta), 0.0, _np.sin(theta)],
                    [0.0, 1.0, 0.0],
                    [-_np.sin(theta), 0.0, _np.cos(theta)]
                ])
                rot = rotY @ rot
            if rotation[2] > 0:
                # rotate in z
                theta = rotation[2] * _np.pi / 2
                rotZ = _np.array([
                    [_np.cos(theta), -_np.sin(theta), 0.0],
                    [_np.sin(theta), _np.cos(theta), 0.0],
                    [0.0, 0.0, 1.0]
                ])
                rot = rotZ @ rot

        # translation
        if len(translation) != 3 or not isinstance(translation, list):
            msg = "translation must be a list of length 3 for [x,y,z] translation"
            raise TypeError(msg)

        rotMat = _np.array(rot)
        transMat = _np.array(translation)

        surfaces_new = []
        holeInfo_new = []

        for s in surfaces:
            s_new = s.transform(rotation=rot, translation=translation)
            surfaces_new.append(s_new)
        geometry = self._makeGeometry(surfaces_new)

        for hi in _np.array(self.holeInfo):
            pos_new = rotMat @ hi[0] + transMat
            dir_new = rotMat @ hi[1]
            holeInfo_new.append([pos_new, dir_new])
        self.holeInfo = holeInfo_new

        super().__init__(surfaces=surfaces_new, geometry=geometry, cellNumber=cellNumber, reg=reg)

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

    def _makeGeometry(self, surfaces):
        # polythene box
        geomBoxX = pyg4ometry.mcnp.Intersection(surfaces[0], pyg4ometry.mcnp.Complement(surfaces[1]))
        geomBoxY = pyg4ometry.mcnp.Intersection(surfaces[2], pyg4ometry.mcnp.Complement(surfaces[3]))
        geomBoxZ = pyg4ometry.mcnp.Intersection(surfaces[4], pyg4ometry.mcnp.Complement(surfaces[5]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBoxZ, pyg4ometry.mcnp.Intersection(geomBoxY, geomBoxX))

        for i in range(6, 20+1):
            # connector holes tube, left, right, front
            geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[i]))

        if self.blockType == "full":
            for i in range(21, 27+1):
                # connector holes back
                geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[i]))

        return geomBlock

    def transform(self, rotation=[0, 0, 0], translation=[0, 0, 0]):
        rot = _np.eye(3)

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

            if rotation[0] > 0:
                # rotate in x
                theta = rotation[0] * _np.pi / 2
                rotX = _np.array([
                    [1.0, 0.0, 0.0],
                    [0.0, _np.cos(theta), -_np.sin(theta)],
                    [0.0, _np.sin(theta), _np.cos(theta)]
                ])
                rot = rotX @ rot
            if rotation[1] > 0:
                # rotate in y
                theta = rotation[1] * _np.pi / 2
                rotY = _np.array([
                    [_np.cos(theta), 0.0, _np.sin(theta)],
                    [0.0, 1.0, 0.0],
                    [-_np.sin(theta), 0.0, _np.cos(theta)]
                ])
                rot = rotY @ rot
            if rotation[2] > 0:
                # rotate in z
                theta = rotation[2] * _np.pi / 2
                rotZ = _np.array([
                    [_np.cos(theta), -_np.sin(theta), 0.0],
                    [_np.sin(theta), _np.cos(theta), 0.0],
                    [0.0, 0.0, 1.0]
                ])
                rot = rotZ @ rot

            # translation
            if len(translation) != 3 or not isinstance(translation, list):
                msg = "translation must be a list of length 3 for [x,y,z] translation"
                raise TypeError(msg)

        rotMat = _np.array(rot)
        transMat = _np.array(translation)

        # new block
        block_p = Block(blockType=self.blockType, cellNumber=self.cellNumber)  # copy block

        surfaces_p = []
        holeInfo_p = []

        for s in self.surfaceList:
            s_p = s.transform(rotation=rot, translation=translation)
            surfaces_p.append(s_p)

        for hi in _np.array(self.holeInfo):
            pos_p = rotMat @ hi[0] + transMat
            dir_p = rotMat @ hi[1]
            holeInfo_p.append([pos_p, dir_p])

        block_p.surfaceList = surfaces_p
        block_p.geometry = block_p._makeGeometry(block_p.surfaceList)
        block_p.holeInfo = holeInfo_p

        return block_p

    def holeToHoleTransform(self):
        pass

    def addToRegistry(self, registry, replace=False):
        registry.addCell(self, replace=replace)

    def makeNewConnectBlock(self, localHole, foreignHole, angle, blockType, cellNumber=None, makeConnector=False):
        """
        localHole: where on THIS block to connect this block to EG. blockThis.getHolePosition(hole=3)
        foreignHole: where on ANOTHER block to connect this block to EG. blockOther.getHolePosition(hole=3)
        angle: rotate around the axis of connection, from zenith
        makeConnector : True
            returns [block, connector]
        makeConnector : False
            returns [block]
        """
        if cellNumber is None:
            cellNumber = self.cellNumber +1

        # work out transform from foreignHole to localHole
        h1Pos = _np.array(foreignHole[0])
        h1Vect = _np.array(foreignHole[1])
        h2Pos = _np.array(localHole[0])
        h2Vect = _np.array(localHole[1])

        a = h1Vect / _np.linalg.norm(h1Vect)
        b = h2Vect / _np.linalg.norm(h2Vect)
        axisIn = _np.cross(a, b)  # rotation axis
        axisInNorm = _np.linalg.norm(axisIn)
        dotProduct = _np.dot(a, b)
        angleRad = _np.arccos(dotProduct)
        angleDeg = _np.degrees(angleRad)
        axisIn = axisIn / axisInNorm
        print(f"{axisIn} {angleDeg}")

        # build block
        foreignBlock = Block(blockType=blockType, cellNumber=cellNumber)

        # transform block
        translation = h2Pos-h1Pos
        print(f"translation {translation}")

        # build connector
        if makeConnector is True:
            print(f"making connector")
            #connnector = Connector(localHole=, foreignHole=, translation=, rotation=, cellNumber=, reg=None)

        # transform connector
        #connector.transform()

        #return [block_p, connector]

    def addConnector(self, hole):
        connector = Connector(localHole=hole)
        return connector


class Connector(pyg4ometry.mcnp.Cell):
    """
    angle from zenith
    """
    def __init__(self, localHole, foreignHole=None, translation=[0, 0, 0], rotation=[0, 0, 0], cellNumber=None, reg=None):
        self.localHole = localHole
        self.foreignHole = foreignHole


        surface = pyg4ometry.mcnp.RCC(self.localHole, self.foreignHole, 0.3)

        self.transform(rotation=rotation, translation=translation)

        super().__init__(surfaces=surface, geometry=surface, cellNumber=cellNumber, reg=reg)

    def transform(self, rotation=[0, 0, 0], translation=[0, 0, 0]):
        rot = _np.eye(3)

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

            if rotation[0] > 0:
                # rotate in x
                theta = rotation[0] * _np.pi / 2
                rotX = _np.array([
                    [1.0, 0.0, 0.0],
                    [0.0, _np.cos(theta), -_np.sin(theta)],
                    [0.0, _np.sin(theta), _np.cos(theta)]
                ])
                rot = rotX @ rot
            if rotation[1] > 0:
                # rotate in y
                theta = rotation[1] * _np.pi / 2
                rotY = _np.array([
                    [_np.cos(theta), 0.0, _np.sin(theta)],
                    [0.0, 1.0, 0.0],
                    [-_np.sin(theta), 0.0, _np.cos(theta)]
                ])
                rot = rotY @ rot
            if rotation[2] > 0:
                # rotate in z
                theta = rotation[2] * _np.pi / 2
                rotZ = _np.array([
                    [_np.cos(theta), -_np.sin(theta), 0.0],
                    [_np.sin(theta), _np.cos(theta), 0.0],
                    [0.0, 0.0, 1.0]
                ])
                rot = rotZ @ rot

            # translation
            if len(translation) != 3 or not isinstance(translation, list):
                msg = "translation must be a list of length 3 for [x,y,z] translation"
                raise TypeError(msg)

        rotMat = _np.array(rot)
        transMat = _np.array(translation)

        # new block
        connector_p = Connector(blockType=self.blockType, cellNumber=self.cellNumber)  # copy block

        surfaces_new = []
        holeInfo_new = []

        surfaces_new = connector_p.surfaces.transform(rotation=rot, translation=translation)

        connector_p.surfaces = surfaces_new
        connector_p.geometry = connector_p._makeGeometry()

        return connector_p
