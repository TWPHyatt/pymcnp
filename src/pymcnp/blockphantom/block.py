import pyg4ometry
import numpy as _np
import utils as _utils
import connector as _connector


class Block(pyg4ometry.mcnp.Cell):
    def __init__(self, blockType, translation=[0, 0, 0], rotationSteps=[0, 0, 0], cellNumber=None, reg=None):
        self.blockType = blockType
        self.dim = [11.0, 16.5, 5.5] if blockType == "full" else [11.0, 16.5, 2.5] if blockType == "half" else None
        if self.dim is None:
            msg = f"Block type can only be 'full' or 'half'"
            raise TypeError(msg)

        self.small = 0.02  # to extend past the planes of the box by 0.01 cm so no inf small surface mesh covering hole
        self.unit = self.dim[1] / (3 * 2)  # hole separation unit on the surface of the block

        # user defined transforms for block
        rotationMatrix = _utils.rotationStepsToMatrix(rotationSteps)
        translationVector = _np.array(translation)

        # define holes in local space
        self.holeInfo = self._defineHoles()

        # define surfaces in local space
        surfaces = self._makeSurfaces()

        # apply transformations to holes and surfaces to make them global space
        rotationMatrix = self._inputToRotationMatrix(rotationSteps)

        surfaces_p = [s.transform(rotation=rotationMatrix.tolist(), translation=translationVector) for s in surfaces]
        geometry = self._makeGeometry(surfaces_p)

        self.holeStatus = {i: {"connected": False, "covered": False, "hasConnector": False} for i in range(len(self.holeInfo))}

        # mcnp cell
        super().__init__(surfaces=surfaces_p, geometry=geometry, cellNumber=cellNumber, reg=reg)

        # if registry add block material to mcnp reg
        if reg:
            material = pyg4ometry.mcnp.Material(materialNumber=1, density=0.92)  # polyethylene
            reg.addMaterial(material)
            self.addMaterial(material)

    def printHoleInfo(self):
        msg = f""
        for i, el in enumerate(self.holeStatus):
            msg = f"{i} {el[i]['name']} : connected {el[i]['connected']} covered {el[i]['covered']} hasConnector {el[i]['hasConnector']} \n"
        return msg

    def _defineHoles(self):
        # define holes [xyz co-ordinates, xyz directions] relative to block center
        holes = [
            [[-self.unit, -(self.dim[1] + self.small) / 2, 0], [0, (self.dim[1] + self.small) / 2, 0]],  # bottom tubeL
            [[+self.unit, -(self.dim[1] + self.small) / 2, 0], [0, (self.dim[1] + self.small) / 2, 0]],  # bottom tubeR
            [[-self.unit, (self.dim[1] + self.small) / 2, 0], [0, -(self.dim[1] + self.small) / 2, 0]],  # top tubeL
            [[+self.unit, (self.dim[1] + self.small) / 2, 0], [0, -(self.dim[1] + self.small) / 2, 0]],  # top tubeR
            [[-self.dim[0] / 2, +self.unit * 2, 0], [1 + (self.small / 2), 0, 0]],  # holeL1
            [[-self.dim[0] / 2, 0, 0], [1 + (self.small / 2), 0, 0]],               # holeL2
            [[-self.dim[0] / 2, -self.unit * 2, 0], [1 + (self.small / 2), 0, 0]],  # holeL3
            [[self.dim[0] / 2, +self.unit * 2, 0], [-1 - (self.small / 2), 0, 0]],  # holeR1
            [[self.dim[0] / 2, 0, 0], [-1 - (self.small / 2), 0, 0]],               # holeR2
            [[self.dim[0] / 2, -self.unit * 2, 0], [-1 - (self.small / 2), 0, 0]],  # holeR3
            [[-self.unit, +self.unit * 2, (self.dim[2] + self.small) / 2], [0, 0, -1 - (self.small / 2)]],  # holeF1
            [[+self.unit, +self.unit * 2, (self.dim[2] + self.small) / 2], [0, 0, -1 - (self.small / 2)]],  # holeF2
            [[-self.unit, 0, (self.dim[2] + self.small) / 2], [0, 0, -1 - (self.small / 2)]],               # holeF3
            [[0, 0, (self.dim[2] + self.small) / 2], [0, 0, -1 - (self.small / 2)]],                        # holeF4
            [[+self.unit, 0, (self.dim[2] + self.small) / 2], [0, 0, -1 - (self.small / 2)]],               # holeF5
            [[-self.unit, -self.unit * 2, (self.dim[2] + self.small) / 2], [0, 0, -1 - (self.small / 2)]],  # holeF6
            [[+self.unit, -self.unit * 2, (self.dim[2] + self.small) / 2], [0, 0, -1 - (self.small / 2)]],  # holeF7
        ]

        if self.blockType == "full":
            holesBack = [
                [[-self.unit, +self.unit * 2, -(self.dim[2] + self.small) / 2], [0, 0, 1 + (self.small / 2)]],  # holeB1
                [[+self.unit, +self.unit * 2, -(self.dim[2] + self.small) / 2], [0, 0, 1 + (self.small / 2)]],  # holeB2
                [[-self.unit, 0, -(self.dim[2] + self.small) / 2], [0, 0, 1 + (self.small / 2)]],               # holeB3
                [[0, 0, -(self.dim[2] + self.small) / 2], [0, 0, 1 + (self.small / 2)]],                        # holeB4
                [[+self.unit, 0, -(self.dim[2] + self.small) / 2], [0, 0, 1 + (self.small / 2)]],               # holeB5
                [[-self.unit, -self.unit * 2, -(self.dim[2] + self.small) / 2], [0, 0, 1 + (self.small / 2)]],  # holeB6
                [[+self.unit, -self.unit * 2, -(self.dim[2] + self.small) / 2], [0, 0, 1 + (self.small / 2)]],  # holeB7
            ]
            holes.extend(holesBack)

        return holes

    def _makeSurfaces(self):
        surfaces = [pyg4ometry.mcnp.PX((-self.dim[0] / 2)),  # px1 (left)
                    pyg4ometry.mcnp.PX((self.dim[0] / 2)),   # px2 (right)
                    pyg4ometry.mcnp.PY((-self.dim[1] / 2)),  # py1 (bottom)
                    pyg4ometry.mcnp.PY((self.dim[1] / 2)),   # py2 (top)
                    pyg4ometry.mcnp.PZ((-self.dim[2] / 2)),  # pz1 (back)
                    pyg4ometry.mcnp.PZ((self.dim[2] / 2)),   # pz2 (front)
        ]
        for holePosition, holeDirection in self.holeInfo:
            surfaces.append(pyg4ometry.mcnp.RCC(*holePosition, *holeDirection, 0.4))
        return surfaces

    def _makeGeometry(self, surfaces):
        # polythene box
        geomBoxX = pyg4ometry.mcnp.Intersection(surfaces[0], pyg4ometry.mcnp.Complement(surfaces[1]))
        geomBoxY = pyg4ometry.mcnp.Intersection(surfaces[2], pyg4ometry.mcnp.Complement(surfaces[3]))
        geomBoxZ = pyg4ometry.mcnp.Intersection(surfaces[4], pyg4ometry.mcnp.Complement(surfaces[5]))
        geomBlock = pyg4ometry.mcnp.Intersection(geomBoxZ, pyg4ometry.mcnp.Intersection(geomBoxY, geomBoxX))

        for i in range(6, 20 + 1):
            # connector holes tube, left, right, front
            geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[i]))

        if self.blockType == "full":
            for i in range(21, 27 + 1):
                # connector holes back
                geomBlock = pyg4ometry.mcnp.Intersection(geomBlock, pyg4ometry.mcnp.Complement(surfaces[i]))

        return geomBlock

    def getHole(self, holeNumber):
        if not (0 <= holeNumber < len(self.holeInfo)):
            msg = f"Invalid hole number."
            raise TypeError(msg)
        return self.holeInfo[holeNumber]

    def _computeRotationMatrix(self, v1, v2):
        v1 = v1 / _np.linalg.norm(v1)
        v2 = v2 / _np.linalg.norm(v2)

        crossProd = _np.cross(v1, v2)
        dotProd = _np.dot(v1, v2)

        if _np.allclose(crossProd, 0):  # vectors are parallel or opposite direction
            if dotProd > 0:
                return _np.eye(3)  # no rotation
            else:
                return -_np.eye(3)  # opposite direction

        crossProdNorm = _np.linalg.norm(crossProd)
        crossProd = crossProd / crossProdNorm  # normalise axis

        angle = _np.arccos(dotProd)

        # rodrigues' formula where k is matrix of cross products
        K = _np.array([
            [0, -crossProd[2], crossProd[1]],
            [crossProd[2], 0, -crossProd[0]],
            [-crossProd[1], crossProd[0], 0]
        ])

        return _np.eye(3) + _np.sin(angle) * K + (1 - _np.cos(angle)) * _np.dot(K, K)

    def transform(self, rotationMatrix, translation):
        rotMat = _np.array(rotationMatrix)
        transMat = _np.array(translation)

        block_p = Block(blockType=self.blockType, cellNumber=self.cellNumber)

        surfaces_p = [s.transform(rotation=rotMat.tolist(), translation=transMat.tolist()) for s in self.surfaceList]

        holeInfo_p = []
        for hi in _np.array(self.holeInfo):
            pos_p = rotMat @ hi[0] + transMat
            dir_p = rotMat @ hi[1]
            holeInfo_p.append([pos_p, dir_p])

        block_p.surfaceList = surfaces_p
        block_p.geometry = block_p._makeGeometry(surfaces_p)
        block_p.holeInfo = holeInfo_p

        return block_p

    def _isPointInsideBlock(self, dimensions, rotation, translation, point):
        """
        calculates the min and max corners of a block (coordinates post transformation) and checks if point is inside the box given by min_corner and max_corner
        """
        cornerCoords = _np.array([
            [0, 0, 0],              # origin corner (min x,y,z)
            [dimensions[0], 0, 0],  # corner in the +x
            [0, dimensions[1], 0],  # corner in the +y
            [0, 0, dimensions[2]],  # corner in the +z
            [dimensions[0], dimensions[1], 0],  # corner in the +x +y
            [dimensions[0], 0, dimensions[2]],  # corner in the +x +z
            [0, dimensions[1], dimensions[2]],  # corner in the +y +z
            [dimensions[0], dimensions[1], dimensions[2]]  # corner in the +x +y +z (max x,y,z)
        ])
        cornerCoords_p = rotation @ cornerCoords.T  # apply rotation to each corner
        cornerCoords_p = cornerCoords_p.T + translation  # apply translation to each corner (now global coordinates)
        minCorner = _np.min(cornerCoords_p, axis=0)
        maxCorner = _np.max(cornerCoords_p, axis=0)
        if _np.all(point >= minCorner) and _np.all(point <= maxCorner):
            return True
        return False

    def makeNewConnectedBlock(self, localHole, newBlockHole, newBlockType, cellNumber=None, makeConnector=False, reg=None):
        if not (0 <= localHole < len(self.holeInfo)) or not (0 <= newBlockHole < len(self.holeInfo)):
            msg = f"Hole numbers must be between 0 and {len(self.holeInfo) - 1}."
            raise TypeError(msg)

        # check if hole is available
        if self.holeStatus[localHole]["connected"] or self.holeStatus[localHole]["covered"]:
            msg = f"Hole {localHole} is not available for connection."
            raise ValueError(msg)

        h1Pos = _np.array(self.getHole(localHole)[0])
        h1Vect = _np.array(self.getHole(localHole)[1])

        block_p = Block(blockType=newBlockType, cellNumber=cellNumber)
        h2Pos = _np.array(block_p.getHole(newBlockHole)[0])
        h2Vect = -_np.array(block_p.getHole(newBlockHole)[1])  # negative for opposite direction of holes

        # transformation
        rotMat = self._computeRotationMatrix(h2Vect, h1Vect)
        h2Pos_p = rotMat @ h2Pos  # rotated h2Pos
        trans = h1Pos - h2Pos_p

        # transform the new block
        block_p = block_p.transform(rotationMatrix=rotMat.tolist(), translation=trans.tolist())

        # update hole status as connected
        self.holeStatus[localHole]["connected"] = True
        block_p.holeStatus[newBlockHole]["connected"] = True

        # update hole status as covered for overlapped holes
        for holeNum, hole in enumerate(self.holeInfo):
            if self._isPointInsideBlock(block_p.dim, rotMat, trans, hole[0]):
                self.holeStatus[holeNum]["covered"] = True  # update hole status for holes under block as 'covered'

        # make connector
        if makeConnector:
            if self.holeStatus[localHole]["hasConnector"]:
                msg = f"Hole {localHole} already has a connector."
                raise ValueError(msg)
            connector = self.addConnector(localHole, cellNumber, reg)
            self.holeStatus[localHole]["hasConnector"] = True
            block_p.holeStatus[newBlockHole]["hasConnector"] = True
            return [block_p, connector]

        return block_p

    def addConnector(self, localHole, cellNumber=None, reg=None):
        if not (0 <= localHole < len(self.holeInfo)):
            msg = f"Hole numbers must be between 0 and {len(self.holeInfo) - 1}."
            raise TypeError(msg)

        if self.holeStatus[localHole]["hasConnector"]:
            msg = f"Hole {localHole} already has connector."
            raise ValueError(msg)

        h1Pos = _np.array(self.getHole(localHole)[0])
        h1Vect = _np.array(self.getHole(localHole)[1])

        # connector creation
        connectorLength = 1.5
        rotMat = self._computeRotationMatrix(_np.array([0, 0, 1]), h1Vect)
        direction = h1Vect / _np.linalg.norm(h1Vect)
        translation = h1Pos - direction * (connectorLength / 2)
        connector = Connector(
            translation=translation.tolist(),
            rotation=rotMat.tolist(),
            length=connectorLength,
            cellNumber=cellNumber,
            reg=reg
        )
        return connector

    def rotateAboutConnection(self, hole, rotation):
        """
        un-transform block so hole is at origin and then apply rotation and translate back
        """

        if self.holeStatus[hole]["connected"] is False:
            msg = f"hole {hole} has no connection to be rotated around."
            raise ValueError(msg)

        holePos = _np.array(self.getHole(hole)[0])
        holeVect = _np.array(self.getHole(hole)[1])
        holeVect = holeVect / _np.linalg.norm(holeVect)

        rotMatIn = self._inputToRotationMatrix(rotation)

        # check that the rotation given is around hole axis
        # the hole vector should be invarient for flipped if the rotation is allowed
        holeVect_p = rotMatIn @ holeVect
        if not (_np.allclose(holeVect_p, holeVect, atol=1e-6) or _np.allclose(holeVect_p, -holeVect, atol=1e-6)):
            msg = f"Rotation must be around the hole's axis of connection."
            raise ValueError(msg)

        # apply rotation around the hole
        transMat = holePos - rotMatIn @ holePos
        block_p = self.transform(rotationMatrix=rotMatIn.tolist(), translation=transMat.tolist())

        return block_p
