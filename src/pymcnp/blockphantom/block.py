import pyg4ometry
import numpy as _np
from ..blockphantom import utils as _utils
from ..blockphantom import connector as _connector


class Block(pyg4ometry.mcnp.Cell):
    def __init__(self, blockType, translation=[0, 0, 0], rotationSteps=[0, 0, 0], cellNumber=None, reg=None):
        self.blockType = blockType
        self.dim = [11.0, 16.5, 5.5] if blockType == "full" else [11.0, 16.5, 2.5] if blockType == "half" else None
        if self.dim is None:
            msg = f"Block type can only be 'full' or 'half'"
            raise TypeError(msg)

        self.small = 0.02  # to extend past the planes of the box by 0.01 cm so no inf small surface mesh covering hole
        self.unit = self.dim[1] / (3 * 2)  # hole separation unit on the surface of the block

        # user inputted block transforms
        rotationMatrix = _utils.rotationStepsToMatrix(rotationSteps)
        translationVector = _np.array(translation)

        # define holes and surfaces in local space
        self.holeInfo = self._defineHoles()
        surfaces = self._makeSurfaces()

        # apply transformations to holes and surfaces to make them global space
        self.holeInfo = self._transformHoles(rotationMatrix, translationVector)
        surfaces_p = [s.transform(translation=translationVector.tolist(), rotation=rotationMatrix.tolist()) for s in surfaces]
        geometry = self._makeGeometry(surfaces_p)

        holeNames = ["Bottom-Left", "Bottom-Right", "Top-Left", "Top-Right", "Left-Top", "Left-Middle",
                     "Left-Bottom", "Right-Top", "Right-Middle", "Right-Bottom", "Front-TopLeft", "Front-TopRight",
                     "Front-MiddleLeft", "Front-MiddleCenter", "Front-MiddleRight", "Front-BottomLeft",
                     "Front-BottomRight", "Back-TopLeft", "Back-TopRight", "Back-MiddleLeft", "Back-MiddleCenter",
                     "Back-MiddleRight", "Back-BottomLeft", "Back-BottomRight"]
        self.holeStatus = {i: {"name": holeNames[i], "connected": False, "covered": False, "hasConnector": False} for i in range(len(self.holeInfo))}

        # a block is a cell
        super().__init__(surfaces=surfaces_p, geometry=geometry, cellNumber=cellNumber, reg=reg)

        # add block material to registry
        if reg:
            material = pyg4ometry.mcnp.Material(materialNumber=1, density=0.92)  # polyethylene
            reg.addMaterial(material)
            self.addMaterial(material)

    def printHoleInfo(self):
        msg = f""
        for i in range(0, len(self.holeStatus)):
            print(f"{i} {self.holeStatus[i]['name']} : connected {self.holeStatus[i]['connected']} covered {self.holeStatus[i]['covered']} hasConnector {self.holeStatus[i]['hasConnector']}")
        return

    def _defineHoles(self):
        # define holes [xyz co-ordinates, xyz directions] relative to block center
        holes = [
            [[-self.unit, -(self.dim[1] + self.small) / 2, 0], [0, self.dim[1]/2 + self.small, 0]],  # bottom tubeL
            [[+self.unit, -(self.dim[1] + self.small) / 2, 0], [0, self.dim[1]/2 + self.small, 0]],  # bottom tubeR
            [[-self.unit, (self.dim[1] + self.small) / 2, 0], [0, -self.dim[1]/2 - self.small, 0]],  # top tubeL
            [[+self.unit, (self.dim[1] + self.small) / 2, 0], [0, -self.dim[1]/2 - self.small, 0]],  # top tubeR
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

    def _transformHoles(self, rotationMatrix, translationVector):
        holeInfo_p = []
        for position, direction in self.holeInfo:
            position_p = rotationMatrix @ position + translationVector
            direction_p = rotationMatrix @ direction
            holeInfo_p.append([position_p, direction_p])

        return holeInfo_p

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
        geom = pyg4ometry.mcnp.Intersection(geomBoxZ, pyg4ometry.mcnp.Intersection(geomBoxY, geomBoxX))
        # connector holes
        for s in surfaces[6:]:
            geom = pyg4ometry.mcnp.Intersection(geom, pyg4ometry.mcnp.Complement(s))

        return geom

    def transform(self, translation=[0, 0, 0], rotation=[0, 0, 0], isRotationMatrix=False):
        if isRotationMatrix:
            rotationMatrix = _np.array(rotation)
        else:
            rotationMatrix = _utils.rotationStepsToMatrix(rotation)
        translationVector = _np.array(translation)

        if hasattr(self, 'reg'):
            reg = self.reg
        else:
            reg = None

        # new block (prime)
        block_p = Block(
            blockType=self.blockType,
            cellNumber=self.cellNumber,
            reg=reg
        )

        # transform surface
        surfaces_p = [s.transform(translation=translationVector.tolist(), rotation=rotationMatrix.tolist()) for s in self.surfaceList]

        # update the new block
        block_p.surfaceList = surfaces_p
        block_p.geometry = block_p._makeGeometry(surfaces_p)
        block_p.holeInfo = block_p._transformHoles(rotationMatrix, translationVector)
        block_p.holeStatus = self.holeStatus.copy()

        return block_p

    def _isPointInsideBlock(self, dimensions, rotationMatrix, translationVector, pointXYZ):
        """
        calculates the min and max corners of a block (coordinates post transformation)
        and checks if point is inside the box given by minCorner and maxCorner
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
        cornerCoords_p = rotationMatrix @ cornerCoords.T  # apply rotation to each corner
        cornerCoords_p = cornerCoords_p.T + translationVector  # apply translation to each corner (now global coords)
        minCorner = _np.min(cornerCoords_p, axis=0)
        maxCorner = _np.max(cornerCoords_p, axis=0)
        if _np.all(pointXYZ >= minCorner) and _np.all(pointXYZ <= maxCorner):
            return True
        else:
            return False

    def makeNewConnectedBlock(self, newBlockType, newBlockHole, localHole, cellNumber=None, makeConnector=False, reg=None):
        # check input is correct
        if not (0 <= localHole < len(self.holeInfo)) or not (0 <= newBlockHole < len(self.holeInfo)):
            msg = f"Hole numbers must be between 0 and {len(self.holeInfo) - 1}"
            raise TypeError(msg)

        # check if hole is available
        if self.holeStatus[localHole]["connected"] or self.holeStatus[localHole]["covered"]:
            msg = f"Hole {localHole} is not available for connection"
            raise ValueError(msg)

        h1Position, h1Direction = self.holeInfo[localHole]
        block_p = Block(blockType=newBlockType, cellNumber=cellNumber)
        h2Position, h2Direction = block_p.holeInfo[newBlockHole]

        # calculate transformation (-h2 to h1)
        rotationMatrix = _utils.computeRotationMatrix(-_np.array(h2Direction), _np.array(h1Direction)) # negative hole 2 direction
        h2Position_rotated = rotationMatrix @ h2Position
        translationVector = h1Position - h2Position_rotated

        # apply transformation to the new block
        block_p = block_p.transform(translation=translationVector.tolist(), rotation=rotationMatrix.tolist(), isRotationMatrix=True)

        # update hole status as connected
        self.holeStatus[localHole]["connected"] = True
        block_p.holeStatus[newBlockHole]["connected"] = True

        # update hole status as covered for overlapped holes
        for holeNum, hole in enumerate(self.holeInfo):
            if self._isPointInsideBlock(block_p.dim, rotationMatrix, translationVector, hole[0]):
                self.holeStatus[holeNum]["covered"] = True  # update hole status for holes under block as 'covered'

        # make connector
        if makeConnector:
            if self.holeStatus[localHole]["hasConnector"]:
                msg = f"Hole {localHole} already has a connector"
                raise ValueError(msg)
            connector = self.addConnector(localHole, cellNumber, reg)
            self.holeStatus[localHole]["hasConnector"] = True
            block_p.holeStatus[newBlockHole]["hasConnector"] = True
            return [block_p, connector]

        return block_p

    def addConnector(self, localHole, cellNumber=None, reg=None):
        # check input is correct
        if not (0 <= localHole < len(self.holeInfo)):
            msg = f"Hole numbers must be between 0 and {len(self.holeInfo) - 1}"
            raise TypeError(msg)

        # check if hole is available for connector
        if self.holeStatus[localHole]["hasConnector"]:
            msg = f"Hole {localHole} already has connector"
            raise ValueError(msg)

        h1Position, h1Direction = self.holeInfo[localHole]

        # connector creation
        length = 1.5
        rotationMatrix = _utils.computeRotationMatrix(_np.array([0, 0, 1]), _np.array(h1Direction))
        direction = h1Direction / _np.linalg.norm(h1Direction)
        translationVector = h1Position - (direction * (length / 2))

        # connector at the origin
        connector = _connector.Connector(
            length=length,
            cellNumber=cellNumber,
            reg=reg
        )
        # apply transformation to move and align to hole, with 50% length outside hole
        connector_p = connector.transform(translation=translationVector, rotation=rotationMatrix, isRotationMatrix=True)

        return connector_p

    def rotateAboutConnection(self, hole, rotationSteps):
        """
        un-transform block so hole is at origin and then apply rotation and translate back
        """

        # check if hole has connection to rotate around
        if self.holeStatus[hole]["connected"] is False:
            msg = f"hole {hole} has no connection to be rotated around."
            raise ValueError(msg)

        holePosition, holeDirection = self.holeInfo[hole]
        holeDirection = holeDirection / _np.linalg.norm(holeDirection)

        rotationMatrix = _utils.rotationStepsToMatrix(rotationSteps)
        holeDirection_rotated = rotationMatrix @ holeDirection

        # check that the rotation given is around hole axis
        # the hole vector should be invariant for flipped if the rotation is allowed
        if not (_np.allclose(holeDirection_rotated, holeDirection, atol=1e-6) or _np.allclose(holeDirection_rotated, -holeDirection, atol=1e-6)):
            msg = f"Rotation must be around the hole's axis of connection"
            raise ValueError(msg)

        # apply rotation around the hole
        translationVector = holePosition - rotationMatrix @ holePosition
        block_p = self.transform(translation=translationVector.tolist(), rotation=rotationMatrix.tolist(), isRotationMatrix=True)

        return block_p
