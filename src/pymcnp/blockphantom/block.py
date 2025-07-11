import pyg4ometry
import numpy as _np
from ..blockphantom import utils as _utils
from ..blockphantom import connector as _connector
import time

fullBlockDim = [11.0, 16.5, 5.5]  # dimensions of a full block
halfBlockDim = [11.0, 16.5, 2.5]  # dimensions of a half block

class Block(pyg4ometry.mcnp.Cell):
    fullBlockCache = None
    halfBlockCache = None
    def __init__(self, blockType, translation=[0, 0, 0], rotationSteps=[0, 0, 0], cellNumber=None, reg=None):
        self._mesh = None
        super().__init__(surfaces=[], reg=reg)  # a block is a cell

        self.blockType = blockType
        self.dim = fullBlockDim if blockType == "full" else halfBlockDim if blockType == "half" else None
        if self.dim is None:
            msg = f"Block type can only be 'full' or 'half'"
            raise TypeError(msg)
        self.small = 0.02  # to extend past the planes of the box by 0.01 cm so no inf small surface mesh covering hole
        self.unit = self.dim[1] / (3 * 2)  # hole separation unit on the surface of the block

        # define holes in local space
        self.holeInfo = self._defineHoles()
        # define surfaces in local space
        surfaces = self._makeSurfaces()
        # set cell surfaces
        self.addSurfaces(surfaces)
        # set cell geometry
        self.addGeometry(self._makeGeometry(self.surfaceList))

        # cash mesh if needed
        if blockType == "full" and Block.fullBlockCache is None:
            start_time = time.time()
            print("caching full block mesh...")
            Block.fullBlockCache = self.mesh()
            print("%s seconds to mesh full block" % (time.time() - start_time))
            print(" > cache complete")
        if blockType == "half" and Block.halfBlockCache is None:
            start_time = time.time()
            print("caching half block mesh...")
            Block.halfBlockCache = self.mesh()
            print("%s seconds to mesh half block" % (time.time() - start_time))
            print(" > cache complete")

        # copy cashed mesh for this block
        if blockType == "full":
            self._mesh = Block.fullBlockCache.clone()
        elif blockType == "half":
            self._mesh = Block.halfBlockCache.clone()

        # user inputted block transforms
        rotationMatrix = _utils.rotationStepsToMatrix(rotationSteps)
        translationVector = _np.array(translation)

        # apply transformations to holes to make them global space
        self.holeInfo = self._transformHoles(rotationMatrix, translationVector)

        holeNames = ["Bottom-Left", "Bottom-Right", "Top-Left", "Top-Right", "Left-Top", "Left-Middle",
                     "Left-Bottom", "Right-Top", "Right-Middle", "Right-Bottom", "Front-TopLeft", "Front-TopRight",
                     "Front-MiddleLeft", "Front-MiddleCenter", "Front-MiddleRight", "Front-BottomLeft",
                     "Front-BottomRight", "Back-TopLeft", "Back-TopRight", "Back-MiddleLeft", "Back-MiddleCenter",
                     "Back-MiddleRight", "Back-BottomLeft", "Back-BottomRight"]

        self.holeStatus = {i: {"name": holeNames[i], "connected": False, "covered": False, "hasConnector": False} for i in range(len(self.holeInfo))}

        # apply transformations to surfaces to make them global space
        surfaces_p = [s.transform(translation=translationVector.tolist(), rotation=rotationMatrix.tolist()) for s in surfaces]  # transformed surfaces

        # update mesh
        axis, angle = _utils.rotationMatrixToAxisAndAngle(rotationMatrix)
        self._mesh = self._mesh
        self._mesh.rotate(axis, 360-angle)
        self._mesh.translate(translationVector)

        # update geometry
        self.addGeometry(self._makeGeometry(self.surfaceList))

        if reg:
            # add s to registry and generate unique surfaceNumbers
            for s_p in surfaces_p:
                if s_p.surfaceNumber in reg.surfaceDict:
                    s_p.surfaceNumber = reg.getNewSurfaceNumber()
                if not s_p.surfaceNumber:
                    s_p.surfaceNumber = reg.getNewSurfaceNumber()
                reg.surfaceDict[s_p.surfaceNumber] = s_p
                self.addSurface(s_p)  # also add to the cell's surfaceList
        else:
            self.surfaceList = surfaces_p  # cell's surfaceList

        if reg:
            m1 = pyg4ometry.mcnp.Material(materialNumber=1, density=0.92, reg=reg)  # polyethylene
            self.addMaterial(m1)



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
        for hi in _np.array(self.holeInfo):
            position_p = rotationMatrix @ hi[0] + translationVector
            direction_p = rotationMatrix @ hi[1]
            holeInfo_p.append([position_p, direction_p])

        return holeInfo_p

    def _makeSurfaces(self):
        surfaces = [pyg4ometry.mcnp.PX((-self.dim[0] / 2), surfaceNumber=1),  # px1 (left)
                    pyg4ometry.mcnp.PX((self.dim[0] / 2), surfaceNumber=2),   # px2 (right)
                    pyg4ometry.mcnp.PY((-self.dim[1] / 2), surfaceNumber=3),  # py1 (bottom)
                    pyg4ometry.mcnp.PY((self.dim[1] / 2), surfaceNumber=4),   # py2 (top)
                    pyg4ometry.mcnp.PZ((-self.dim[2] / 2), surfaceNumber=5),  # pz1 (back)
                    pyg4ometry.mcnp.PZ((self.dim[2] / 2), surfaceNumber=6),   # pz2 (front)
                    ]
        i = len(surfaces)
        for holePosition, holeDirection in self.holeInfo:
            i = i + 1
            surfaces.append(pyg4ometry.mcnp.RCC(*holePosition, *holeDirection, 0.4, surfaceNumber=i))

        return surfaces

    def _makeGeometry(self, surfaces, reg=None):
        # polythene box
        geomBoxX = pyg4ometry.mcnp.Intersection(surfaces[0], pyg4ometry.mcnp.Complement(surfaces[1]))
        geomBoxY = pyg4ometry.mcnp.Intersection(surfaces[2], pyg4ometry.mcnp.Complement(surfaces[3]))
        geomBoxZ = pyg4ometry.mcnp.Intersection(surfaces[4], pyg4ometry.mcnp.Complement(surfaces[5]))
        geom = pyg4ometry.mcnp.Intersection(geomBoxZ, pyg4ometry.mcnp.Intersection(geomBoxY, geomBoxX))
        # connector holes
        for s in surfaces[6:]:
            #geom = pyg4ometry.mcnp.Intersection(geom, s)
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
            cellNumber=self.cellNumber
        )

        # transform surface
        surfaces_p = [s.transform(translation=translationVector.tolist(), rotation=rotationMatrix.tolist()) for s in self.surfaceList]

        # update the new block
        block_p.surfaceList = surfaces_p
        block_p.geometry = block_p._makeGeometry(surfaces_p)
        block_p.holeInfo = self._transformHoles(rotationMatrix, translationVector)
        block_p.holeStatus = self.holeStatus.copy()

        # update mesh
        axis, angle = _utils.rotationMatrixToAxisAndAngle(rotationMatrix)
        block_p._mesh = self._mesh
        block_p._mesh.rotate(axis, 360-angle)
        block_p._mesh.translate(translationVector)

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
        print(f"({localHole}) h1: {h1Position} , {h1Direction}")
        block_p = Block(blockType=newBlockType, cellNumber=cellNumber, reg=reg)
        h2Position, h2Direction = block_p.holeInfo[newBlockHole]
        print(f"({newBlockHole}) h2: {h2Position} , {h2Direction}")

        # calculate transformation (-h2 to h1)
        rotationMatrix = _utils.computeRotationMatrix(-_np.array(h2Direction), _np.array(h1Direction)) # negative hole 2 direction
        h2Position_rotated = rotationMatrix @ h2Position
        translationVector = h1Position - h2Position_rotated
        print(f"tr: {translationVector}")
        print(f"h2Pos_p: {h2Position_rotated}")

        # apply transformation to the new block
        block_p = block_p.transform(translation=translationVector.tolist(), rotation=rotationMatrix.tolist(), isRotationMatrix=True)

        # update connected hole status for new block
        self.holeStatus[localHole]["connected"] = True
        # update connected hole status for local block
        block_p.holeStatus[newBlockHole]["connected"] = True

        # update covered hole status for new block and local block
        # ToDo

        # make connector
        if makeConnector:
            if self.holeStatus[localHole]["hasConnector"]:
                msg = f"Hole {localHole} already has a connector"
                raise ValueError(msg)
            if cellNumber is not None:
                cellNumber = cellNumber+1
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
        rotationMatrix = _utils.computeRotationMatrix(_np.array([0, 0, 1]), _np.array(h1Direction))
        direction = h1Direction / _np.linalg.norm(h1Direction)
        translationVector = h1Position - (direction * (_connector.length / 2))

        # connector at the origin
        connector = _connector.Connector(
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

        holePosition, holeDirection = self.holeInfo[hole]

        # check if hole has connection to rotate around
        if self.holeStatus[hole]["connected"] is False:
            msg = f"hole {hole} has no connection to be rotated around."
            raise ValueError(msg)

        # check that the rotation given is around hole axis
        # the hole vector should be invariant for flipped if the rotation is allowed
        holeDirection = holeDirection / _np.linalg.norm(holeDirection)
        rotationMatrix = _utils.rotationStepsToMatrix(rotationSteps)
        holeDirection_rotated = rotationMatrix @ holeDirection
        if not (_np.allclose(holeDirection_rotated, holeDirection, atol=1e-6) or _np.allclose(holeDirection_rotated, -holeDirection, atol=1e-6)):
            msg = f"Rotation must be around the hole's axis of connection"
            raise ValueError(msg)

        # move block from connection so hole is at origin
        block_p = self.transform(translation=-holePosition, rotation=[0, 0, 0])
        # rotate by rotationSteps
        block_p = block_p.transform(translation=[0, 0, 0], rotation=rotationSteps)
        # move block back to connection
        block_p = block_p.transform(translation=holePosition, rotation=[0, 0, 0])

        # todo update hole status  - some holes will become uncovered and some covered

        return block_p

    def mesh(self):
        if self._mesh is not None:
            return self._mesh
        else:
            return super().mesh()
