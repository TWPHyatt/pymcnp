import pyg4ometry
import numpy as _np

class Block(pyg4ometry.mcnp.Cell):
    def __init__(self, blockType, translation=[0, 0, 0], rotation=[0, 0, 0], cellNumber=None, reg=None):
        self.blockType = blockType
        self.dim = [11.0, 16.5, 5.5] if blockType == "full" else [11.0, 16.5, 2.5] if blockType == "half" else None
        if self.dim is None:
            raise TypeError("Block type can only be 'full' or 'half'")

        self.small = 0.02  # to extend past the planes of the box by 0.01 cm so no inf small surface mesh covering hole
        self.unit = self.dim[1] / (3 * 2)  # hole separation unit on the surface of the block
        # HoleInfo is the [xyz co-ordinates, xyz axis-vector] for each hole FROM THE CENTER OF EACH BLOCK
        self.holeInfo = [[[-self.unit, -(self.dim[1] + self.small) / 2, 0], [0, (self.dim[1] + self.small) / 2, 0]],  # bottom tubeL
                         [[+self.unit, -(self.dim[1] + self.small) / 2, 0], [0, (self.dim[1] + self.small) / 2, 0]],  # bottom tubeR
                         [[-self.unit, (self.dim[1] + self.small) / 2, 0], [0, -(self.dim[1] + self.small) / 2, 0]],  # top tubeL
                         [[+self.unit, (self.dim[1] + self.small) / 2, 0], [0, -(self.dim[1] + self.small) / 2, 0]],  # top tubeR
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

        self.holeStatus = {
            0: {"connected": False, "covered": False, "hasConnector": False},
            1: {"connected": False, "covered": False, "hasConnector": False},
            2: {"connected": False, "covered": False, "hasConnector": False},
            3: {"connected": False, "covered": False, "hasConnector": False},
            4: {"connected": False, "covered": False, "hasConnector": False},
            5: {"connected": False, "covered": False, "hasConnector": False},
            6: {"connected": False, "covered": False, "hasConnector": False},
            7: {"connected": False, "covered": False, "hasConnector": False},
            8: {"connected": False, "covered": False, "hasConnector": False},
            9: {"connected": False, "covered": False, "hasConnector": False},
            10: {"connected": False, "covered": False, "hasConnector": False},
            11: {"connected": False, "covered": False, "hasConnector": False},
            12: {"connected": False, "covered": False, "hasConnector": False},
            13: {"connected": False, "covered": False, "hasConnector": False},
            14: {"connected": False, "covered": False, "hasConnector": False},
            15: {"connected": False, "covered": False, "hasConnector": False},
            16: {"connected": False, "covered": False, "hasConnector": False},
            17: {"connected": False, "covered": False, "hasConnector": False},
            18: {"connected": False, "covered": False, "hasConnector": False},
            19: {"connected": False, "covered": False, "hasConnector": False},
            20: {"connected": False, "covered": False, "hasConnector": False},
            21: {"connected": False, "covered": False, "hasConnector": False},
            22: {"connected": False, "covered": False, "hasConnector": False},
            23: {"connected": False, "covered": False, "hasConnector": False}
        }

        surfaces = self._makeSurfaces()

        rotMat = self._inputToRotationMatrix(rotation)
        transMat = _np.array(translation)

        surfaces_p = [s.transform(rotation=rotMat.tolist(), translation=translation) for s in surfaces]
        geometry = self._makeGeometry(surfaces_p)

        holeInfo_new = []
        for hi in _np.array(self.holeInfo):
            pos_new = rotMat @ hi[0] + transMat
            dir_new = rotMat @ hi[1]
            holeInfo_new.append([pos_new, dir_new])
        self.holeInfo = holeInfo_new

        super().__init__(surfaces=surfaces_p, geometry=geometry, cellNumber=cellNumber, reg=reg)

    def _inputToRotationMatrix(self, rotation):
        """
        Converts a list of 90-degree step rotations [x, y, z] into a 3x3 rotation matrix
        """
        if len(rotation) != 3 or not all(isinstance(i, int) for i in rotation):
            raise TypeError("rotation must be a list of 3 integers [x, y, z] representing 90° steps")

        rotMat = _np.eye(3)

        # Rotation around x-axis
        if rotation[0] != 0:
            theta = rotation[0] * _np.pi / 2
            rotX = _np.array([
                [1, 0, 0],
                [0, _np.cos(theta), -_np.sin(theta)],
                [0, _np.sin(theta), _np.cos(theta)]
            ])
            rotMat = rotX @ rotMat

        # Rotation around y-axis
        if rotation[1] != 0:
            theta = rotation[1] * _np.pi / 2
            rotY = _np.array([
                [_np.cos(theta), 0, _np.sin(theta)],
                [0, 1, 0],
                [-_np.sin(theta), 0, _np.cos(theta)]
            ])
            rotMat = rotY @ rotMat

        # Rotation around z-axis
        if rotation[2] != 0:
            theta = rotation[2] * _np.pi / 2
            rotZ = _np.array([
                [_np.cos(theta), -_np.sin(theta), 0],
                [_np.sin(theta), _np.cos(theta), 0],
                [0, 0, 1]
            ])
            rotMat = rotZ @ rotMat

        return rotMat

    def _makeSurfaces(self):
        surfaces = [pyg4ometry.mcnp.PX((-self.dim[0] / 2)),  # px1 (left)
                    pyg4ometry.mcnp.PX((self.dim[0] / 2)),   # px2 (right)
                    pyg4ometry.mcnp.PY((-self.dim[1] / 2)),  # py1 (bottom)
                    pyg4ometry.mcnp.PY((self.dim[1] / 2)),   # py2 (top)
                    pyg4ometry.mcnp.PZ((-self.dim[2] / 2)),  # pz1 (back)
                    pyg4ometry.mcnp.PZ((self.dim[2] / 2)),   # pz2 (front)
                    pyg4ometry.mcnp.RCC(*self.holeInfo[0][0], *self.holeInfo[0][1], 0.4),  # bottom-top left tube
                    pyg4ometry.mcnp.RCC(*self.holeInfo[1][0], *self.holeInfo[1][1], 0.4),  # bottom-top right tube
                    pyg4ometry.mcnp.RCC(*self.holeInfo[2][0], *self.holeInfo[2][1], 0.4),  # left-side top
                    pyg4ometry.mcnp.RCC(*self.holeInfo[3][0], *self.holeInfo[3][1], 0.4),  # left-side middle
                    pyg4ometry.mcnp.RCC(*self.holeInfo[4][0], *self.holeInfo[4][1], 0.4),  # left-side bottom
                    pyg4ometry.mcnp.RCC(*self.holeInfo[5][0], *self.holeInfo[5][1], 0.4),  # right-side top
                    pyg4ometry.mcnp.RCC(*self.holeInfo[6][0], *self.holeInfo[6][1], 0.4),  # right-side middle
                    pyg4ometry.mcnp.RCC(*self.holeInfo[7][0], *self.holeInfo[7][1], 0.4),  # right-side bottom
                    pyg4ometry.mcnp.RCC(*self.holeInfo[8][0], *self.holeInfo[8][1], 0.4),  # front-side top left
                    pyg4ometry.mcnp.RCC(*self.holeInfo[9][0], *self.holeInfo[9][1], 0.4),  # front-side top right
                    pyg4ometry.mcnp.RCC(*self.holeInfo[10][0], *self.holeInfo[10][1], 0.4),  # front-side middle left
                    pyg4ometry.mcnp.RCC(*self.holeInfo[11][0], *self.holeInfo[11][1], 0.4),  # front-side middle cent
                    pyg4ometry.mcnp.RCC(*self.holeInfo[12][0], *self.holeInfo[12][1], 0.4),  # front-side middle right
                    pyg4ometry.mcnp.RCC(*self.holeInfo[13][0], *self.holeInfo[13][1], 0.4),  # front-side bottom left
                    pyg4ometry.mcnp.RCC(*self.holeInfo[14][0], *self.holeInfo[14][1], 0.4)]  # front-side bottom right


        if self.blockType == "full":
            # connector holes back-side
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
            msg = "Invalid hole number."
            raise TypeError(msg)
        return self.holeInfo[holeNumber]

    def _computeRotationMatrix(self, v1, v2):
        v1 = v1 / _np.linalg.norm(v1)  # Normalize v1
        v2 = v2 / _np.linalg.norm(v2)  # Normalize v2

        cross_prod = _np.cross(v1, v2)  # Axis of rotation
        dot_prod = _np.dot(v1, v2)  # Cosine of angle between vectors

        if _np.allclose(cross_prod, 0):  # Vectors are parallel or anti-parallel
            if dot_prod > 0:
                return _np.eye(3)  # No rotation needed
            else:
                return -_np.eye(3)  # 180° rotation (opposite direction)

        cross_prod_norm = _np.linalg.norm(cross_prod)
        cross_prod = cross_prod / cross_prod_norm  # Normalize axis

        angle = _np.arccos(dot_prod)  # Angle between vectors

        # Skew-symmetric matrix for cross product
        K = _np.array([
            [0, -cross_prod[2], cross_prod[1]],
            [cross_prod[2], 0, -cross_prod[0]],
            [-cross_prod[1], cross_prod[0], 0]
        ])

        return _np.eye(3) + _np.sin(angle) * K + (1 - _np.cos(angle)) * _np.dot(K, K)

    def transform(self, rotation_matrix, translation):
        import numpy as np
        rotMat = np.array(rotation_matrix)
        transMat = np.array(translation)

        block_p = Block(blockType=self.blockType, cellNumber=self.cellNumber)

        surfaces_p = [s.transform(rotation=rotMat.tolist(), translation=transMat.tolist()) for s in self.surfaceList]

        holeInfo_p = []
        for hi in np.array(self.holeInfo):
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

    def makeNewConnectedBlock(self, localHole, newBlockHole, newBlockType, cellNumber=None, makeConnector=False):
        if not (0 <= localHole < len(self.holeInfo)) or not (0 <= newBlockHole < len(self.holeInfo)):
            raise TypeError(f"Hole numbers must be between 0 and {len(self.holeInfo) - 1}.")

        if self.holeStatus[localHole]["connected"] or self.holeStatus[localHole]["covered"]:
            raise ValueError(f"Hole {localHole} is not available for connection.")

        if cellNumber is None:
            cellNumber = self.cellNumber + 1

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
        block_p = block_p.transform(rotation_matrix=rotMat.tolist(), translation=trans.tolist())

        # update hole status as connected
        self.holeStatus[localHole]["connected"] = True
        block_p.holeStatus[newBlockHole]["connected"] = True

        for holeNum, hole in enumerate(self.holeInfo):
            if self.holeStatus[holeNum]["connected"] is True:
                continue  # if already connected skip this hole
            elif self._isPointInsideBlock(block_p.dim, rotMat, trans, hole[0]):
                self.holeStatus[holeNum]["covered"] = True  # update hole status for holes under block as 'covered'

        if makeConnector:
            if self.holeStatus[localHole]["hasConnector"]:
                raise ValueError(f"Hole {localHole} already has a connector.")
            self.holeStatus[localHole]["hasConnector"] = True
            block_p.holeStatus[newBlockHole]["hasConnector"] = True

            # connector creation
            connector = Connector(
                translation=trans.tolist(),
                rotation=rotMat.tolist(),
                cellNumber=cellNumber,
                reg=None
            )
            return [block_p, connector]
        return block_p

    def addConnector(self, localHole, foreignBlock, foreignBlockHole, cellNumber=None, reg=None):
        if not (0 <= localHole < len(self.holeInfo)) or not (0 <= foreignBlockHole < len(self.holeInfo)):
            raise TypeError(f"Hole numbers must be between 0 and {len(self.holeInfo) - 1}.")

        if self.holeStatus[localHole]["connected"] or self.holeStatus[localHole]["covered"]:
            raise ValueError(f"Hole {localHole} is not available for connection.")

        if cellNumber is None:
            cellNumber = self.cellNumber + 1

        h1Pos = _np.array(self.getHole(localHole)[0])
        h1Vect = _np.array(self.getHole(localHole)[1])

        h2Pos = _np.array(foreignBlock.getHole(foreignBlockHole)[0])
        h2Vect = -_np.array(foreignBlock.getHole(foreignBlockHole)[1])  # negative for opposite direction of holes

        if h1Vect != h2Vect:
            msg = f"hole {localHole} and {foreignBlock} are not aligned for a connector"

        # connector creation
        connector = Connector(
            translation=trans.tolist(),
            rotation=rotMat.tolist(),
            cellNumber=cellNumber,
            reg=None
        )
        return connector

    def addIsolatedConnector(self):

        # connector creation
        connector = Connector(
            translation=trans.tolist(),
            rotation=rotMat.tolist(),
            cellNumber=cellNumber,
            reg=None
        )
        return connector




class Connector(pyg4ometry.mcnp.Cell):
    def __init__(self, translation=[0, 0, 0], rotation=_np.eye(3), cellNumber=None,
                 reg=None):
        self.position =
        self.direction =

        surface = pyg4ometry.mcnp.RCC(self.localHole, self.newBlockHole, 0.3)

        rotMat = self._inputToRotationMatrix(rotation)
        surface_p = surface.transform(rotation=rotMat.tolist(), translation=translation)

        super().__init__(surfaces=surface_p, geometry=surface_p, cellNumber=cellNumber, reg=reg)

