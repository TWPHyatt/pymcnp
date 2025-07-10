import pyg4ometry
import numpy as _np
from ..blockphantom import utils as _utils

length = 1.5
radius = 0.3

class Connector(pyg4ometry.mcnp.Cell):
    connectorCache = None
    def __init__(self, translation=[0, 0, 0], rotationSteps=[0, 0, 0], cellNumber=None, reg=None):
        self._mesh = None
        super().__init__(surfaces=[], reg=reg)  # a connector is a cell

        rotationMatrix = _utils.rotationStepsToMatrix(rotationSteps)
        translationVector = _np.array(translation)

        surface = pyg4ometry.mcnp.RCC(0, 0, 0, 0, 0, length, radius)

        if Connector.connectorCache is None:
            print("caching connector mesh...")
            Connector.connectorCache = self.mesh()
            print(" > cache complete")

        self._mesh = Connector.connectorCache.clone()

        surface_p = surface.transform(translation=translationVector.tolist(), rotation=rotationMatrix.tolist())

        if reg:
            # add s to registry and generate unique surfaceNumbers
            if surface_p.surfaceNumber in reg.surfaceDict:
                surface_p.surfaceNumber = reg.getNewSurfaceNumber()
            if not surface_p.surfaceNumber:
                surface_p.surfaceNumber = reg.getNewSurfaceNumber()
            reg.surfaceDict[surface_p.surfaceNumber] = surface_p
            self.addSurface(surface_p)  # also add to the cell's surfaceList
        else:
            self.surfaceList = [surface_p]  # cell's surfaceList

        if reg:
            m1 = pyg4ometry.mcnp.Material(materialNumber=1, density=0.92, reg=reg)  # polyethylene
            self.addMaterial(m1)

        super().__init__(surfaces=[surface_p], geometry=surface_p, cellNumber=cellNumber, materialNumber=2, reg=reg)

        if reg:
            m2 = pyg4ometry.mcnp.Material(materialNumber=2, density=2.699)  # aluminium
            reg.addMaterial(m2)
            self.addMaterial(m2)

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

        # new connector (prime)
        connector_p = Connector(
            cellNumber=self.cellNumber,
        )

        # transform surface
        surfaces_p = [s.transform(rotation=rotationMatrix.tolist(), translation=translationVector.tolist()) for s in self.surfaceList]

        # update the new connector
        connector_p.surfaceList = surfaces_p
        connector_p.geometry = surfaces_p[0]

        # update mesh
        axis, angle = _utils.rotationMatrixToAxisAndAngle(rotationMatrix)
        connector_p._mesh = self._mesh
        connector_p._mesh.rotate(axis, 360-angle)
        connector_p._mesh.translate(translationVector)

        return connector_p

    def mesh(self):
        if self._mesh is not None:
            return self._mesh
        else:
            surface = pyg4ometry.mcnp.RCC(0, 0, 0, 0, 0, length, radius)
            super().__init__(surfaces=[surface], geometry=surface)  # a block is a cell
            return super().mesh()

