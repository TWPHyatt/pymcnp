import pyg4ometry
import numpy as _np
from ..blockphantom import utils as _utils

class Connector(pyg4ometry.mcnp.Cell):
    def __init__(self, translation=[0, 0, 0], rotationSteps=[0, 0, 0], length=1.5, cellNumber=None, reg=None):

        rotationMatrix = _utils.rotationStepsToMatrix(rotationSteps)
        translationVector = _np.array(translation)

        surface = pyg4ometry.mcnp.RCC(0, 0, 0, 0, 0, length, 0.3)

        surface_p = surface.transform(translation=translationVector.tolist(), rotation=rotationMatrix.tolist())

        super().__init__(surfaces=[surface_p], geometry=surface_p, cellNumber=cellNumber, reg=reg)

        if reg:
            m2 = pyg4ometry.mcnp.Material(materialNumber=2, density=2.699)  # aluminium
            reg.addMaterial(m2)

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
            reg=reg
        )

        # transform surface
        surfaces_p = [s.transform(rotation=rotationMatrix.tolist(), translation=translationVector.tolist()) for s in self.surfaceList]

        # update the new connector
        connector_p.surfaceList = surfaces_p
        connector_p.geometry = surfaces_p[0]

        return connector_p
