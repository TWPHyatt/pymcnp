import pyg4ometry
import numpy as _np

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