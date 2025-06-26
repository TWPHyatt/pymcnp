import numpy as _np


def rotationStepsToMatrix(stepsIn):
    """
    Converts a list of 90-degree step rotations [x, y, z] into a 3x3 rotation matrix
    """

    steps = _np.array(stepsIn)

    if len(steps) != 3 or not all(isinstance(i, int) for i in steps):
        raise TypeError("rotation_steps must be a list of 3 integers [x, y, z]")

    def rotX(n):
        theta = n * _np.pi / 2
        return _np.array([
            [1, 0, 0],
            [0, _np.cos(theta), -_np.sin(theta)],
            [0, _np.sin(theta), _np.cos(theta)]
        ])

    def rotY(n):
        theta = n * _np.pi / 2
        return _np.array([
            [_np.cos(theta), 0, _np.sin(theta)],
            [0, 1, 0],
            [-_np.sin(theta), 0, _np.cos(theta)]
        ])

    def rotZ(n):
        theta = n * _np.pi / 2
        return _np.array([
            [_np.cos(theta), -_np.sin(theta), 0],
            [_np.sin(theta), _np.cos(theta), 0],
            [0, 0, 1]
        ])

    rotationMatrix = _np.eye(3)
    for axis, step in zip([rotX, rotY, rotZ], steps):
        if step != 0:
            rotationMatrix = axis(step) @ rotationMatrix

    return rotationMatrix
