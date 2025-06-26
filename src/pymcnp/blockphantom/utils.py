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

    def computeRotationMatrix(v1, v2):
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