import numpy as _np


def rotationStepsToMatrix(stepsIn):
    """
    converts a list of 90-degree step rotations [x, y, z] into a 3x3 rotation matrix
    """

    #if len(stepsIn) != 3 or not all(isinstance(i, int) for i in stepsIn):
    #    msg = f"rotation steps must be a list of 3 integers [x, y, z]"
    #    raise TypeError(msg)

    steps = _np.array(stepsIn)

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
    """
    computes the rotation matrix that aligns vector v1 to vector v2 using Rodrigues' rotation formula
    """
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

def rotationMatrixToAxisAndAngle(R):
    """
    Converts a rotation matrix to axis-angle representation.
    Returns (axis, angle in degrees)
    """
    angle = _np.arccos((_np.trace(R) - 1) / 2)
    if _np.isclose(angle, 0):
        return _np.array([1, 0, 0]), 0.0  # No rotation

    rx = R[2, 1] - R[1, 2]
    ry = R[0, 2] - R[2, 0]
    rz = R[1, 0] - R[0, 1]
    axis = _np.array([rx, ry, rz])
    axis = axis / _np.linalg.norm(axis)

    return axis, _np.degrees(angle)

def rotationAroundAxis(axis, angleRad):
    """
    Rodrigues’ rotation formula for rotation about an arbitrary axis.
    """
    axis = axis / _np.linalg.norm(axis)
    K = _np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    return _np.eye(3) + _np.sin(angleRad) * K + (1 - _np.cos(angleRad)) * (K @ K)

