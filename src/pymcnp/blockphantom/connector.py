import pyg4ometry


class Connector(pyg4ometry.mcnp.Cell):
    def __init__(self, translation=[0, 0, 0], rotation=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], length=1.5, cellNumber=None, reg=None):

        surface = pyg4ometry.mcnp.RCC(0, 0, 0, 0, 0, length, 0.3)

        surface_p = surface.transform(rotation=rotation, translation=translation)

        super().__init__(surfaces=surface_p, geometry=surface_p, cellNumber=cellNumber, reg=reg)

        if reg:
            m2 = pyg4ometry.mcnp.Material(materialNumber=2, density=2.699)  # aluminium
            reg.addMaterial(m2)
