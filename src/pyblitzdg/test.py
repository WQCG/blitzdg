import pyblitzdg
import pprint

print(dir(pyblitzdg))
nodes1d = pyblitzdg.Nodes1DProvisioner(2, 50, -1.0, -1.0)

print(nodes1d.numLocalPoints)
