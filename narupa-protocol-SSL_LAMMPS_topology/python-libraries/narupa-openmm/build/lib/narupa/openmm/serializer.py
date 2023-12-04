# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Serialize and deserialize OpenMM simulations to and from XML files.

A simulation is described as the concatenation of a starting structure as a PDB
file, an OpenMM serialized system, an OpenMM serialized integrator, and,
optionally, an OpenMM serialized state. The resulting XML file looks like:

::

    <OpenMMSimulation>
        <pdbx>
            // pasted content of the PDBx file
        </pdbx>
        <System ...>
            // XML content of the OpenMM serialized system
        </System>
        <Integrator ...>
            // XML content of the OpenMM serialized integrator
        </Integrator>
    </OpenMMSimulation>

The ``System`` and ``Integrator`` tags are the roots of the serialized system
and integrator, respectively. The ``pdbx`` tag can be replaced by a ``pdb``
one for backward compatibility.

This module provides a function :func:`serialize_simulation` that generates an
XML file from an existing instance of :class:`openmm.app.Simulation`, and
a function :func:`deserialize_simulation` that creates an instance of simulation
from an XML file.
"""
from typing import Optional, List, Tuple
from io import StringIO
from tempfile import TemporaryDirectory
import os
from xml.dom.minidom import getDOMImplementation, parseString, Document, Element

from openmm import app, XmlSerializer, CustomExternalForce

from .imd import populate_imd_force

ROOT_TAG = 'OpenMMSimulation'


def serialize_simulation(simulation: app.Simulation) -> str:
    """
    Generate an XML string from a simulation.

    :param simulation: The simulation to serialize.
    :return: A string with the content of an XML file describing the simulation.
    """
    implementation = getDOMImplementation()
    if implementation is None:
        raise TypeError
    document = implementation.createDocument(None, ROOT_TAG, None)

    # Extract the PDB
    positions = simulation.context.getState(getPositions=True).getPositions()
    pdb_content = StringIO()
    app.PDBxFile.writeFile(simulation.topology, positions, pdb_content)
    pdb_node = document.createElement('pdbx')
    pdb_node.appendChild(document.createTextNode(pdb_content.getvalue()))

    # Extract the system
    system_xml_str = XmlSerializer.serialize(simulation.system)
    system_document = parseString(system_xml_str)

    # Extract the integrator
    integrator_xml_str = XmlSerializer.serialize(simulation.integrator)
    integrator_document = parseString(integrator_xml_str)

    # Combine the element in a single
    root = document.documentElement
    root.appendChild(pdb_node)
    root.appendChild(system_document.documentElement)
    root.appendChild(integrator_document.documentElement)

    return root.toprettyxml()


def deserialize_simulation(
        xml_content: str,
        imd_force: Optional[CustomExternalForce] = None,
) -> app.Simulation:
    """
    Create an OpenMM simulation from XML.

    :param xml_content: The content of an XML file as a string.
    :param imd_force: Optionally, an imd force to populate and add to the
        system. The force must be created by
        :func:`narupa.openmm.potentials.create_imd_force`.
    :return: An instance of the simulation.
    """
    document = parseString(xml_content)

    tag, pdb_node = _get_one_exclusive(document, ['pdbx', 'pdb'])
    node = pdb_node.firstChild
    if node is None:
        raise IOError('No structure content.')
    # Mypy fails on the next line with the following error:
    # "Node" has no attribute "nodeValue"  [attr-defined]
    # Because we know the attribute actually exists, we ignore the error.
    pdb_content = StringIO(node.nodeValue)  # type: ignore[attr-defined]
    with TemporaryDirectory() as tmp_dir:
        pdb_path = os.path.join(tmp_dir, f'configuration.{tag}')
        with open(str(pdb_path), 'w') as outfile:
            outfile.write(pdb_content.getvalue())
        if tag == 'pdb':
            pdb = app.PDBFile(str(pdb_path))
        elif tag == 'pdbx':
            pdb = app.PDBxFile(str(pdb_path))
        else:
            raise IOError('Invalid structure tag: {tag}')

    system_node = _get_node_and_raise_if_more_than_one(document, 'System')
    system_content = system_node.toprettyxml()
    system = XmlSerializer.deserialize(system_content)

    if imd_force is not None:
        populate_imd_force(imd_force, system)
        system.addForce(imd_force)

    integrator_node = _get_node_and_raise_if_more_than_one(document, 'Integrator')
    integrator_content = integrator_node.toprettyxml()
    integrator = XmlSerializer.deserialize(integrator_content)

    simulation = app.Simulation(
        topology=pdb.topology,
        system=system,
        integrator=integrator,
    )
    simulation.context.setPositions(pdb.positions)
    return simulation


def _get_node_and_raise_if_more_than_one(document: Document, tag_name: str) -> Element:
    nodes = document.getElementsByTagName(tag_name)
    if not nodes:
        raise IOError('No {} tag defined in the XML.'.format(tag_name))
    if len(nodes) != 1:
        raise IOError('More than one {} tag defined in the XML.'.format(tag_name))
    return nodes[0]


def _get_one_exclusive(document: Document, tag_names: List[str]) -> Tuple[str, Element]:
    content: List[Tuple[str, Element]] = []
    for name in tag_names:
        content.extend((name, node) for node in document.getElementsByTagName(name))
    if not content:
        raise IOError(f'No data for any of these tags: {tag_names}')
    if len(content) > 1:
        raise IOError(f'More than one of these tags defined in the XML: {tag_names}')
    return content[0]
