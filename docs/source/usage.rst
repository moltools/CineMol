Usage
=====
The easiest way to use CineMol is to use the functions from the `cinemol.api` module.

This module provides a high-level API for drawing SVGs of molecular structures.

Using the Python API
---------------------

.. code-block:: python

    from cinemol.parsers import parse_sdf
    from cinemol.api import Style, Look, draw_molecule

    # Parse atoms and bonds from an SDF file ...
    with open("path/to/your/molecule.sdf", "r") as f:
        sdf_str = f.read()

    atoms, bonds = parse_sdf(sdf_str)

    # ... or create your own atom and bond objects.
    atoms = [
        Atom(
            index=0,
            symbol="C",
            coordinates=(0.0, 0.0, 0.0)
        ),
        Atom(
            index=1,
            symbol="N",
            coordinates=(1.0, 0.0, 0.0)
        )
    ]
    bonds = [
        Bond(
            start_index=0,
            end_index=1,
            order=3
        )
    ]

    # Draw molecule to SVG string.
    svg = draw_molecule(
        atoms=atoms,
        bonds=bonds,
        style=Style.Tube,  # See cinemol.api.Style for more options.
        look=Look.Glossy,  # See cinemol.api.Look for more options.
        resolution=100,
        scale=10.0
    )

    svg_str = svg.to_svg()

Fine-grained control
---------------------
See :class:`Scene` in `cinemol.model` module if you want direct access to
the drawing scene. The scene can be given various objects for drawing (e.g.,
:class:`ModelSphere`, :class:`ModelCylinder`, :class:`ModelWire`). The :func:`draw_molecule`
function is a convenience function on top of :class:`Scene` that creates a scene,
adds the molecule to the scene, and draws the scene to an SVG string.

Examples:

#. Highlighted amino acids in daptomycin conformer (`examples/draw_substructure_highlights.py <https://github.com/moltools/CineMol/blob/main/examples/draw_substructure_highlights.py>`_)

#. Three aligned conformers of benzylphenol (`examples/draw_superimposed_conformers.py <https://github.com/moltools/CineMol/blob/main/examples/draw_superimposed_conformers.py>`_)

#. Wireframe model of lysozome `9LYZ <https://www.rcsb.org/structure/9lyz>`_ with space-filling model of bound ligand trisaccharide NAM-NAG-NAM (`examples/draw_protein_with_ligands.py <https://github.com/moltools/CineMol/blob/main/examples/draw_protein_with_ligands.py>`_)
