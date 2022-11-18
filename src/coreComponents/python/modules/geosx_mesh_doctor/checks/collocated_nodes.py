from collections import defaultdict
from dataclasses import dataclass
import logging
from typing import Tuple

import numpy

from vtkmodules.vtkCommonCore import (
    reference,
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    vtkIncrementalOctreePointLocator, )

from . import vtk_utils


@dataclass(frozen=True)
class Options:
    tolerance: float


@dataclass(frozen=True)
class Result:
    nodes_buckets: Tuple[Tuple[int]]


def __check(mesh, options: Options) -> Result:
    points = mesh.GetPoints()

    locator = vtkIncrementalOctreePointLocator()
    locator.SetTolerance(options.tolerance)
    output = vtkPoints()
    locator.InitPointInsertion(output, points.GetBounds())

    # original ids to/from filtered ids.
    # original_to_filtered = numpy.ones(points.GetNumberOfPoints(), dtype=int) * -1
    filtered_to_original = numpy.ones(points.GetNumberOfPoints(), dtype=int) * -1

    rejected_points = defaultdict(list)
    point_id = reference(0)
    for i in range(points.GetNumberOfPoints()):
        is_inserted = locator.InsertUniquePoint(points.GetPoint(i), point_id)
        if not is_inserted:
            # If it's not inserted, `point_id` contains the node that was already at that location.
            # But in that case, `point_id` is the new numbering in the destination points array.
            # It's more useful for the user to get the old index in the original mesh, so he can look for it in his data.
            logging.debug(
                f"Point {i} at {points.GetPoint(i)} has been rejected, point {filtered_to_original[point_id.get()]} is already inserted."
            )
            rejected_points[point_id.get()].append(i)
        else:
            # If it's inserted, `point_id` contains the new index in the destination array.
            # We store this information to be able to connect the source and destination arrays.
            # original_to_filtered[i] = point_id.get()
            filtered_to_original[point_id.get()] = i

    tmp = []
    for n, ns in rejected_points.items():
        tmp.append((n, *ns))

    return Result(nodes_buckets=tmp)


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)