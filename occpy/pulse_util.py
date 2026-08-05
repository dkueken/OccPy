"""
Pulse- and return-related helpers for OccPy
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Tuple

import laspy
import numpy as np
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

class ReturnMode(Enum):
    """What the return-number fields of a dataset actually contain."""

    #: At least one pulse has more than one echo.
    MULTI = "multi"

    #: Every echo is return 1 of 1 (max return_number == 1, max number_of_returns == 1).
    SINGLE = "single"

    #: return_number and number_of_returns are 0 everywhere. Seen with some
    #: sensors (e.g. GeoSLAM ZebHorizon). The fields carry no information, so the data can
    #: only be treated as single-return.
    DEGENERATE_ZERO = "degenerate_zero"

    @property
    def implies_single_return(self) -> bool:
        return self is not ReturnMode.MULTI


class TraceMode(Enum):
    """How pulses are handed to the ray tracer."""

    #: Every echo is its own pulse; traced chunk-by-chunk, pulse dataset unused.
    SINGLE_RETURN = "single_return"

    #: Multi-return, GPS-time sorted: trace each chunk immediately, then clear.
    STREAMING = "streaming"

    #: Multi-return, not sorted: accumulate everything, trace once at the end.
    DEFERRED = "deferred"


@dataclass(frozen=True)
class ReturnModeResult:
    """Outcome of return-mode detection, with provenance for logging."""

    mode: ReturnMode
    source: str          # "header" or "scan"
    detail: str          # human-readable reason, for the log
    points_inspected: int = 0


@dataclass
class ScanJob:
    """One LAZ file plus the sensor position(s) that belong to it."""

    laz_file: str
    name: str
    #: Fixed sensor position for static acquisitions, or None for mobile
    #: acquisitions, where positions are interpolated per GPS time.
    sensor_xyz: Optional[Tuple[float, float, float]]

    @property
    def is_mobile(self) -> bool:
        return self.sensor_xyz is None


@dataclass
class Chunk:
    """A block of echoes read from a LAZ file, normalised for the C++ side."""

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    gps_time: np.ndarray
    return_number: np.ndarray
    number_of_returns: np.ndarray

    def __len__(self) -> int:
        return int(self.x.size)


# ---------------------------------------------------------------------------
# Return-mode detection
# ---------------------------------------------------------------------------

def return_mode_from_header(header) -> Optional[ReturnModeResult]:
    """
    Determine the return mode from the LAS header alone -- no point data read.

    LAS stores ``number_of_points_by_return``. Many writers get it wrong, though: some leave it zeroed,
    and it is often not updated after a file has been filtered. It is therefore
    only trusted when the per-return counts add up to the header point count.

    Returns None when the header cannot be trusted, meaning the caller should
    fall back to scanning the points.

    Parameters
    ----------
    header: laspy.header

    Returns
    -------

    """

    counts = np.asarray(getattr(header, "number_of_points_by_return", []), dtype=np.int64)
    total = int(getattr(header, "point_count", 0))

    if counts.size == 0:
        return None
    if total <= 0:
        return None

    counted = int(counts.sum())
    if counted == 0:
        # Field never populated by the writer.
        return None
    if counted != total:
        # Stale after filtering, or the writer disagrees with itself. Either way
        # the counts describe a file that is not the one in front of us.
        return None

    if int(counts[1:].sum()) > 0:
        return ReturnModeResult(
            mode=ReturnMode.MULTI,
            source="header",
            detail=f"number_of_points_by_return reports {int(counts[1:].sum())} echoes "
                   f"beyond return 1 (max return index {int(np.max(np.nonzero(counts)[0])) + 1})",
            points_inspected=total,
        )

    return ReturnModeResult(
        mode=ReturnMode.SINGLE,
        source="header",
        detail=f"number_of_points_by_return places all {total} echoes in return 1",
        points_inspected=total,
    )

def scan_return_mode(path: str, points_per_iter: int, progress: bool = True) -> ReturnModeResult:
    """
    Determine the return mode by reading the point data.

    Only ``return_number`` and ``number_of_returns`` are retained, so memory
    stays flat regardless of file size. The scan short-circuits as soon as a
    return number above 1 appears, which for genuine multi-return data is
    almost always inside the first chunk. Only a true single-return file costs
    a full pass.

    Parameters
    ----------
    path: str
        path to input laz
    points_per_iter: int
        points to be used per iteration
    progress: bool
        whether progressbar should be shown

    Returns
    -------
    ReturnModeResult

    """
    max_return_number = 0
    max_number_of_returns = 0
    inspected = 0

    with laspy.open(path) as file:
        total = file.header.point_count
        with tqdm(
            total=total,
            desc=f"Checking returns in {os.path.basename(path)}",
            unit="pts",
            disable=not progress,
            leave=False,
        ) as pbar:
            for points in file.chunk_iterator(points_per_iter):
                n = len(points)
                if n == 0:
                    continue

                return_number = np.asarray(points.return_number)
                number_of_returns = np.asarray(points.number_of_returns)

                max_return_number = max(max_return_number, int(return_number.max()))
                max_number_of_returns = max(max_number_of_returns, int(number_of_returns.max()))
                inspected += n
                pbar.update(n)

                if max_return_number > 1 or max_number_of_returns > 1:
                    return ReturnModeResult(
                        mode=ReturnMode.MULTI,
                        source="scan",
                        detail=f"found number_of_returns up to {max_number_of_returns} "
                               f"after inspecting {inspected} of {total} echoes",
                        points_inspected=inspected,
                    )

    if inspected == 0:
        raise ValueError(f"{path} contains no points, cannot determine the return mode.")

    if max_return_number == 0 and max_number_of_returns == 0:
        return ReturnModeResult(
            mode=ReturnMode.DEGENERATE_ZERO,
            source="scan",
            detail=f"return_number and number_of_returns are 0 for all {inspected} echoes",
            points_inspected=inspected,
        )

    return ReturnModeResult(
        mode=ReturnMode.SINGLE,
        source="scan",
        detail=f"all {inspected} echoes are return 1 of 1",
        points_inspected=inspected,
    )


def detect_return_mode(path: str, points_per_iter: int, progress: bool = True) -> ReturnModeResult:
    """
    detect whether input file is single return or multi-return data. We check header first
    and fall back to the actual point cloud if header cannot be trused
    Parameters
    ----------
    path: str
        path to input laz
    points_per_iter: int
        points to be used per iteration
    progress: bool
        whether progressbar should be shown

    Returns
    -------
    ReturnModeResult

    """

    """Header first, point scan only if the header cannot be trusted."""
    with laspy.open(path) as file:
        from_header = return_mode_from_header(file.header)

    if from_header is not None:
        return from_header

    return scan_return_mode(path, points_per_iter, progress=progress)

# ---------------------------------------------------------------------------
# Chunk handling
# ---------------------------------------------------------------------------

def normalise_chunk(points, source_name: str, require_gps_time: bool) -> Chunk:
    """
    Copy laspy views into contiguous arrays.

    The C++ bindings need contiguous memory, and laspy's sub-field views are neither contiguous nor writable
    Parameters
    ----------
    points: LaspyPoints
    source_name: str
    require_gps_time: bool

    Returns
    -------

    """

    x = np.ascontiguousarray(points.x, dtype=np.float64)
    y = np.ascontiguousarray(points.y, dtype=np.float64)
    z = np.ascontiguousarray(points.z, dtype=np.float64)

    if "gps_time" in set(points.point_format.dimension_names):
        gps_time = np.ascontiguousarray(points.gps_time, dtype=np.float64)
    elif require_gps_time:
        raise ValueError(
            f"{source_name} has no gps_time dimension, which is required to group echoes "
            f"into pulses for multi-return data."
        )
    else:
        gps_time = np.zeros(x.size, dtype=np.float64)

    return Chunk(
        x=x,
        y=y,
        z=z,
        gps_time=gps_time,
        return_number=np.ascontiguousarray(np.asarray(points.return_number)),
        number_of_returns=np.ascontiguousarray(np.asarray(points.number_of_returns)),
    )


def chunk_sortedness(chunk: Chunk, prev_gps_max: Optional[float]) -> Tuple[bool, str]:
    """
    Check GPS-time monotonicity within the chunk and across the chunk boundary.

    Returns (is_sorted, reason).
    Parameters
    ----------
    chunk: Chunk
    prev_gps_max: float

    Returns
    -------
    Tuple[bool, str]

    """
    if len(chunk) == 0:
        return True, ""

    if not bool(np.all(chunk.gps_time[:-1] <= chunk.gps_time[1:])):
        return False, "gps_time decreases within a chunk"

    if prev_gps_max is not None and prev_gps_max > float(chunk.gps_time[0]):
        return False, "gps_time decreases across a chunk boundary"

    return True, ""


def constant_sensor_positions(
    n: int, sensor_xyz: Tuple[float, float, float]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-echo sensor position arrays for a static (TLS) acquisition."""
    sx, sy, sz = sensor_xyz
    return (
        np.full(n, sx, dtype=np.float64),
        np.full(n, sy, dtype=np.float64),
        np.full(n, sz, dtype=np.float64),
    )


def single_return_pulse_ids(chunk: Chunk, offset: int) -> np.ndarray:
    """
    Pulse identifiers for single-return tracing.

    Each echo is its own pulse, so the identifier only has to be unique within
    the batch being traced. Real GPS times are used when available, otherwise a
    running counter.
    """
    if np.any(chunk.gps_time != 0.0):
        return chunk.gps_time
    return np.arange(offset + 1, offset + len(chunk) + 1, dtype=np.float64)
