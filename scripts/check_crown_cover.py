"""Independent scalar/geometry audit; does not replay matrix residuals."""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def point(value):
    require(len(value) == 2, 'point dimension')
    return tuple(map(Q, value))


def check(report):
    require(report['contour_coverage'] is True, 'coverage not claimed')
    require(report['lower_boundary_by_real_conjugation'] is True, 'missing conjugation')
    tiles = report['tiles']
    require(bool(tiles), 'empty tiles')
    defects, inverse_bounds = [], []
    for tile in tiles:
        a, v, r, defect = map(Q, (tile['center_defect'], tile['variation'],
                                  tile['radius'], tile['uniform_defect']))
        require(a >= 0 and v >= 0 and r > 0, 'invalid scalar bound')
        require(defect == a+r*v and defect < 1, 'invalid contraction')
        require(tile['accepted'] is True, 'rejected tile')
        defects.append(defect)
        if 'parametrix_norm' in tile:
            norm = Q(tile['parametrix_norm'])
            require(norm > 0, 'invalid parametrix norm')
            inverse_bounds.append(norm/(1-defect))
    vertices = [(Q(1, 8), Q(0)), (Q(1, 8), Q(1, 4)),
                (Q(1), Q(1, 4)), (Q(1), Q(0))]
    current, leg, used = vertices[0], 0, set()
    for segment in report['upper_boundary_segments']:
        require(leg < 3, 'extra segment')
        start, end = point(segment['start']), point(segment['end'])
        index = segment['tile']
        require(type(index) is int and 0 <= index < len(tiles), 'invalid tile index')
        require(start == current == point(tiles[index]['center']), 'gap or wrong center')
        axis = 0 if leg == 1 else 1
        direction = -1 if leg == 2 else 1
        require(end[1-axis] == start[1-axis], 'off boundary')
        step = direction*(end[axis]-start[axis])
        remaining = direction*(vertices[leg+1][axis]-start[axis])
        require(0 < step <= remaining, 'wrong direction or corner overshoot')
        require(step <= Q(tiles[index]['radius']), 'segment outside tile')
        used.add(index)
        current = end
        if current == vertices[leg+1]:
            leg += 1
    require(leg == 3, 'incomplete boundary')
    require(used == set(range(len(tiles))), 'unused tiles')
    inverse = max(inverse_bounds) if len(inverse_bounds) == len(tiles) else None
    return dict(scalar_geometry_verified=True, tiles=len(tiles),
                max_uniform_defect=str(max(defects)),
                max_uniform_defect_decimal=float(max(defects)),
                uniform_inverse_bound=str(inverse) if inverse is not None else None,
                uniform_inverse_bound_decimal=float(inverse) if inverse is not None else None,
                matrix_residuals_replayed=False, physical_certificate=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    args = parser.parse_args()
    print(json.dumps(check(json.loads(args.report.read_text())), indent=2))
