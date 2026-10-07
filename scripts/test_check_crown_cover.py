import copy
import unittest
from check_crown_cover import check


def fixture():
    points = [['1/8', '0'], ['1/8', '1/4'], ['1', '1/4'], ['1', '0']]
    return dict(contour_coverage=True, lower_boundary_by_real_conjugation=True,
                tiles=[dict(center=p, center_defect='1/8', variation='1/8',
                            radius='1', uniform_defect='1/4', accepted=True,
                            parametrix_norm='3') for p in points[:-1]],
                upper_boundary_segments=[dict(start=points[i], end=points[i+1], tile=i)
                                         for i in range(3)])


class CoverTests(unittest.TestCase):
    def test_valid_and_optional_inverse(self):
        data = fixture()
        self.assertEqual(check(data)['uniform_inverse_bound'], '4')
        del data['tiles'][0]['parametrix_norm']
        self.assertIsNone(check(data)['uniform_inverse_bound'])

    def test_tampering_rejected(self):
        for field, value in [('uniform_defect', '1/8'), ('radius', '1/16'),
                             ('variation', '-1/8'), ('accepted', False)]:
            data = fixture()
            data['tiles'][0][field] = value
            with self.assertRaises(ValueError):
                check(data)
        data = fixture()
        for segments in (data['upper_boundary_segments'][:-1],
                         list(reversed(data['upper_boundary_segments']))):
            bad = copy.deepcopy(data)
            bad['upper_boundary_segments'] = segments
            with self.assertRaises(ValueError):
                check(bad)

    def test_geometry_and_conjugation_rejected(self):
        for end in (['1/8', '1/3'], ['1/7', '1/4'], ['1/8', '-1/4']):
            data = fixture()
            data['upper_boundary_segments'][0]['end'] = end
            with self.assertRaises(ValueError):
                check(data)
        data = fixture()
        data['lower_boundary_by_real_conjugation'] = False
        with self.assertRaises(ValueError):
            check(data)
        data = fixture()
        data['upper_boundary_segments'][0]['tile'] = -1
        with self.assertRaises(ValueError):
            check(data)


if __name__ == '__main__':
    unittest.main()
