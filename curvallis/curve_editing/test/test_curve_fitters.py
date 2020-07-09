import inspect, unittest as ut
from collections import namedtuple

import numpy as np

from curvallis.run import CurveInteractor
from curvallis.curve_editing import curve_fitters as cf


class TestCurveFitters(ut.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def _make_args(self, **kwargs):
        Args = namedtuple('Args', kwargs.keys())
        return Args(*kwargs.values())

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
    #   Poly_Original tests
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

    def _get_default_poly_args(self):
        return {'derivative_scale': 1,
                'second_derivative_scale': 1,
                'integral_scale': 1,
                'x_integral_ref': 0,
                'y_integral_ref': 0}

    def _make_poly_args(self, **kwargs):
        """ use **kwargs to overwrite or add to defaults from _get_default_poly_args """
        poly_args = self._get_default_poly_args()
        poly_args.update(**kwargs)
        # loc = locals()  # capture locals here so python doesn't freak out - TODO: figure out why and make this better
        # argnames = inspect.getfullargspec(self._make_poly_args).args
        # kwargs = {n: loc[n] for n in argnames if n != 'self'}
        return self._make_args(**poly_args)

    def _make_poly_and_points(self, degree, x_min=0, x_max=100, num_pts=0):
        """ Generates a random polynomial of the specified degree and a number of points along that curve

        :param num_pts: the number of points to generate.
                    If <=0, will default to degree * (degree + 1)
        :return: poly, points
                    where poly is a numpy.poly1d object used to generate the points
                    and points is a list of (x,y) tuples with points on the curve represented by poly
        """
        num_pts = int(num_pts)
        if num_pts <= 0:
            num_pts = degree ** 2 + degree
        x = np.linspace(x_min, x_max, num_pts)
        poly = np.poly1d(np.random.randint(10, size=degree + 1))
        return poly, list(zip(x, poly(x)))

    def test_Poly_name_prefix_class(self):
        self.assertEqual(cf.Poly_Original.name_prefix, 'poly')

    def test_Poly_name_prefix_obj(self):
        deg = np.random.default_rng().integers(1, 12, endpoint=True)
        p = cf.Poly_Original(self._make_poly_args(), 'poly{0}'.format(deg))
        self.assertEqual(p.name_prefix, 'poly')

    # TODO: every once in a while this fails just like the noisy version - maybe remove the random nature of the target?
    def test_Poly_fit_to_points(self):
        deg = 5
        actual_f, pts = self._make_poly_and_points(deg)
        poly5 = cf.Poly_Original(self._make_poly_args(), 'poly{0}'.format(deg))
        poly5.fit_to_points(pts)
        np.testing.assert_array_almost_equal(poly5._f.c, actual_f.c, decimal=2)

    @ut.skip('fitter often comes up with wildly different coefficients than the ideal case')
    def test_Poly_fit_to_points_noisy(self):
        deg = 5
        actual_f, pts = self._make_poly_and_points(deg, num_pts=200)
        # add ~2% noise to the y values of the points
        noisy_pts = [(x, y * n) for ((x, y), n) in zip(pts, np.random.normal(1, .02, len(pts)))]
        poly5 = cf.Poly_Original(self._make_poly_args(), 'poly{0}'.format(deg))
        poly5.fit_to_points(noisy_pts)
        np.testing.assert_array_almost_equal(poly5._f.c, actual_f.c, decimal=1)

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
    #   GammaPoly tests
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

    def _get_default_gammapoly_args(self):
        gpoly_args = self._get_default_poly_args()
        gpoly_args['rho0'] = 1
        return gpoly_args

    def _make_gammapoly_args(self, **kwargs):
        gpoly_args = self._get_default_gammapoly_args()
        gpoly_args.update(**kwargs)
        return self._make_args(**gpoly_args)

    def test_GammaPoly_name_prefix_class(self):
        self.assertEqual(cf.GammaPoly.name_prefix, 'gammapoly')

    def test_GammaPoly_name_prefix_obj(self):
        deg = np.random.default_rng().integers(1, 12, endpoint=True)
        p = cf.GammaPoly(self._make_gammapoly_args(), 'gammapoly{0}'.format(deg))
        self.assertEqual(p.name_prefix, 'gammapoly')


if __name__ == '__main__':
    ut.main()
