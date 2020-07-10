import warnings, unittest as ut
from collections import namedtuple

import numpy as np

from curvallis.run import CurveInteractor
from curvallis.curve_editing import curve_fitters as cf


class TestCurveFitters(ut.TestCase):
    def setUp(self):
        warnings.filterwarnings('ignore', message='divide by zero encountered in true_divide', category=RuntimeWarning)

    def tearDown(self):
        warnings.resetwarnings()

    def _get_random_degree(self, min_degree=1, max_degree=12):
        return np.random.default_rng().integers(min_degree, max_degree, endpoint=True)

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
        deg = self._get_random_degree()
        p = cf.Poly_Original(self._make_poly_args(), 'poly{0}'.format(deg))
        self.assertEqual(p.name_prefix, 'poly')

    # TODO: every once in a while this fails just like the noisy version - maybe remove the random nature of the target?
    def test_Poly5_fit_to_points(self):
        deg = 5
        actual_f, pts = self._make_poly_and_points(deg)
        poly5 = cf.Poly_Original(self._make_poly_args(), 'poly{0}'.format(deg))
        poly5.fit_to_points(pts)
        np.testing.assert_array_almost_equal(poly5._f.c, actual_f.c, decimal=2)

    @ut.skip('fitter often comes up with wildly different coefficients than the ideal case')
    def test_Poly5_fit_to_points_noisy(self):
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

    def _make_random_GammaPoly(self, rho_is_density=True, **kwargs):
        deg = self._get_random_degree()
        return cf.GammaPoly(self._make_gammapoly_args(**kwargs), 'gammapoly{0}'.format(deg), rho_is_density)

    def test_GammaPoly_name_prefix_class(self):
        self.assertEqual(cf.GammaPoly.name_prefix, 'gammapoly')

    def test_GammaPoly_name_prefix_obj(self):
        p = self._make_random_GammaPoly()
        self.assertEqual(p.name_prefix, 'gammapoly')

    ## if using density (rho_is_density=True), higher pressure should correspond to higher density

    def test_GammaPoly__get_highP_lowP_x_indices_density_sans_rho0(self):
        p = self._make_random_GammaPoly(rho0=3.5)
        x = np.arange(10)
        hiPi, loPi = p._get_highP_lowP_x_indices(x)
        np.testing.assert_equal(hiPi, np.arange(4, 10), 'High pressure indices are incorrect')
        np.testing.assert_equal(loPi, np.arange(4), 'Low pressure indices are incorrect')

    def test_GammaPoly__get_highP_lowP_x_indices_density_with_rho0(self):
        p = self._make_random_GammaPoly(rho0=3)
        x = np.arange(10)
        hiPi, loPi = p._get_highP_lowP_x_indices(x)
        # confirm rho0 point is in both sets of indices
        np.testing.assert_equal(hiPi, np.arange(3, 10), 'High pressure indices are incorrect')
        np.testing.assert_equal(loPi, np.arange(4), 'Low pressure indices are incorrect')

    def test_GammaPoly__get_highP_lowP_x_indices_density_all_lowP(self):
        p = self._make_random_GammaPoly(rho0=11)
        x = np.arange(10)
        hiPi, loPi = p._get_highP_lowP_x_indices(x)
        self.assertEqual(hiPi.size, 0, 'There should be no high pressure points')
        np.testing.assert_equal(loPi, x, 'All indices should be in the low pressure regime')

    def test_GammaPoly__get_highP_lowP_x_indices_density_all_highP(self):
        p = self._make_random_GammaPoly(rho0=-1)
        x = np.arange(10)
        hiPi, loPi = p._get_highP_lowP_x_indices(x)
        np.testing.assert_equal(hiPi, x, 'All indices should be in the high pressure regime')
        self.assertEqual(loPi.size, 0, 'There should be no low pressure points')

    def test_GammaPoly__get_highP_fit_x_density(self):
        r0 = 3
        p = self._make_random_GammaPoly(rho0=r0)
        x = np.arange(4, 10)
        fit_x = p._get_highP_fit_x(x)
        np.testing.assert_equal(fit_x, np.divide(x, r0), 'High pressure density should do polynomial fit in x/rho0')

    def test_GammaPoly__get_lowP_fit_x_density(self):
        r0 = 5
        p = self._make_random_GammaPoly(rho0=r0)
        x = np.arange(r0)
        fit_x = p._get_lowP_fit_x(x)
        np.testing.assert_equal(fit_x, np.divide(r0, x), 'Low pressure density should do polynomial fit in rho0/x')

    def test_GammaPoly__get_highP_lowP_fit_points_density(self):
        r0 = 5
        p = self._make_random_GammaPoly(rho0=r0)
        data = [(0, 1), (2, 4), (3, 8), (6, 7), (8, 6)]
        hiP_points, loP_points = p._get_highP_lowP_fit_points(data)
        np.testing.assert_equal(hiP_points, np.array([(6/r0, 7), (8/r0, 6)]), 'High pressure points are incorrect')
        np.testing.assert_equal(loP_points, np.array([(np.inf, 1), (r0/2, 4), (r0/3, 8)]), 'Low pressure points are incorrect')


    ## if using volume (rho_is_density=False, higher pressure should correspond to lower unit volume

    def test_GammaPoly__get_highP_lowP_x_indices_volume_sans_rho0(self):
        # if using volume, high pressure should correspond to low unit volume
        p = self._make_random_GammaPoly(rho_is_density=False, rho0=3.5)
        x = np.arange(10)
        hiPi, loPi = p._get_highP_lowP_x_indices(x)
        np.testing.assert_equal(hiPi, np.arange(4), 'High pressure indices are incorrect')
        np.testing.assert_equal(loPi, np.arange(4, 10), 'Low pressure indices are incorrect')

    def test_GammaPoly__get_highP_lowP_x_indices_volume_with_rho0(self):
        p = self._make_random_GammaPoly(rho_is_density=False, rho0=3)
        x = np.arange(10)
        hiPi, loPi = p._get_highP_lowP_x_indices(x)
        # confirm rho0 point is in both sets of indices
        np.testing.assert_equal(hiPi, np.arange(4), 'High pressure indices are incorrect')
        np.testing.assert_equal(loPi, np.arange(3, 10), 'Low pressure indices are incorrect')

    def test_GammaPoly__get_highP_lowP_x_indices_volume_all_lowP(self):
        p = self._make_random_GammaPoly(rho_is_density=False, rho0=-1)
        x = np.arange(10)
        hiPi, loPi = p._get_highP_lowP_x_indices(x)
        self.assertEqual(hiPi.size, 0, 'There should be no high pressure points')
        np.testing.assert_equal(loPi, x, 'All indices should be in the low pressure regime')

    def test_GammaPoly__get_highP_lowP_x_indices_volume_all_highP(self):
        p = self._make_random_GammaPoly(rho_is_density=False, rho0=11)
        x = np.arange(10)
        hiPi, loPi = p._get_highP_lowP_x_indices(x)
        np.testing.assert_equal(hiPi, x, 'All indices should be in the high pressure regime')
        self.assertEqual(loPi.size, 0, 'There should be no low pressure points')

    def test_GammaPoly__get_highP_fit_x_volume(self):
        r0 = 5
        p = self._make_random_GammaPoly(rho_is_density=False, rho0=r0)
        x = np.arange(0, 5)
        fit_x = p._get_highP_fit_x(x)
        np.testing.assert_equal(fit_x, np.divide(r0, x), 'High pressure volume should do polynomial fit in V0/x')

    def test_GammaPoly__get_highP_lowP_fit_points_volume(self):
        r0 = 5
        p = self._make_random_GammaPoly(rho_is_density=False, rho0=r0)
        data = [(0, 1), (2, 4), (3, 8), (6, 7), (8, 6)]
        hiP_points, loP_points = p._get_highP_lowP_fit_points(data)
        np.testing.assert_equal(hiP_points, np.array([(np.inf, 1), (r0/2, 4), (r0/3, 8)]), 'High pressure points are incorrect')
        np.testing.assert_equal(loP_points, np.array([(6/r0, 7), (8/r0, 6)]), 'Low pressure points are incorrect')








if __name__ == '__main__':
    ut.main()
