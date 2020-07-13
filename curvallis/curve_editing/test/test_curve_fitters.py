import re, warnings
import unittest as ut
from collections import namedtuple

import numpy as np

from curvallis.curve_editing import curve_fitters as cf


class TestCurveFitters(ut.TestCase):
    def setUp(self):
        warnings.filterwarnings('ignore', message='divide by zero encountered in true_divide', category=RuntimeWarning)

    def tearDown(self):
        warnings.resetwarnings()

    def _make_args(self, **kwargs):
        Args = namedtuple('Args', kwargs.keys())
        return Args(*kwargs.values())

    def _replace_if_zero(self, val):
        return val if val != 0 else np.spacing(1)

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
    #   Factory tests
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

    def test_Factory_get_sorted_fit_names_Poly_Original(self):
        fit_names = cf.factory.get_sorted_fit_names()
        # all expected degrees are present
        for i in range(1, 13):
            self.assertTrue(cf.Poly_Original.name_prefix+str(i) in fit_names)
        # no unexpected degree is present
        self.assertFalse(cf.Poly_Original.name_prefix in fit_names)
        r = re.compile(cf.Poly_Original.name_prefix+'(1|2|3|4|5|6|7|8|9|10|11|12)')
        self.assertEqual(len(list(filter(r.match, fit_names))), 12)

    def test_Factory_get_sorted_fit_names_GammaPoly(self):
        fit_names = cf.factory.get_sorted_fit_names()
        for i in range(1, 13):
            self.assertTrue(cf.GammaPoly.name_prefix + str(i) in fit_names)
        # no unexpected degree is present
        self.assertFalse(cf.GammaPoly.name_prefix in fit_names)
        r = re.compile(cf.GammaPoly.name_prefix + '(1|2|3|4|5|6|7|8|9|10|11|12)')
        self.assertEqual(len(list(filter(r.match, fit_names))), 12)

    def test_Factory_get_sorted_refine_fit_names_Poly_Original(self):
        refine_fit_names = cf.factory.get_sorted_refine_fit_names()
        for i in range(1, 13):
            self.assertTrue(cf.Poly_Original.name_prefix + str(i) in refine_fit_names)
        # no unexpected degree is present
        self.assertFalse(cf.Poly_Original.name_prefix in refine_fit_names)
        r = re.compile(cf.Poly_Original.name_prefix + '(1|2|3|4|5|6|7|8|9|10|11|12)')
        self.assertEqual(len(list(filter(r.match, refine_fit_names))), 12)

    def test_Factory_make_object_of_class_Poly_Original(self):
        degree = self._get_random_degree()
        poly = cf.factory.make_object_of_class(cf.Poly_Original.name_prefix+str(degree), self._make_poly_args())
        self.assertEqual(type(poly), cf.Poly_Original)
        self.assertEqual(poly._order, degree)

    def test_Factory_make_object_of_class_GammaPoly(self):
        degree = self._get_random_degree()
        gpoly = cf.factory.make_object_of_class(cf.GammaPoly.name_prefix+str(degree), self._make_poly_args(rho0=5))
        self.assertEqual(type(gpoly), cf.GammaPoly)
        self.assertEqual(gpoly._order, degree, 'GammaPoly object has wrong degree')
        self.assertEqual(gpoly._hiP_fitter._order, degree, 'GammaPoly high pressure fitter has wrong degree')
        self.assertEqual(gpoly._loP_fitter._order, degree, 'GammaPoly low pressure fitter has wrong degree')



    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
    #   Poly_Original tests
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

    def _get_random_degree(self, min_degree=1, max_degree=12):
        return np.random.default_rng().integers(min_degree, max_degree, endpoint=True)

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
        x = np.logspace(np.log10(self._replace_if_zero(x_min)), np.log10(self._replace_if_zero(x_max)), num_pts)
        return self._make_poly_and_setpoints(degree, x)

    def _make_random_poly1d(self, degree):
        return np.poly1d(np.random.randint(10, size=degree + 1).astype(np.float))

    def _make_poly_and_setpoints(self, degree, x):
        poly = self._make_random_poly1d(degree)
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

    @ut.skip('fitter usually comes up with wildly different coefficients than the ideal case')
    def test_Poly5_fit_to_points_noisy(self):
        deg = 5
        actual_f, pts = self._make_poly_and_points(deg, num_pts=200)
        # add ~2% noise to the y values of the points
        noisy_pts = [(x, y * n) for ((x, y), n) in zip(pts, np.random.normal(1, .01, len(pts)))]
        poly5 = cf.Poly_Original(self._make_poly_args(), 'poly{0}'.format(deg))
        poly5.fit_to_points(noisy_pts)
        np.testing.assert_array_almost_equal(poly5._f.c, actual_f.c, decimal=0)

    def _test_Poly__set_poly(self):
        deg = self._get_random_degree()
        f = self._make_random_poly1d(deg)
        df = np.polyder(f)
        poly = cf.Poly_Original(self._make_poly_args(), 'poly{0}'.format(deg))
        poly._set_poly(f)
        np.testing.assert_array_equal(poly._f.c, f.c)
        np.testing.assert_array_equal(poly._der.c, df.c)
        np.testing.assert_array_equal(poly._scnd_der.c, np.polyder(df).c)
        np.testing.assert_array_equal(poly._int.c, np.polyint(f).c)



    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
    #   GammaPoly tests
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

    def _get_default_gammapoly_args(self):
        gpoly_args = self._get_default_poly_args()
        gpoly_args['rho0'] = 1
        return gpoly_args

    def _make_gammapoly_args(self, **kwargs):
        gpoly_args = self._get_default_gammapoly_args()
        gpoly_args.update(**kwargs)     # will overwrite any defaults
        return self._make_args(**gpoly_args)

    def _make_GammaPoly(self, degree, rho_is_density=True, **kwargs):
        return cf.GammaPoly(self._make_gammapoly_args(**kwargs), 'gammapoly{0}'.format(degree), rho_is_density)

    def _make_random_GammaPoly(self, rho_is_density=True, **kwargs):
        return self._make_GammaPoly(self._get_random_degree(), rho_is_density, **kwargs)

    def _make_gammapoly_and_points(self, degree, rho0, x_min=0.1, x_max=100, num_pts=0, rho_is_density=True):
        """ Generates two 1-D polynomials of the specified degree (one for high pressure, the other for low)
            as well as a number of points to fit to those curves
            most parameters are analogous to those in _make_poly_and_points
        :param rho0: the reference density if rho_is_density, or unit volume otherwise -- this serves as the
                    break point between high and low pressure regimes
        :param num_pts: the number of points generated for each of the high pressure and low pressure regimes
        :param rho_is_density: controls whether independent variable is in density (True) or unit volume (False)
        :return: gammapoly, hiP_poly, hiP_points, loP_poly, loP_points
                    gammapoly is the fitter object
                    the [hi|lo]P_poly variables are numpy.poly1d objects used to generate the points.  Note that the
                        independent variable in these is actually some ratio of x and rho0.  See GammaPoly.
                    the [hi|lo]P_points variables are lists of (x,y) tuples, as expected by GammaPoly.fit_to_points
        """
        def get_regime_x(xmin, xmax, numpts):
            return np.logspace(np.log10(self._replace_if_zero(xmin)), np.log10(self._replace_if_zero(xmax)), numpts)

        num_pts = int(num_pts)
        if num_pts <= 0:
            num_pts = degree ** 2 + degree
        # this object will be discarded -- creating it just to use existing well-tested logic re: hi/lo pressure regimes
        gammapoly = self._make_GammaPoly(degree, rho_is_density, rho0=rho0)
        x = np.concatenate((get_regime_x(x_min, rho0, num_pts), get_regime_x(rho0, x_max, num_pts)))
        hiP_xi, loP_xi = gammapoly._get_highP_lowP_x_indices(x)
        hiP_poly, hiP_fitxy = self._make_poly_and_setpoints(degree, gammapoly._get_highP_fit_x(x[hiP_xi]))
        loP_poly, loP_fitxy = self._make_poly_and_setpoints(degree, gammapoly._get_lowP_fit_x(x[loP_xi]))
        hiP_y = [y for (fitx, y) in hiP_fitxy]
        loP_y = [y for (fitx, y) in loP_fitxy]
        return gammapoly, hiP_poly, list(zip(x[hiP_xi], hiP_y)), loP_poly, list(zip(x[loP_xi], loP_y))

    # TODO: figure out a nice way to combine these next 4 functions
    def _test_gammapoly_func(self, degree, rho0, x, rho_is_density=True):
        gammapoly = self._make_GammaPoly(degree, rho_is_density=rho_is_density, rho0=rho0)
        hiP_poly, _ = self._make_poly_and_setpoints(degree, [])
        loP_poly, _ = self._make_poly_and_setpoints(degree, [])
        gammapoly._hiP_fitter._set_poly(hiP_poly)  # artificially set these as the fits
        gammapoly._loP_fitter._set_poly(loP_poly)
        hiP_xi, loP_xi = gammapoly._get_highP_lowP_x_indices(x)
        y = np.empty(x.shape, np.float)
        y[hiP_xi] = hiP_poly(gammapoly._get_highP_fit_x(x[hiP_xi]))
        y[loP_xi] = loP_poly(gammapoly._get_lowP_fit_x(x[loP_xi]))
        np.testing.assert_array_equal(gammapoly.func(x), y)

    def _test_gammapoly_derivative(self, degree, rho0, x, rho_is_density=True):
        gammapoly = self._make_GammaPoly(degree, rho_is_density=rho_is_density, rho0=rho0)
        hiP_poly, _ = self._make_poly_and_setpoints(degree, [])
        loP_poly, _ = self._make_poly_and_setpoints(degree, [])
        gammapoly._hiP_fitter._set_poly(hiP_poly)  # artificially set these as the fits
        gammapoly._loP_fitter._set_poly(loP_poly)
        hiP_xi, loP_xi = gammapoly._get_highP_lowP_x_indices(x)
        y = np.empty(x.shape, np.float)
        y[hiP_xi] = np.polyder(hiP_poly)(gammapoly._get_highP_fit_x(x[hiP_xi]))
        y[loP_xi] = np.polyder(loP_poly)(gammapoly._get_lowP_fit_x(x[loP_xi]))
        np.testing.assert_array_equal(gammapoly.derivative(x), y)

    def _test_gammapoly_second_derivative(self, degree, rho0, x, rho_is_density=True):
        gammapoly = self._make_GammaPoly(degree, rho_is_density=rho_is_density, rho0=rho0)
        hiP_poly, _ = self._make_poly_and_setpoints(degree, [])
        loP_poly, _ = self._make_poly_and_setpoints(degree, [])
        gammapoly._hiP_fitter._set_poly(hiP_poly)  # artificially set these as the fits
        gammapoly._loP_fitter._set_poly(loP_poly)
        hiP_xi, loP_xi = gammapoly._get_highP_lowP_x_indices(x)
        y = np.empty(x.shape, np.float)
        y[hiP_xi] = np.polyder(hiP_poly, 2)(gammapoly._get_highP_fit_x(x[hiP_xi]))
        y[loP_xi] = np.polyder(loP_poly, 2)(gammapoly._get_lowP_fit_x(x[loP_xi]))
        np.testing.assert_array_equal(gammapoly.second_derivative(x), y)

    def _test_gammapoly_integral(self, degree, rho0, x, rho_is_density=True):
        gammapoly = self._make_GammaPoly(degree, rho_is_density=rho_is_density, rho0=rho0)
        hiP_poly, _ = self._make_poly_and_setpoints(degree, [])
        loP_poly, _ = self._make_poly_and_setpoints(degree, [])
        gammapoly._hiP_fitter._set_poly(hiP_poly)  # artificially set these as the fits
        gammapoly._loP_fitter._set_poly(loP_poly)
        hiP_xi, loP_xi = gammapoly._get_highP_lowP_x_indices(x)
        y = np.empty(x.shape, np.float)
        y[hiP_xi] = np.polyint(hiP_poly)(gammapoly._get_highP_fit_x(x[hiP_xi]))
        y[loP_xi] = np.polyint(loP_poly)(gammapoly._get_lowP_fit_x(x[loP_xi]))
        np.testing.assert_array_equal(gammapoly.integral(x), y)

    def test_GammaPoly_name_prefix_class(self):
        self.assertEqual(cf.GammaPoly.name_prefix, 'gammapoly')

    def test_GammaPoly_name_prefix_obj(self):
        p = self._make_random_GammaPoly()
        self.assertEqual(p.name_prefix, 'gammapoly')

    def test_GammaPoly___init___order(self):
        degree = self._get_random_degree()
        gpoly = self._make_GammaPoly(degree)
        self.assertEqual(gpoly._order, degree, 'order is incorrect')
        self.assertEqual(gpoly._hiP_fitter._order, degree, 'High pressure order is incorrect')
        self.assertEqual(gpoly._loP_fitter._order, degree, 'Low pressure order is incorrect')

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

    def test_GammaPoly_fit_to_points_density(self):
        r0 = 5
        degree = 3
        gammapoly, hiP_poly, hiP_pts, loP_poly, loP_pts = self._make_gammapoly_and_points(degree, r0, x_min=1, x_max=10)
        points = np.concatenate((loP_pts, hiP_pts))
        gammapoly.fit_to_points(points)
        np.testing.assert_array_almost_equal(gammapoly._hiP_fitter._f.c, hiP_poly.c, decimal=0, err_msg='Bad high pressure fit')
        np.testing.assert_array_almost_equal(gammapoly._loP_fitter._f.c, loP_poly.c, decimal=0, err_msg='Bad low pressure fit')

    def test_GammaPoly_func_density_sans_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        x = x[x != r0]  # remove rho0 if it's there
        self._test_gammapoly_func(degree, r0, x)


    def test_GammaPoly_func_density_with_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        ri = np.searchsorted(x, r0)
        if x[ri] != r0:      # add rho0 if it's not there
            np.insert(x, ri, r0)
        self._test_gammapoly_func(degree, r0, x)

    def test_GammaPoly_derivative_density_sans_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        x = x[x != r0]  # remove rho0 if it's there
        self._test_gammapoly_derivative(degree, r0, x)

    def test_GammaPoly_derivative_density_with_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        ri = np.searchsorted(x, r0)
        if x[ri] != r0:  # add rho0 if it's not there
            np.insert(x, ri, r0)
        self._test_gammapoly_derivative(degree, r0, x)

    def test_GammaPoly_second_derivative_density_sans_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        x = x[x != r0]  # remove rho0 if it's there
        self._test_gammapoly_second_derivative(degree, r0, x)

    def test_GammaPoly_second_derivative_density_with_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        ri = np.searchsorted(x, r0)
        if x[ri] != r0:  # add rho0 if it's not there
            np.insert(x, ri, r0)
        self._test_gammapoly_second_derivative(degree, r0, x)

    def test_GammaPoly_integral_density_sans_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        x = x[x != r0]  # remove rho0 if it's there
        self._test_gammapoly_integral(degree, r0, x)

    def test_GammaPoly_integral_density_with_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        ri = np.searchsorted(x, r0)
        if x[ri] != r0:  # add rho0 if it's not there
            np.insert(x, ri, r0)
        self._test_gammapoly_integral(degree, r0, x)


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

    def test_GammaPoly_fit_to_points_volume(self):
        r0 = 5; degree = 3
        gammapoly, hiP_poly, hiP_pts, loP_poly, loP_pts = \
            self._make_gammapoly_and_points(degree, r0, x_min=1, x_max=10, rho_is_density=False)
        points = np.concatenate((loP_pts, hiP_pts))  # order of points shouldn't matter
        gammapoly.fit_to_points(points)
        np.testing.assert_array_almost_equal(gammapoly._hiP_fitter._f.c, hiP_poly.c, decimal=0, err_msg='Bad high pressure fit')
        np.testing.assert_array_almost_equal(gammapoly._loP_fitter._f.c, loP_poly.c, decimal=0, err_msg='Bad low pressure fit')

    def test_GammaPoly_fit_to_points_func_volume_sans_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        x = x[x != r0]      # remove rho0 if it's there
        self._test_gammapoly_func(degree, r0, x, rho_is_density=False)


    def test_GammaPoly_fit_to_points_func_volume_with_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        ri = np.searchsorted(x, r0)
        if x[ri] != r0:      # add rho0 if it's not there
            np.insert(x, ri, r0)
        self._test_gammapoly_func(degree, r0, x, rho_is_density=False)

    def test_GammaPoly_func_volume_sans_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        x = x[x != r0]  # remove rho0 if it's there
        self._test_gammapoly_func(degree, r0, x, rho_is_density=False)


    def test_GammaPoly_func_volume_with_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        ri = np.searchsorted(x, r0)
        if x[ri] != r0:      # add rho0 if it's not there
            np.insert(x, ri, r0)
        self._test_gammapoly_func(degree, r0, x, rho_is_density=False)

    def test_GammaPoly_derivative_volume_sans_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        x = x[x != r0]  # remove rho0 if it's there
        self._test_gammapoly_derivative(degree, r0, x, rho_is_density=False)

    def test_GammaPoly_derivative_volume_with_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        ri = np.searchsorted(x, r0)
        if x[ri] != r0:  # add rho0 if it's not there
            np.insert(x, ri, r0)
        self._test_gammapoly_derivative(degree, r0, x, rho_is_density=False)

    def test_GammaPoly_second_derivative_volume_sans_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        x = x[x != r0]  # remove rho0 if it's there
        self._test_gammapoly_second_derivative(degree, r0, x, rho_is_density=False)

    def test_GammaPoly_second_derivative_volume_with_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        ri = np.searchsorted(x, r0)
        if x[ri] != r0:  # add rho0 if it's not there
            np.insert(x, ri, r0)
        self._test_gammapoly_second_derivative(degree, r0, x, rho_is_density=False)

    def test_GammaPoly_integral_volume_sans_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        x = x[x != r0]  # remove rho0 if it's there
        self._test_gammapoly_integral(degree, r0, x, rho_is_density=False)

    def test_GammaPoly_integral_volume_with_rho0(self):
        r0 = 5
        degree = 3
        x = np.linspace(0.1, 10, 1000)
        ri = np.searchsorted(x, r0)
        if x[ri] != r0:  # add rho0 if it's not there
            np.insert(x, ri, r0)
        self._test_gammapoly_integral(degree, r0, x, rho_is_density=False)





if __name__ == '__main__':
    ut.main()
