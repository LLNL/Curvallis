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

    @staticmethod
    def makeArgs(names, values):
        Args = namedtuple('Args', names)
        return Args(*values)

    @staticmethod
    def makePolyArgs(derivative_scale=1, second_derivative_scale=1, integral_scale=1,
                     x_integral_ref=1, y_integral_ref=1):
        argnames = inspect.getfullargspec(TestCurveFitters.makePolyArgs).args
        l = locals()  # for some reason, Python won't let me use locals() directly in the next statement
        argvalues = [l[a] for a in argnames]
        return TestCurveFitters.makeArgs(argnames, argvalues)

    @staticmethod
    def make_poly_and_points(degree, x_min=0, x_max=100, num_pts=0):
        """ Generates a random polynomial of the specified degree and a number of points along that curve

        :param num_pts: the number of points to generate.
                    If <=0, will default to degree * (degree + 1)
        :return: poly, points
                    where poly is a numpy.poly1d object used to generate the points
                    and points is a list of (x,y) tuples with points on the curve represented by poly
        """
        num_pts = int(num_pts)
        if num_pts <= 0:
            num_pts = degree**2 + degree
        x = np.linspace(x_min, x_max, num_pts)
        poly = np.poly1d(np.random.randint(10, size=degree+1))
        return poly, list(zip(x, poly(x)))

    def test_Poly5_fit_to_points(self):
        deg = 5
        actual_f, pts = self.make_poly_and_points(deg)
        poly5 = cf.Poly_Original(self.makePolyArgs(), 'poly{0}'.format(deg))
        poly5.fit_to_points(pts)
        np.testing.assert_array_almost_equal(poly5._f.c, actual_f.c, decimal=5)

    @ut.skip('fitter often comes up with completely different polynomial with positive and negative coefficients')
    def test_Poly5_fit_to_points_noisy(self):
        deg = 5
        actual_f, pts = self.make_poly_and_points(deg, num_pts=200)
        # add ~2% noise to the y values of the points
        noisy_pts = [(x, y*n) for ((x, y), n) in zip(pts, np.random.normal(1, .02, len(pts)))]
        poly5 = cf.Poly_Original(self.makePolyArgs(), 'poly{0}'.format(deg))
        poly5.fit_to_points(noisy_pts)
        np.testing.assert_array_almost_equal(poly5._f.c, actual_f.c, decimal=1)



if __name__ == '__main__':
    ut.main()


