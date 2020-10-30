import unittest as ut
from unittest.mock import patch

import curvallis.curve_editing.test.test_utilities as util
from curvallis.curve_editing import io,  curve_fitters as cf
from curvallis.curve_editing.configargparse import ArgumentParser



class MyTestCase(ut.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def _get_args(self, parser, **kwargs):
        args = dict(parser._defaults)
        args.update({'no_config_file': False,
                     'input_file': 'foo',
                     'rho0': None,
                     'region_bound': None,
                     'numpoints': 1})
        args.update(kwargs)
        return util.make_args(**args)

    def test_process_args_gammapoly_no_rho0_guess(self):
        with patch.object(ArgumentParser, 'error', return_value=None) as mock_error:
            parser = ArgumentParser()
            io.define_args(parser)
            args = self._get_args(parser, **{'fit_type': [cf.GammaPoly.name_prefix+'2']})
            io.process_args(parser, args)
        mock_error.assert_called_once_with('If using fitter "sandiapc", "gammapoly", or "gammapolyv", you must give a value for "rho0_guess".')


    def test_process_args_gammapoly_with_rho0_guess(self):
        with patch.object(ArgumentParser, 'error', return_value=None) as mock_error:
            parser = ArgumentParser()
            io.define_args(parser)
            args = self._get_args(parser, **{'fit_type': [cf.GammaPoly.name_prefix+'2'], 'rho0': 5})
            io.process_args(parser, args)
        mock_error.assert_not_called()

    def test_process_args_gammapolyv_no_rho0_guess(self):
        with patch.object(ArgumentParser, 'error', return_value=None) as mock_error:
            parser = ArgumentParser()
            io.define_args(parser)
            args = self._get_args(parser, **{'fit_type': [cf.GammaPolyV.name_prefix + '2']})
            io.process_args(parser, args)
        mock_error.assert_called_once_with(
            'If using fitter "sandiapc", "gammapoly", or "gammapolyv", you must give a value for "rho0_guess".')

    def test_process_args_gammapolyv_with_rho0_guess(self):
        with patch.object(ArgumentParser, 'error', return_value=None) as mock_error:
            parser = ArgumentParser()
            io.define_args(parser)
            args = self._get_args(parser, **{'fit_type': [cf.GammaPolyV.name_prefix + '2'], 'rho0': 5})
            io.process_args(parser, args)
        mock_error.assert_not_called()







if __name__ == '__main__':
    ut.main()
