import numpy as np

from astropy.modeling.core import Fittable1DModel
from astropy.modeling.parameters import Parameter
from astropy.units import Quantity


TWOPI = 2 * np.pi


class FourierSeries1D(Fittable1DModel):
    _param_names = ()

    def __init__(self, frequencies,
                 n_models=None, model_set_axis=None, name=None, meta=None):
        self._frequencies = np.array(frequencies)
        self._param_names = self._generate_coeff_names()

        for param_name in self._param_names:
            self._parameters_[param_name] = Parameter(param_name, default=np.zeros(()))

        super().__init__(
            n_models=n_models, model_set_axis=model_set_axis, name=name,
            meta=meta)

    @property
    def param_names(self):
        return self._param_names

    @property
    def n_terms(self):
        return len(self._frequencies)

    def _generate_coeff_names(self):
        names = []

        for index in range(self.n_terms):
            names.append(f"a{index}")
            names.append(f"b{index}")

        return tuple(names)

    @staticmethod
    def _evaluate_term(x, freq, coeff_a, coeff_b):

        argument = TWOPI * (freq * x)
        if isinstance(argument, Quantity):
            argument = argument.value

        return coeff_a * np.sin(argument) + coeff_b * np.cos(argument)

    def evaluate(self, x, *coeffs):
        coefficients = list(coeffs)

        value = 0
        for index in range(self.n_terms):
            coeff_a = coefficients.pop(0)
            coeff_b = coefficients.pop(0)

            value += self._evaluate_term(x, self._frequencies[index], coeff_a, coeff_b)

        return value
