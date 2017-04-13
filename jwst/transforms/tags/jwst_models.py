
from __future__ import absolute_import, division, unicode_literals, print_function
import numpy as np
from asdf import yamlutil
from asdf.tags.transform.basic import TransformType

from ..models import (WavelengthFromGratingEquation, AngleFromGratingEquation,
                      NRSZCoord, Unitless2DirCos, DirCos2Unitless, Rotation3DToGWA,
                      LRSWavelength, Gwa2Slit, Slit2Msa, Logical, NirissSOSSModel, V23ToSky,
                      RefractionIndexFromPrism, Snell, MIRI_AB2Slice)


__all__ = ['GratingEquationType', 'CoordsType', 'RotationSequenceType', 'LRSWavelengthType',
           'Gwa2SlitType', 'Slit2MsaType', 'LogicalType', 'NirissSOSSType', 'V23ToSky',
           'RefractionIndexType', 'SnellType', 'MIRI_AB2SliceType']


class RotationSequenceType(TransformType):
    name = "rotation_sequence"
    types = [Rotation3DToGWA]
    standard = "jwst_pipeline"
    version = "0.7.0"

    @classmethod
    def from_tree_transform(cls, node, ctx):
        angles = node['angles']
        axes_order = node['axes_order']
        return Rotation3DToGWA(angles, axes_order)

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'angles': list(model.angles.value)}
        node['axes_order'] = model.axes_order
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)


class V23ToSkyType(TransformType):
    name = "v23tosky"
    types = [V23ToSky]
    standard = "jwst_pipeline"
    version = "0.7.0"

    @classmethod
    def from_tree_transform(cls, node, ctx):
        angles = node['angles']
        axes_order = node['axes_order']
        return V23ToSky(angles, axes_order=axes_order)

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'angles': list(model.angles.value)}
        node['axes_order'] = model.axes_order
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)


class CoordsType(TransformType):
    name = "coords"
    types = [Unitless2DirCos, DirCos2Unitless, NRSZCoord]
    standard = "jwst_pipeline"
    version = "0.7.0"

    @classmethod
    def from_tree_transform(cls, node, ctx):
        model_type = node['model_type']
        if model_type == 'nrszcoord':
            return NRSZCoord()
        elif model_type == 'to_dircos':
            return Unitless2DirCos()
        elif model_type == 'from_dircos':
            return DirCos2Unitless()
        else:
            raise TypeError("Unknown model_type")

    @classmethod
    def to_tree_transform(cls, model, ctx):
        if isinstance(model, NRSZCoord):
            model_type = 'nrszcoord'
        elif isinstance(model, DirCos2Unitless):
            model_type = 'from_dircos'
        elif isinstance(model, Unitless2DirCos):
            model_type = 'to_dircos'
        else:
            raise TypeError("Model of type {0} i snot supported.".format(model.__class__))
        node = {'model_type': model_type}
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)


class GratingEquationType(TransformType):
    name = "grating_equation"
    types = [WavelengthFromGratingEquation, AngleFromGratingEquation]
    standard = "jwst_pipeline"
    version = "0.7.0"

    @classmethod
    def from_tree_transform(cls, node, ctx):
        groove_density = node['groove_density']
        order = node['order']
        output = node['output']
        if output == "wavelength":
            model = WavelengthFromGratingEquation(groove_density, order)
        elif output == "angle":
            model = AngleFromGratingEquation(groove_density, order)
        else:
            raise ValueError("Can't create a GratingEquation model with output {0}".format(output))
        return model

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'order': model.order.value,
                'groove_density': model.groove_density.value
                }
        if isinstance(model, AngleFromGratingEquation):
            node['output'] = 'angle'
        elif isinstance(model, WavelengthFromGratingEquation):
            node['output'] = 'wavelength'
        else:
            raise TypeError("Can't serialize an instance of {0}".format(model.__class__.__name__))
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)

    #@classmethod
    #def assert_equal(cls, a, b):
        #from astropy import modeling

        #TransformType.assert_equal(a, b)
        #assert (isinstance(a, modeling.models.Shift) and
                #isinstance(b, modeling.models.Shift))
        #assert_array_equal(a.offset.value, b.offset.value)


class LRSWavelengthType(TransformType):
    name = "lrs_wavelength"
    types = [LRSWavelength]
    standard = "jwst_pipeline"
    version = "0.7.0"

    @classmethod
    def from_tree_transform(cls, node, ctx):
        wavetable = node['wavetable']
        zero_point = node['zero_point']
        return LRSWavelength(wavetable, zero_point)

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'wavetable': model.wavetable,
                'zero_point': model.zero_point}
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)


class Gwa2SlitType(TransformType):
    name = "gwa_to_slit"
    types = [Gwa2Slit]
    standard = "jwst_pipeline"
    version = "0.7.0"

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return Gwa2Slit(node['slits'], node['models'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'slits': model._slits,
                'models': model.models}
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)


class Slit2MsaType(TransformType):
    name = "slit_to_msa"
    types = [Slit2Msa]
    standard = "jwst_pipeline"
    version = "0.7.0"

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return Slit2Msa(node['slits'], node['models'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'slits': model._slits,
                'models': model.models}
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)


class LogicalType(TransformType):
    name = "logical"
    types = [Logical]
    standard = "jwst_pipeline"
    version = "0.7.0"

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return Logical(node['condition'], node['compareto'], node['value'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'condition': model.condition,
                'compareto': model.compareto,
                'value': model.value}
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)


class NirissSOSSType(TransformType):
    name = "niriss_soss"
    types = [NirissSOSSModel]
    standard = "jwst_pipeline"
    version = "0.7.0"

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return NirissSOSSModel(node['spectral_orders'], node['models'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'spectral_orders': list(model.models.keys()),
                'models': list(model.models.values())
                }
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)


class RefractionIndexType(TransformType):
    name = "refraction_index_from_prism"
    types = [RefractionIndexFromPrism]
    standard = "jwst_pipeline"
    version = "0.7.0"

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return RefractionIndexFromPrism(node['prism_angle'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'prism_angle': model.prism_angle.value
                }
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)

    @classmethod
    def assert_equal(cls, a, b):
        TransformType.assert_equal(a, b)
        assert (isinstance(a, RefractionIndexFromPrism) and
                isinstance(b, RefractionIndexFromPrism))
        assert_array_equal(a.prism_angle, b.prism_angle)


class SnellType(TransformType):
    name = "snell"
    types = [Snell]
    standard = "jwst_pipeline"
    version = "0.7.0"

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return Snell(node['prism_angle'], node['kcoef'], node['lcoef'], node['tcoef'],
                     node['ref_temp'], node['ref_pressure'], node['temp'], node['pressure'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'prism_angle': model.prism_angle,
                'kcoef': model.kcoef.tolist(),
                'lcoef': model.lcoef.tolist(),
                'tcoef': model.tcoef.tolist(),
                'ref_temp': model.tref,
                'ref_pressure': model.pref,
                'temp': model.temp,
                'pressure': model.pressure
                }
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)

    @classmethod
    def assert_equal(cls, a, b):
        TransformType.assert_equal(a, b)
        assert (isinstance(a, Snell) and
                isinstance(b, Snell))
        assert_array_equal(a.prism_angle, b.prism_angle)
        assert_array_equal(a.kcoef, b.kcoef)
        assert_array_equal(a.lcoef, b.lcoef)
        assert_array_equal(a.tcoef, b.tcoef)
        assert_array_equal(a.tref, b.tref)
        assert_array_equal(a.pref, b.pref)
        assert_array_equal(a.temp, b.temp)
        assert_array_equal(a.pressure, b.pressure)


class MIRI_AB2SliceType(TransformType):
    name = "miri_ab2slice"
    types = [MIRI_AB2Slice]
    standard = "jwst_pipeline"
    version = "0.7.0"

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return MIRI_AB2Slice(node['beta_zero'], node['beta_del'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'beta_zero': model.beta_zero.value, 'beta_del': model.beta_del.value}
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)

    @classmethod
    def assert_equal(cls, a, b):
        TransformType.assert_equal(a, b)
        assert (isinstance(a, MIRI_AB2Slice) and
                isinstance(b, MIRI_AB2Slice))
        assert_array_equal(a.beta_zero, b.beta_zero)
        assert_array_equal(a.beta_del, b.beta_del)
