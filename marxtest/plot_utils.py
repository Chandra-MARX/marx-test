import numpy as np
from numpy import ma
from matplotlib import scale as mscale
from matplotlib import transforms as mtransforms
from matplotlib.ticker import AutoLocator, ScalarFormatter


class PowerScale(mscale.ScaleBase):
    """
    Scales values v with with v^pow.
    """
    name = 'power'

    def __init__(self, axis, **kwargs):
        """
        power: The p in x**p.
        """
        mscale.ScaleBase.__init__(self)
        self.power = kwargs.pop("power", 1.)

    def get_transform(self):
        """
        Override this method to return a new instance that does the
        actual transformation of the data.
        """
        return self.PowerTransform(self.power)

    def set_default_locators_and_formatters(self, axis):
        """
        Override to set up the locators and formatters to use with the
        scale.
        """
        axis.set_major_locator(AutoLocator())
        axis.set_major_formatter(ScalarFormatter())

    def limit_range_for_scale(self, vmin, vmax, minpos):
        """
        Override to limit the bounds of the axis to the domain of the
        transform. For power < 1, negative values are not allowed.
        Unlike the
        autoscaling provided by the tick locators, this range limiting
        will always be adhered to, whether the axis range is set
        manually, determined automatically or changed through panning
        and zooming.
        """
        if self.power >= 1:
            return vmin, vmax
        else:
            return max(vmin, 0), max(vmax, 0)

    class PowerTransform(mtransforms.Transform):
        # There are two value members that must be defined.
        # ``input_dims`` and ``output_dims`` specify number of input
        # dimensions and output dimensions to the transformation.
        # These are used by the transformation framework to do some
        # error checking and prevent incompatible transformations from
        # being connected together.  When defining transforms for a
        # scale, which are, by definition, separable and have only one
        # dimension, these members should always be set to 1.
        input_dims = 1
        output_dims = 1
        is_separable = True

        def __init__(self, power):
            mtransforms.Transform.__init__(self)
            self.power = power

        def transform_non_affine(self, a):
            """
            This transform takes an Nx1 ``numpy`` array and returns a
            transformed copy.
            """
            return a**self.power

        def inverted(self):
            """
            Override this method so matplotlib knows how to get the
            inverse transform for this transform.
            """
            return PowerScale.PowerTransform(1./self.power)

# Now that the Scale class has been defined, it must be registered so
# that ``matplotlib`` can find it.
mscale.register_scale(PowerScale)
