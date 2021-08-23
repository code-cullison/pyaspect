'''
This module adds functionality to the MomentTensor Class
in the pyrocko Python package.  I added featuress to make
it easier to use the moment-tensors in conjunction with
SPECFEM3D modeling, inversion, and for training neural-
network training
'''

__version__ = '0.5'
__author__ = 'Thomas Cullison @ Utrecht University'


import numpy as np

import pyrocko.moment_tensor as pmt
from pyrocko.moment_tensor import MomentTensor as RockoMT


################################################
## MomentTensor class helper functions

def aki_moment_to_magnitude(moment):

    '''
    Convert scalar moment to [Aki & Richards] magnitude (Mw).
    
    :param moment: scalar moment [Nm]
    :returns: moment magnitude Mw

    Moment magnitude is defined as
    .. math::
        M_\\mathrm{w} = {\\frac{2}{3}}\\(log_{10}(M_0) - 9.1)
    where :math:`M_0` is the scalar moment given in [Nm].
    .. note::
    '''

    return (np.log10(moment)-9.1)/1.5


def aki_dcm_moment_to_magnitude(moment):

    '''
    Convert scalar moment to [Aki & Richards] magnitude (Mw).
    
    :param moment: scalar moment [dyne-cm]
    :returns: moment magnitude Mw

    Moment magnitude is defined as
    .. math::
        M_\\mathrm{w} = {\\frac{2}{3}}\\(log_{10}(M_0) - 9.1)
    where :math:`M_0` is the scalar moment given in [Nm].
    .. note::
    '''

    return (np.log10(moment)-16.1)/1.5


def aki_magnitude_to_moment(magnitude):
    '''
    Convert [Aki & Richards] magnitude moment (Mw) to scalar moment.
    
    :param magnitude: moment magnitude
    :returns: scalar moment [Nm]

    See :py:func:`aki_moment_to_magnitude`.
    '''

    return 10.0**(1.5*(magnitude)+9.1)


def aki_dcm_magnitude_to_moment(magnitude):
    '''
    Convert [Aki & Richards] magnitude moment (Mw) to scalar moment.
    
    :param magnitude: moment magnitude
    :returns: scalar moment [dyne-cm]

    See :py:func:`aki_moment_to_magnitude`.
    '''

    #return 10.0**(1.5*(magnitude)+16.1)
    return aki_magnitude_to_moment(magnitude)*10.0**7


#####################################################################
#                                                
#  Class Definition:                           
#  Note: inheriting MomentTensor from pyrocko. However, 
#  pyrocko uses a slightly different Mw to M0 transformation
#  and a different naming convention for tensors and components.
#                                              
#####################################################################

class MomentTensor(RockoMT):

    """
    Here we wrap the MomentTensor object from the pyrocko package.

    We just need a few features which is why we do not inherit the object.
    
    See documentation at https://pyrocko.org/

    :param strike,dip,rake: angles in [degrees]
    :param scalar_moment: scalar moment (i.e. M0) in as N*m
    """
    

    #def __init__(self,mw=None,strike=None,dip=None,rake=None):
    def __init__(self,**kwargs):

        self._arg_strike = 0.
        if 'strike' in kwargs.keys():
            self._arg_strike = kwargs['strike']

        self._arg_dip = 0.
        if 'dip' in kwargs.keys():
            self._arg_dip = kwargs['dip']

        self._arg_rake = 0.
        if 'rake' in kwargs.keys():
            self._arg_rake = kwargs['rake']

        super(MomentTensor, self).__init__(**kwargs) #pyrocko.moment_tensor.MomentTensor() 


    @classmethod
    def from_values(cls, values):
        pyrocko_m = pmt.values_to_matrix(values) #little trick
        return MomentTensor(m=pyrocko_m)
    


    ######################################
    ## Accessors

    @property
    def strike(self): # the strike as given by argument
        return self._arg_strike

    @property
    def dip(self): # the dip as given by argument
        return self._arg_dip

    @property
    def rake(self): # the rake as given by argument
        return self._arg_rake

    '''
    @property
    def Mw(self): # the rake as given by argument
        return self._arg_Mw
    '''
    

    ######################
    #                    #
    # Aki & Richards CMT #
    #                    #
    ######################

    
    @property
    def aki_moment(self): # same as for the harvard CMT
        return aki_magnitude_to_moment(self.magnitude)

    @property 
    def aki_dcm_moment(self): # same as for the harvard CMT
        return aki_dcm_magnitude_to_moment(self.magnitude)

    #for some reason IPython doesn't recognize super() function
    def m_east_north_up(self): 
        return self._to_east_north_up.T * self._m * self._to_east_north_up

    #for some reason IPython doesn't recognize super() function
    def m6_east_north_up(self):
        return pmt.to6(self.m_east_north_up())


    def aki_richards_unit_m6(self):
        """
          returns a [unit] moment-tensor array based on Aki & Richards Mw :-> M0
        """
        return self.m6_east_north_up()/self.moment


    def aki_richards_m6(self):
        """
          returns a [Nm] moment-tensor array based on Aki & Richards Mw :-> M0 
        """
        return self.aki_moment*self.aki_richards_unit_m6()


    def aki_richards_dcm_m6(self):
        """
          returns a [dyne-cm] moment-tensor array based on Aki & Richards Mw :-> M0 
        """
        return self.aki_dcm_moment*self.aki_richards_unit_m6()


    def aki_richards_unit_matrix(self):
        """
          returns a [unit] moment-tensor matrix based on Aki & Richards Mw :-> M0 
        """
        return  self.m()/self.moment


    def aki_richards_matrix(self):
        """
          returns a [Nm] moment-tensor matrix based on Aki & Richards Mw :-> M0 
        """
        return self.aki_moment*self.aki_richards_unit_matrix()


    def aki_richards_dcm_matrix(self):
        """
          returns a [dyne-cm] moment-tensor matrix based on Aki & Richards Mw :-> M0 
        """
        return self.aki_dcm_moment*self.aki_richards_unit_matrix()



    ###############
    #             #
    # Harvard CMT #
    #             #
    ###############
    

    @property
    def harvard_moment(self): # also used for the harvard CMT
        return aki_magnitude_to_moment(self.magnitude)

    @property 
    def harvard_dcm_moment(self): # also used for the harvard CMT
        return aki_dcm_magnitude_to_moment(self.magnitude)
    

    def harvard_unit_m6(self):
        """
          returns a [unit] Harvard CMT array based on Aki & Richards Mw :-> M0 
        """
        return self.m6_up_south_east()/self.moment


    def harvard_m6(self):
        """
          returns a [Nm] Harvard CMT array based on Aki & Richards Mw :-> M0 
        """
        return self.harvard_moment*self.harvard_unit_m6()


    def harvard_dcm_m6(self):
        """
          returns a [dyne-cm] Harvard CMT array based on Aki & Richards Mw :-> M0 
        """
        return self.harvard_dcm_moment*self.harvard_unit_m6()


    def harvard_unit_matrix(self):
        """
          returns a [unit] Harvard CMT matrix based on Aki & Richards Mw :-> M0 
        """
        return self.m_up_south_east()/self.moment
        

    def harvard_matrix(self):
        """
          returns a [Nm] Harvard CMT matrix based on Aki & Richards Mw :-> M0 
        """
        return self.harvard_moment*self.harvard_unit_matrix()


    def harvard_dcm_matrix(self):
        """
          returns a [dyne-cm] Harvard CMT matrix based on Aki & Richards Mw :-> M0 
        """
        return self.harvard_dcm_moment*self.harvard_unit_matrix()
