from matplotlib import pyplot

import numpy

from scipy.special import expi

class Transient():
    """Calculates well bottomhole pressure in field units for the given
    constant flow rate oil production, single phase flow."""

    def __init__(self,rrock,phase,well,tcomp=None,immob=None):
        """Initializes reservoir, oil and well properties."""

        self.rrock = rrock
        self.phase = phase
        self.tcomp = tcomp

        self.well  = well

        self.immob = immob

    def press(self,time=None,size=50,scale='linear'):
        """It calculates delta bottom hole flowing pressure for constant
        flow rate production at transient state, exponential integral solution."""

        if time is not None:
            self.time = time[time>=self.twell]

        else:

            if scale == "linear":
                self.time = numpy.linspace(self.twell,self.tp,size)
            elif scale == "log":
                self.time = numpy.logspace(*numpy.log10([self.twell,self.tp]),size)

            self._scale = scale

        Ei = expi(-39.5*(self.rw**2)/(self.eta*self.time))

        self._delta = -self.flowterm*(1/2*Ei-self.skin)

    def set_boundary(self,radius=None,area=None):
        """
        Sets the radius of outer circular boundary, ft and calculates
        pseudo steady state solution correcting pressure for boundary effects.
        
        If area is not None, the radius is calculated from area [acre].
        """

        if area is not None:
            radius = numpy.sqrt(area*43560/numpy.pi)

        self._radius = radius

        pseudo_cond = self.time>=self.tbound

        pseudo_time = self.time[pseudo_cond]

        term_time = -0.012648*(self.eta*pseudo_time)/self.radius**2

        term_boundary = -numpy.log(self.radius/self.rw)

        delta = -self.flowterm*(term_time+term_boundary+3/4-self.skin)

        self._delta[pseudo_cond] = delta

    def view(self,axis=None,pinitial=None,scale=None):

        showFlag = True if axis is None else False

        if axis is None:
            figure,axis = pyplot.subplots(nrows=1,ncols=1)

        if pinitial is None:
            yaxis,ylabel = -self.delta,f"Wellbore Pressure Drop [psi]"
        else:
            yaxis,ylabel = pinitial+self.delta,"Wellbore Pressure [psi]"

        axis.plot(self.time,yaxis)

        if scale is None:
            try:
                scale = self._scale
            except AttributeError:
                scale = "linear"

        axis.set_xscale(scale)

        if pinitial is None:
            axis.invert_yaxis()

        axis.set_xlabel("Time [days]")
        axis.set_ylabel(ylabel)

        if showFlag:
            pyplot.show()

    def set_dimensionless(self,pinitial):

        self.pD = (self.rrock.perm*self.rrock.height)/(141.2*self.well.rate/self.fluid.mobil)*(pinitial-self.delta)

        self.tD = (0.00632829*self.rrock.perm*self.time)/(self.rrock.poro*self.fluid.visc*self.ccomp*self.well.radius**2)

        self.CD = (0.8936*self.Cs)/(self.rrock.poro*self.ccomp*self.rrock.height*self.well.radius**2)

    @property
    def ccomp(self):
        """Sets the total compressibility. If None,
        it will try to calculate it from rock and fluid compressibilities."""

        if self.tcomp is not None:
            return self.tcomp

        irrcomp,irrsat = 0.,0.

        if self.immob is not None:
            irrcomp,irrsat = self.immob.comp,self.immob.sat
            
        return self.rrock.comp+self.phase.comp*(1-irrsat)+irrcomp*irrsat

    @property
    def eta(self):
        return self.rrock.perm/(self.rrock.poro*self.fluid.visc*self.ccomp)

    @property
    def storage(self):
        return (144*self.well.Awb)/(5.615*self.fluid.rho)

    @property
    def flowterm(self):
        return -141.2*(self.well.rate)/(self.fluid.mobil*self.rrock.perm*self.rrock.height)

    @property
    def radius(self):
        return self._radius

    @property
    def delta(self):
        return self._delta

    @property
    def twell(self):
        return 15802*self.well.radius**2/self.eta

    @property
    def tedge(self):
        return 39.5*self.rrock.radius**2/self.eta

if __name__ == "__main__":

    print(finite.setnodes(10,r=2.5))

    # tD = 1e-2

    # umin = numpy.sqrt(1/tD)/1e6
    # umax = numpy.sqrt(1/tD)*1e5

    # u = numpy.logspace(numpy.log10(umin),numpy.log10(umax),2000)

    # y1 = everdingen.pressure_integrand(u,tD)
    # y2 = everdingen.pressure_integrand(u,1e8)
    # y1 = everdingen.pressure_integrand(u,0.0001)

    # y1 = agarwal.pressure_integrand(u,1e8,1e5,0)
    # y2 = agarwal.pressure_integrand(u,1e8,1e5,0)
    # y2 = agarwal.pressure_lineSource_integrand(u,1e3,10000,0)
    # y3 = agarwal.pressure_lineSource_integrand(u,1e3,10000,5)
    # y4 = agarwal.pressure_lineSource_integrand(u,1e3,10000,10)
    # y5 = agarwal.pressure_lineSource_integrand(u,1e3,10000,20)

    # print(f"{agarwal.pressure(1e8,1e5,0)[0]:8.5f}")

    # pyplot.loglog(u,y1)
    # pyplot.loglog(u,y1)
    # pyplot.loglog(u,y2)
    # pyplot.loglog(u,y3)
    # pyplot.loglog(u,y4)
    # pyplot.loglog(u,y5)

    # pyplot.show()



