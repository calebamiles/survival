'''
Created on Oct 31, 2012

@author: jason
'''
import unittest
import pymc
import numpy
from src.survival import ConstantHazardEventObservation, EventProposer, MultiEventProposer,\
    ConstantHazardMultiEventObservation

pymc.StepMethodRegistry.append(EventProposer)
pymc.StepMethodRegistry.append(MultiEventProposer)

numpy.random.seed(0)

class Test(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass


    def testConstantHazardEventObservation(self):
#        params = {'iter':100000, 'burn':30000, 'thin':500}
        hazard_val = 0.1
        params = {'iter':10000, 'burn':300, 'thin':50}
        hazard = pymc.Uniform('hazard',lower=0,upper=100,value=hazard_val,observed=True)
        observation_time = pymc.Uniform('observation_time',lower=0,upper=100,value=100*numpy.ones(1000),observed=True)
        event = ConstantHazardEventObservation('event',hazard=hazard,period=observation_time,observed=False)
        vars = {hazard, observation_time, event}
        mcmc = pymc.MCMC(vars)
        mcmc.sample(**params)
        print ''
        events = mcmc.trace('event')[:]
        hazard = pymc.Uniform('hazard',lower=0,upper=100,value=1,observed=False)
        observation_time = pymc.Uniform('observation_time',lower=0,upper=100,value=100*numpy.ones(1000),observed=True)
        event = ConstantHazardEventObservation('event',hazard=hazard,period=observation_time,value=events,observed=True)
        vars = {hazard, observation_time, event}
        mcmc = pymc.MCMC(vars)
        mcmc.sample(**params)
        fitted_hazard_val = numpy.mean(mcmc.trace('hazard')[:])
        self.assertTrue(abs(fitted_hazard_val - hazard_val) < 0.01)
        #Uncomment below to plot trace to hazard.png
#        pymc.Matplot.plot(mcmc.trace('hazard'))

    def testConstantHazardMultiEventObservation(self):
#        params = {'iter':100000, 'burn':30000, 'thin':500}
        params = {'iter':1000, 'burn':300, 'thin':50}#TODO: The number of iterations should be much higher in practice.
        hazard_val = 0.1
        hazard = pymc.Uniform('hazard',lower=0,upper=100,value=hazard_val,observed=True)
        observation_time = pymc.Uniform('observation_time',lower=0,upper=100,value=100*numpy.ones(1000),observed=True)
        event = ConstantHazardMultiEventObservation('event',hazard=hazard,period=observation_time,observed=False)
        vars = {hazard, observation_time, event}
        mcmc = pymc.MCMC(vars)
        mcmc.sample(**params)
        events = mcmc.trace('event')[:]
        hazard = pymc.Uniform('hazard',lower=0,upper=100,value=1,observed=False)
        observation_time = pymc.Uniform('observation_time',lower=0,upper=100,value=100*numpy.ones(1000),observed=True)
        event = ConstantHazardMultiEventObservation('event',hazard=hazard,period=observation_time,value=events,observed=True)
        vars = {hazard, observation_time, event}
        mcmc = pymc.MCMC(vars)
        mcmc.sample(**params)
        fitted_hazard_val = numpy.mean(mcmc.trace('hazard')[:])
        self.assertTrue(abs(fitted_hazard_val - hazard_val) < 0.01)
        #Uncomment below to plot trace to hazard.png
        pymc.Matplot.plot(mcmc.trace('hazard'))
        self.assertTrue(abs(fitted_hazard_val - hazard_val) < 0.01)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testConstantHazard']
    unittest.main()