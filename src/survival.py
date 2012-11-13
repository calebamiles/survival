'''
Created on Oct 31, 2012

@author: jason
'''

import pymc
from math import log
import random as rand
import numpy
from itertools import izip


class Event(object):
    '''Represents a single event that can happen once.  The typical example is death.
        time - (float) The time to death or censorship
        censored - (bool) True if the event was right censored; False otherwise. 
    '''
    def __init__(self, time, censored = False):
        self.time = time
        self.censored = censored
        
class MultiEvent(object):
    '''
    Represents a repeatable event, such as a heart attack.
        uncensored_events - (list) A list of Event objects, all of which have censored=False.
    '''
    def __init__(self, uncensored_events, decumulate_times=False):
        self.uncensored_events = uncensored_events
        assert(not any([event.censored for event in self.uncensored_events]))
        if decumulate_times:
            self.uncensored_events = self.uncensored_events.sort(key=lambda event: event.time)
            last_time = 0.0
            for event in self.uncensored_events:
                event.time = event.time - last_time
            
    def __iter__(self):
        return iter(self.uncensored_events)
    
    def __len__(self):
        return len(self.uncensored_events)
    
    @property
    def time(self):
        if self.uncensored_events:
            return sum([event.time for event in self.uncensored_events])
        else:
            return 0.0

class MultiCost(object):
    def __init__(self, costs=[]):
        self.costs = costs
    
    def __iter__(self):
        return iter(self.costs)
    
    def __len__(self):
        return len(self.costs)
    
    @property
    def cost(self):
        if self.costs:
            return sum([cost for cost in self.costs if cost is not None])
        else:
            return 0.0

class EventProposer(pymc.Metropolis):
    '''
    Simple proposal distribution for Event objects for use in Metropolis samplers.
    '''
    def __init__(self, stochastic, *args, **kwargs):
        return super(self.__class__, self).__init__(stochastic, *args, **kwargs)
    
    def propose(self):
        tau = 1./(self.adaptive_scale_factor * self.proposal_sd)**2
        time = pymc.rnormal(self.stochastic.value.time, tau)
        censored = rand.random() > 0.5
        self.stochastic.value = Event(time, censored)
    
    def competence(self, stochastic):
        if stochastic.dtype == Event:
            return 3
        else:
            return 0
    
class MultiEventProposer(pymc.Metropolis):
    '''
    Simple proposal distribution for MultiEvent objects for use in Metropolis samplers.
    '''
    def __init__(self, stochastic, *args, **kwargs):
        return super(self.__class__, self).__init__(stochastic, *args, **kwargs)
    
    def propose(self):
        tau = 1./(self.adaptive_scale_factor * self.proposal_sd)**2
        time = pymc.rnormal(self.stochastic.value.time, tau)
        n = pymc.rnormal(len(self.stochastic.value), tau)
        if n <= 0:
            n = 0
        times = [rand.random() for _ in range(n)]
        total = float(sum(times))
        times = [item*time/total for item in times]
        events = [Event(time=item, censored=False) for item in times]
        self.stochastic.value = MultiEvent(events)
    
    def competence(self, stochastic):
        if stochastic.dtype == MultiEvent:
            return 3
        else:
            return 0
    

class NegativeHazardError(Exception):
    pass


def ConstantHazardMultiEventObservation(name, hazard, period, observed = False, value = None):
    '''
    Create a pymc Stochastic representing a Poisson process.  The values produced are 
    MultiEvent objects.
    '''
    def multievent_logp(value, hazard, period):
        '''Scalar valued function of scalars.
            value - (MultiEvent) The event times
            hazard - (float) The hazard rate
            period - (float) The duration of the observation period
        '''
        result = 0
        remaining_time = period
        for event in value:
            result += log(hazard) - hazard * event.time
            remaining_time = remaining_time - event.time
        if remaining_time < 0:
            return 0
        result += -1 * hazard * remaining_time
        return result
    
    __logp = numpy.vectorize(multievent_logp,otypes=[MultiEvent])
    
    def logp(value, hazard, period):
        '''
        Array valued function of arrays.
            value - array of MultiEvent objects
            hazard - array of hazard rates
            period - array of observation period durations
        '''
        if type(value) is list:
            value = numpy.array(value)
        return numpy.sum(__logp(value, hazard, period))
    
    def multievent_random(hazard, period):
        '''
        Scalar valued function of scalars.  Generates a random MultiEvent from the
        distribution.
            hazard - (float) The hazard rate
            period - (float) The duration of the observation period
        '''
        if hazard == 0:
            return MultiEvent([])
        if hazard < 0.0:
            raise NegativeHazardError
        remaining_time = period
        events = []
        while True:
            F = rand.random()
            time = -1 * log(1-F) / float(hazard)
            if time > remaining_time:
                break
            remaining_time -= time
            events.append(Event(time=time,censored=False))
        return MultiEvent(events)

    __random = numpy.vectorize(multievent_random,otypes=[MultiEvent])
    
    def random(hazard, period):
        '''
        Array valued function of arrays.  Generated random MultiEvents from the 
        distribution.
            hazard - Array of floats representing hazard rates
            period - Array of floats representing observation period durations
        '''
        return  __random(hazard, period)

    dtype = MultiEvent
    
    #Create the pymc Stochastic object
    result = pymc.Stochastic(logp = logp,
                             doc = 'A constant hazard survival time node for multiple events.',
                             name = name,
                             parents = {'hazard':hazard,'period':period},
                             random = random,
                             trace = True,
                             value = value,
                             dtype = dtype,
                             rseed = 1,
                             observed = observed,
                             cache_depth = 2,
                             plot = True,
                             verbose = 0)
    return result

def ConstantHazardEventObservation(name, hazard, period, observed = False, value = None):
    '''
    Create a pymc Stochastic representing a constant hazard failure process.  The values
    produced are Event objects.
    '''
    
    def logp(value, hazard, period):
        '''
        Array valued function of arrays.
            value - Array of event objects
            hazard - Array of floats representing hazard rates
            period - Array of floats representing observation period durations
        '''
        if type(value) is list:
            value = numpy.array(value)
        time = numpy.empty_like(value,dtype=numpy.float)
        time.flat = [val.time for val in value.flat]
        censored = numpy.empty_like(value,dtype=numpy.float)
        censored.flat = [float(val.censored) for val in value.flat]
        
        return numpy.sum((time < period) * ((1.-censored)*numpy.log(hazard) - time*hazard))

    def random(hazard, period):
        '''
        Array valued function of arrays.
            hazard - Array of floats representing hazard rates
            period - Array of floats representing observation period durations
        '''
        if len(hazard.shape) > len(period.shape):
            F = numpy.empty_like(hazard.shape)
        else:
            F = numpy.empty_like(period.shape)
        #TODO: What if a dimension has length 1?
        F = numpy.random.uniform(size = F.shape)
        time = -1. * numpy.log(1-F) / hazard
        censored = time >= period
        return [Event(t, c) for t, c in izip(time,censored)]
    
    dtype = Event
    
    result = pymc.Stochastic(logp = logp,
                             doc = 'A constant hazard survival time node for single events.',
                             name = name,
                             parents = {'hazard':hazard,'period':period},
                             random = random,
                             trace = True,
                             value = value,
                             dtype = dtype,
                             rseed = 1,
                             observed = observed,
                             cache_depth = 2,
                             plot = True,
                             verbose = 0)
    return result
