import xarray

class TimeObj(object):

    def __init__(self, nc):
        self.timeVarName = ''
        self.timeVar = []
        for vName, var in nc.items():
            if getattr(var, 'standard_name', '') == 'time' or \
               getattr(var, 'long_name', '') == 'Time axis':
                self.timeName = vName
                self.timeVar = var


    def getValues(self):
        return self.timeVar[:]


    def getSize(self):
        return len(self.timeVar)


    def getTimeAsString(self, timeIndex):
        timeVal =  self.timeVar.values[timeIndex]
        year, month, day = timeVal.year, timeVal.month, timeVal.day
        return f'{year}-{month}-{day}'


###############################################################################
def test():
    nc = xarray.open_dataset('../data/nz/U.nc')
    to = TimeObj(nc)
    print(to.timeVarName)
    print(to.timeVar)
    print(f'number of steps: {to.getSize()}')
    for timeIndex in range(to.getSize()):
        print(f'index: {timeIndex} date: {to.getTimeAsString(timeIndex)}')

if __name__ == '__main__':
    test()