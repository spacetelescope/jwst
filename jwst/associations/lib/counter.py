class Counter():
    """Like itertools.count but access to the current value"""
    def __init__(self, start=0, step=1, end=None):
        self.value = start
        self.step = step
        self.end = end

    def __iter__(self):
        return self

    def __next__(self):
        if self.end is not None and \
           abs(self.value) > abs(self.end):
            raise StopIteration
        self.value += self.step
        return self.value

    def next(self):
        """python2 compatibility"""
        return self.__next__()

    def set(self, value):
        """Set new value for counter"""
        self.value = value
        return self.value
