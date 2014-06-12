class ConvergenceError(Exception):
    def __init__(self, value=None):
        self.value = value
    def __str(self):
        return repr(self.value)

