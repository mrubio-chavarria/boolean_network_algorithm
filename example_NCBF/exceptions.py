
class NoSolutionException(Exception):
    """
    DESCRIPTION:
    An exception to handle where there are no terms to choose when
    performing inference. This exception is launched when there is
    no solution.
    """

    # Methods
    def __init__(self):
        """
        DESCRIPTION:
        Constructor of the class
        """
        super().__init__('No solution found')


class InputModificationException(Exception):
    """
    DESCRIPTION:
    An exception to handle when the conflicts solving procedure 
    demands to modify the input.
    """

    # Methods
    def __init__(self):
        """
        DESCRIPTION:
        Constructor of the class
        """
        super().__init__('Input modification detected')