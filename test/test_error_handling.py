from unittest import main, TestCase
from cinnamol import error_handler, Fail, Success 


class TestErrorHandling(TestCase):
    def test_assert_fail_simple_function(self):
        """
        Check if failing a simple function results in Fail object result.
        """
        @error_handler
        def func(x): return x / 2
        result = func(None)
        self.assertTrue(isinstance(result, Fail))

    def test_assert_success_simple_function(self):
        """
        Check if succeeding a simple function results in Success object result.
        """
        @error_handler
        def func(x): return x / 2
        result = func(2)
        self.assertTrue(isinstance(result, Success))


if __name__  == "__main__":
    main()
    