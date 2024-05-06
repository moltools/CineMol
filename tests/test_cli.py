# -*- coding: utf-8 -*-

"""Contains unit tests for the cinemol.cli module."""

import os
import unittest
from pathlib import Path

from click.testing import CliRunner

from cinemol.cli import main

FIXTURES = Path(__file__).parent / "fixtures"


class TestCli(unittest.TestCase):
    """Test the Command Line Interface."""

    def test_setup(self):
        """Test the setup."""
        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        self.assertEqual(result.exit_code, 0)
        self.assertTrue(result.output.startswith("Usage: main [OPTIONS] INPUT OUTPUT"))

    def test_main_with_default_options(self):
        """Test the main function with default options."""
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                os.path.join(FIXTURES, "test_input_001.sdf"),
                os.path.join(FIXTURES, "test_output_001.svg"),
            ],
        )
        self.assertEqual(result.exit_code, 0)

    def test_main_throws_error_without_specifying_output_path(self):
        """Test the main function throws error without specifying output path."""
        runner = CliRunner()
        result = runner.invoke(main, [os.path.join(FIXTURES, "test_input_001.sdf")])
        self.assertEqual(result.exit_code, 2)

    def test_main_throws_error_when_input_path_does_not_exist(self):
        """Test the main function throws error when input path does not exist."""
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                os.path.join(FIXTURES, "non_existent.sdf"),
                os.path.join(FIXTURES, "test_output_001.svg"),
            ],
        )
        self.assertEqual(result.exit_code, 2)
