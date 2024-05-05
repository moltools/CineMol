# -*- coding: utf-8 -*-

"""Contains unit tests for the cinemol.version module."""

import unittest

from cinemol.version import get_git_hash, get_version


class TestGitHash(unittest.TestCase):
    """Trivially test a git hash."""

    def test_git_hash_type(self):
        """Test the git hash is a string.

        This is only meant to be an example test.
        """
        git_hash = get_git_hash()
        self.assertIsInstance(git_hash, str)


class TestVersion(unittest.TestCase):
    """Trivially test a version."""

    def test_version_type(self):
        """Test the version is a string.

        This is only meant to be an example test.
        """
        version = get_version()
        self.assertIsInstance(version, str)
