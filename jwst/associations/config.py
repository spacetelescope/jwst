"""Global configuration parameters."""

DEBUG = False
"""
Debug mode.

In conjunction with logging, this mode changes the behavior
of some functions. Common use case is where associations would be
removed, they are no longer removed. This assists in determining
why they were an issue.
"""

__all__ = ["DEBUG"]
