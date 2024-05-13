import unittest as ut
from tools import verify_type
from typing import Any, Type


class TestCase(ut.TestCase):
    def create_valid(self, ty: Type, value: Any) -> Type:
        self.assertTrue(verify_type(value, ty), "Value of wrong type, expected: " + str(ty))
        return value
