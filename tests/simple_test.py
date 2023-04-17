"""
Simple test for cythonbiogeme

:author: Michel Bierlaire
:data: 

"""

import unittest
import pandas as pd
import cythonbiogeme.cythonbiogeme as cc
from idmanager import IdManager
from expressions import Variable

df = pd.DataFrame(
    {
        'Person': [1, 1, 1, 2, 2],
        'Exclude': [0, 0, 1, 0, 1],
        'Variable1': [10, 20, 30, 40, 50],
        'Variable2': [100, 200, 300, 400, 500],
        'Choice': [2, 2, 3, 1, 2],
        'Av1': [0, 1, 1, 1, 1],
        'Av2': [1, 1, 1, 1, 1],
        'Av3': [0, 1, 1, 1, 1],
    }
)


class test_cythonbiogeme(unittest.TestCase):
    def test_simple(self):
        variable_1 = Variable('Variable1')
        variable_2 = Variable('Variable2')
        the_sum = variable_1 + variable_2
        the_sum.setIdManager(None)
        id_manager = IdManager([the_sum], df, 0)
        the_sum.setIdManager(id_manager)
        variable_1.setIdManager(id_manager)
        variable_2.setIdManager(id_manager)

        cpp = cc.pyEvaluateOneExpression()
        cpp.setData(df)
        cpp.setExpression(the_sum.getSignature())
        cpp.setFreeBetas(the_sum.id_manager.free_betas_values)
        cpp.setFixedBetas(the_sum.id_manager.fixed_betas_values)
        cpp.setMissingData(the_sum.missingData)
        cpp.calculate(
            gradient=False,
            hessian=False,
            bhhh=False,
            aggregation=False,
        )
        correct_results = [110, 220, 330, 440, 550]
        f, _, _, _ = cpp.getResults()
        self.assertListEqual(list(f), correct_results)


if __name__ == '__main__':
    unittest.main()
