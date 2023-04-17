"""
Simple test for cythonbiogeme

:author: Michel Bierlaire
:data: 

"""

import unittest

import cythonbiogeme.cythonbiogeme as cc
from idmanager import IdManager
from test_data import get_data
from expressions import Variable


class test_cythonbiogeme(unittest.TestCase):
    def simple_test(self):
        variable_1 = Variable('Variable1')
        variable_2 = Variable('Variable2')
        the_sum = variable_1 + variable_2
        the_sum.setIdManager(None)
        the_data = get_data(2)
        id_manager = IdManager([the_sum], the_data, 0)
        the_sum.setIdManager(id_manager)
        variable_1.setIdManager(id_manager)
        variable_2.setIdManager(id_manager)

        cpp = cc.pyEvaluateOneExpression()
        cpp.setData(the_data.data)
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
        self.assertListEqual(f, correct_results)
    
if __name__ == '__main__':
    unittest.main()
