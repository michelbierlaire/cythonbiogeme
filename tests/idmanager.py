"""Simplified version of the Biogeme implementation for the sake of testing the Cython interface

:author: Michel Bierlaire
:date: Mon Apr 17 10:55:08 2023
"""
import logging
from collections import namedtuple
from enum import Enum, auto


class TypeOfElementaryExpression(Enum):
    """
    Defines the types of elementary expressions
    """

    VARIABLE = auto()
    BETA = auto()
    FREE_BETA = auto()
    FIXED_BETA = auto()
    RANDOM_VARIABLE = auto()
    DRAWS = auto()


ElementsTuple = namedtuple('ElementsTuple', 'expressions indices names')

logger = logging.getLogger(__name__)


class IdManager:
    """Class combining managing the ids of an arithmetic expression."""

    def __init__(self, expressions, dataframe, number_of_draws):
        """Ctor"""
        self.expressions = expressions
        self.dataframe = dataframe
        self.number_of_draws = number_of_draws
        self.elementary_expressions = None
        self.free_betas = None
        self.free_betas_values = None
        self.number_of_free_betas = 0
        self.fixed_betas = None
        self.fixed_betas_values = None
        self.bounds = None
        self.random_variables = None
        self.draws = None
        self.variables = None
        self.requires_draws = False
        for f in self.expressions:
            the_variables = f.set_of_elementary_expression(
                the_type=TypeOfElementaryExpression.VARIABLE
            )
        self.prepare()

    def changeInitValues(self, betas):
        """Modifies the values of the pameters

        :param betas: dictionary where the keys are the names of the
                      parameters, and the values are the new value for
                      the parameters.
        :type betas: dict(string:float)
        """

        def get_value(name):
            v = betas.get(name)
            if v is None:
                return self.free_betas.expressions[name].initValue
            return v

        self.free_betas_values = [get_value(x) for x in self.free_betas.names]

    def expressions_names_indices(self, dict_of_elements):
        """Assigns consecutive indices to expressions

        :param dict_of_elements: dictionary of expressions. The keys
            are the names.
        :type dict_of_elements: dict(str: biogeme.expressions.Expression)

        :return: a tuple with the original dictionary, the indices,
            and the sorted names.
        :rtype: ElementsTuple
        """
        indices = {}
        names = {}
        names = sorted(dict_of_elements)
        for i, v in enumerate(names):
            indices[v] = i

        return ElementsTuple(expressions=dict_of_elements, indices=indices, names=names)

    def prepare(self):
        """Extract from the formulas the literals (parameters,
        variables, random variables) and decide a numbering convention.

        The numbering is done in the following order:

        (i) free betas,
        (ii) fixed betas,
        (iii) random variables for numerical integration,
        (iv) random variables for Monte-Carlo integration,
        (v) variables

        The numbering convention will be performed for all expressions
        together, so that the same elementary expressions in several
        expressions will have the same index.


        """

        # Free parameters (to be estimated), sorted by alphabetical order
        expr = {}
        for f in self.expressions:
            d = f.dict_of_elementary_expression(
                the_type=TypeOfElementaryExpression.FREE_BETA
            )
            expr = dict(expr, **d)

        self.free_betas = self.expressions_names_indices(expr)

        self.bounds = [
            (
                self.free_betas.expressions[b].lb,
                self.free_betas.expressions[b].ub,
            )
            for b in self.free_betas.names
        ]
        self.number_of_free_betas = len(self.free_betas.names)
        # Fixed parameters (not to be estimated), sorted by alphatical order.
        expr = {}
        for f in self.expressions:
            d = f.dict_of_elementary_expression(
                the_type=TypeOfElementaryExpression.FIXED_BETA
            )
            expr = dict(expr, **d)
        self.fixed_betas = self.expressions_names_indices(expr)

        # Random variables for numerical integration
        expr = {}
        for f in self.expressions:
            d = f.dict_of_elementary_expression(
                the_type=TypeOfElementaryExpression.RANDOM_VARIABLE
            )
            expr = dict(expr, **d)
        self.random_variables = self.expressions_names_indices(expr)

        # Draws
        expr = {}
        for f in self.expressions:
            d = f.dict_of_elementary_expression(
                the_type=TypeOfElementaryExpression.DRAWS
            )
            expr = dict(expr, **d)
        self.draws = self.expressions_names_indices(expr)

        # Variables
        # Here, we do not extract the variables from the
        # formulas. Instead, we use all the variables in the database.
        if self.dataframe is not None:
            variables_names = list(self.dataframe.columns.values)
            variables_indices = {}
            for i, v in enumerate(variables_names):
                variables_indices[v] = i
            self.variables = ElementsTuple(
                expressions=None,
                indices=variables_indices,
                names=variables_names,
            )
        else:
            self.variables = ElementsTuple(expressions=None, indices=None, names=[])

        # Merge all the names
        elementary_expressions_names = (
            self.free_betas.names
            + self.fixed_betas.names
            + self.random_variables.names
            + self.draws.names
            + self.variables.names
        )

        elementary_expressions_indices = {
            v: i for i, v in enumerate(elementary_expressions_names)
        }

        self.elementary_expressions = ElementsTuple(
            expressions=None,
            indices=elementary_expressions_indices,
            names=elementary_expressions_names,
        )

        self.free_betas_values = [
            self.free_betas.expressions[x].initValue for x in self.free_betas.names
        ]
        self.fixed_betas_values = [
            self.fixed_betas.expressions[x].initValue for x in self.fixed_betas.names
        ]

    def setData(self, sample):
        """Specify the sample

        :param sample: map of the panel data (see
            :func:`biogeme.database.Database.buildPanelMap`)
        :type sample: pandas.DataFrame

        """
        for f in self.expressions:
            f.cpp.setData(sample)
