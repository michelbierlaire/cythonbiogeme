""" Minimal version of biogeme expressions, just for testing purposes

:author: Michel Bierlaire

:date: Mon Apr 17 08:55:21 2023
"""

import logging
from itertools import chain

logger = logging.getLogger(__name__)

from idmanager import TypeOfElementaryExpression

def isNumeric(obj):
    """Checks if an object is numeric
    :param obj: obj to be checked
    :type obj: Object

    """
    return isinstance(obj, (int, float, bool))


class Expression:
    """This is the general arithmetic expression in biogeme.
    It serves as a base class for concrete expressions.
    """

    def __init__(self):
        """Constructor"""

        self.children = []  #: List of children expressions

        self.id_manager = None  #: in charge of the IDs
        self.keep_id_manager = None  #: a copy of the ID manager

        self.fixedBetaValues = None
        """values of the Beta that are not estimated
        """

        self.numberOfDraws = None
        """number of draws for Monte Carlo integration
        """

        self._row = None
        """Row of the database where the values of the variables are found
        """

        self.missingData = 99999
        """ Value interpreted as missing data
        """

        self.locked = False
        """ Meaningful only for multiple expressions (Catalogs). If
        True, it is not possible to change the specification.
        """

    def setIdManager(self, id_manager):
        """The ID manager contains the IDs of the elementary expressions.

        It is externally created, as it may nee to coordinate the
        numbering of several expressions. It is stored only in the
        expressions of type Elementary.

        :param id_manager: ID manager to be propagated to the
            elementary expressions. If None, all the IDs are set to None.
        :type id_manager: class IdManager
        """
        self.id_manager = id_manager
        for e in self.get_children():
            e.setIdManager(id_manager)

    def set_of_elementary_expression(self, the_type):
        """Extract a dict with all elementary expressions of a specific type

        :param the_type: the type of expression
        :type  the_type: TypeOfElementaryExpression

        :return: returns a set with the names of the elementary expressions
        :rtype: set(string.Expression)

        """
        return set(self.dict_of_elementary_expression(the_type).keys())

    def dict_of_elementary_expression(self, the_type):
        """Extract a dict with all elementary expressions of a specific type

        :param the_type: the type of expression
        :type  the_type: TypeOfElementaryExpression

        :return: returns a dict with the variables appearing in the
               expression the keys being their names.
        :rtype: dict(string:biogeme.expressions.Expression)

        """
        return dict(
            chain(
                *(
                    e.dict_of_elementary_expression(the_type).items()
                    for e in self.children
                )
            )
        )

            
    def __add__(self, other):
        """
        Operator overloading. Generate an expression for addition.

        :param other: expression to be added
        :type other: biogeme.expressions.Expression

        :return: self + other
        :rtype: biogeme.expressions.Expression

        :raise BiogemeError: if one of the expressions is invalid, that is
            neither a numeric value or a
            biogeme.expressions.Expression object.
        """
        return Plus(self, other)

    def __radd__(self, other):
        """
        Operator overloading. Generate an expression for addition.

        :param other: expression to be added
        :type other: biogeme.expressions.Expression

        :return: other + self
        :rtype: biogeme.expressions.Expression

        :raise BiogemeError: if one of the expressions is invalid, that is
            neither a numeric value or a
            biogeme.expressions.Expression object.

        """
        return Plus(other, self)

    def getClassName(self):
        """
        Obtain the name of the top class of the expression structure

        :return: the name of the class
        :rtype: string
        """
        n = type(self).__name__
        return n

    def getSignature(self):
        """The signature of a string characterizing an expression.

        This is designed to be communicated to C++, so that the
        expression can be reconstructed in this environment.

        The list contains the following elements:

            1. the signatures of all the children expressions,
            2. the name of the expression between < >
            3. the id of the expression between { }
            4. the number of children between ( )
            5. the ids of each children, preceeded by a comma.

        Consider the following expression:

        .. math:: 2 \\beta_1  V_1 -
            \\frac{\\exp(-\\beta_2 V_2) }
            { \\beta_3  (\\beta_2 \\geq \\beta_1)}.

        It is defined as::

            2 * beta1 * Variable1 - expressions.exp(-beta2*Variable2) /
                 (beta3 * (beta2 >= beta1))

        And its signature is::

            [b'<Numeric>{4780527008},2',
             b'<Beta>{4780277152}"beta1"[0],0,0',
             b'<Times>{4780526952}(2),4780527008,4780277152',
             b'<Variable>{4511837152}"Variable1",5,2',
             b'<Times>{4780527064}(2),4780526952,4511837152',
             b'<Beta>{4780277656}"beta2"[0],1,1',
             b'<UnaryMinus>{4780527120}(1),4780277656',
             b'<Variable>{4511837712}"Variable2",6,3',
             b'<Times>{4780527176}(2),4780527120,4511837712',
             b'<exp>{4780527232}(1),4780527176',
             b'<Beta>{4780277264}"beta3"[1],2,0',
             b'<Beta>{4780277656}"beta2"[0],1,1',
             b'<Beta>{4780277152}"beta1"[0],0,0',
             b'<GreaterOrEqual>{4780527288}(2),4780277656,4780277152',
             b'<Times>{4780527344}(2),4780277264,4780527288',
             b'<Divide>{4780527400}(2),4780527232,4780527344',
             b'<Minus>{4780527456}(2),4780527064,4780527400']

        :return: list of the signatures of an expression and its children.
        :rtype: list(string)

        """
        listOfSignatures = []
        for e in self.get_children():
            listOfSignatures += e.getSignature()
        mysignature = f'<{self.getClassName()}>'
        mysignature += f'{{{self.get_id()}}}'
        mysignature += f'({len(self.get_children())})'
        for e in self.get_children():
            mysignature += f',{e.get_id()}'
        listOfSignatures += [mysignature.encode()]
        return listOfSignatures

    def get_id(self):
        """Retrieve the id of the expression used in the signature

        :return: id of the object
        :rtype: int
        """
        return id(self)

    def get_children(self):
        """Retrieve the list of children

        :return: list of children
        :rtype: list(Expression)
        """
        return self.children


class BinaryOperator(Expression):
    """
    Base class for arithmetic expressions that are binary operators.
    This expression is the result of the combination of two expressions,
    typically addition, substraction, multiplication or division.
    """

    def __init__(self, left, right):
        """Constructor

        :param left: first arithmetic expression
        :type left: biogeme.expressions.Expression

        :param right: second arithmetic expression
        :type right: biogeme.expressions.Expression

        :raise BiogemeError: if one of the expressions is invalid, that is
            neither a numeric value or a
            biogeme.expressions.Expression object.

        """
        Expression.__init__(self)
        self.left = left
        self.right = right
        self.children.append(self.left)
        self.children.append(self.right)


class Plus(BinaryOperator):
    """
    Addition expression
    """

    def __init__(self, left, right):
        """Constructor

        :param left: first arithmetic expression
        :type left: biogeme.expressions.Expression

        :param right: second arithmetic expression
        :type right: biogeme.expressions.Expression
        """
        BinaryOperator.__init__(self, left, right)


class Elementary(Expression):
    """Elementary expression.

    It is typically defined by a name appearing in an expression. It
    can be a variable (from the database), or a parameter (fixed or to
    be estimated using maximum likelihood), a random variable for
    numrerical integration, or Monte-Carlo integration.

    """

    def __init__(self, name):
        """Constructor

        :param name: name of the elementary experession.
        :type name: string

        """
        Expression.__init__(self)
        self.name = name  #: name of the elementary expressiom

        self.elementaryIndex = None
        """The id should be unique for all elementary expressions
        appearing in a given set of formulas.
        """


class Variable(Elementary):
    """Explanatory variable

    This represents the explanatory variables of the choice
    model. Typically, they come from the data set.
    """

    def __init__(self, name):
        """Constructor

        :param name: name of the variable.
        :type name: string
        """
        Elementary.__init__(self, name)
        # Index of the variable
        self.variableId = None

    def setIdManager(self, id_manager=None):
        """The ID manager contains the IDs of the elementary expressions.

        It is externally created, as it may need to coordinate the
        numbering of several expressions. It is stored only in the
        expressions of type Elementary.

        :param id_manager: ID manager to be propagated to the
            elementary expressions. If None, all the IDs are set to None.
        :type id_manager: class IdManager
        """

        self.id_manager = id_manager
        if id_manager is None:
            self.elementaryIndex = None
            self.variableId = None
            return
        self.elementaryIndex = self.id_manager.elementary_expressions.indices[self.name]
        self.variableId = self.id_manager.variables.indices[self.name]

    def getSignature(self):
        """The signature of a string characterizing an expression.

        This is designed to be communicated to C++, so that the
        expression can be reconstructed in this environment.

        The list contains the following elements:

            1. the name of the expression between < >
            2. the id of the expression between { }
            3. the name of the variable,
            4. the unique ID, preceeded by a comma.
            5. the variabvle ID, preceeded by a comma.

        Consider the following expression:

        .. math:: 2 \\beta_1  V_1 -
         \\frac{\\exp(-\\beta_2 V_2) }{ \\beta_3  (\\beta_2 \\geq \\beta_1)}.

        It is defined as::

            2 * beta1 * Variable1 - expressions.exp(-beta2*Variable2) /
                (beta3 * (beta2 >= beta1))

        And its signature is::

            [b'<Numeric>{4780527008},2',
             b'<Beta>{4780277152}"beta1"[0],0,0',
             b'<Times>{4780526952}(2),4780527008,4780277152',
             b'<Variable>{4511837152}"Variable1",5,2',
             b'<Times>{4780527064}(2),4780526952,4511837152',
             b'<Beta>{4780277656}"beta2"[0],1,1',
             b'<UnaryMinus>{4780527120}(1),4780277656',
             b'<Variable>{4511837712}"Variable2",6,3',
             b'<Times>{4780527176}(2),4780527120,4511837712',
             b'<exp>{4780527232}(1),4780527176',
             b'<Beta>{4780277264}"beta3"[1],2,0',
             b'<Beta>{4780277656}"beta2"[0],1,1',
             b'<Beta>{4780277152}"beta1"[0],0,0',
             b'<GreaterOrEqual>{4780527288}(2),4780277656,4780277152',
             b'<Times>{4780527344}(2),4780277264,4780527288',
             b'<Divide>{4780527400}(2),4780527232,4780527344',
             b'<Minus>{4780527456}(2),4780527064,4780527400']

        :return: list of the signatures of an expression and its children.
        :rtype: list(string)

        :raise biogeme.exceptions.BiogemeError: if no id has been defined for
            elementary expression
        :raise biogeme.exceptions.BiogemeError: if no id has been defined for
            variable
        """
        signature = f'<{self.getClassName()}>'
        signature += f'{{{self.get_id()}}}'
        signature += f'"{self.name}",{self.elementaryIndex},{self.variableId}'
        return [signature.encode()]

    def dict_of_elementary_expression(self, the_type):
        """Extract a dict with all elementary expressions of a specific type

        :param the_type: the type of expression
        :type  the_type: TypeOfElementaryExpression
        """
        if the_type == TypeOfElementaryExpression.VARIABLE:
            return {self.name: self}
        return {}
    
