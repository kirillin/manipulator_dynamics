"""
    The module contains stuff for make elements of Xi regressor
"""
from libs.initialization import n, nL


args_list = lambda a, b: ['{}={}'.format(a[i], b[i]) for i in range(len(a))]


class Class:

    def __init__(self, title, constructor, imports=None):
        self.__class = ''
        if imports:
            self.__class = str(imports) + '\n'
        self.__class += 'class {0}:\n\n'.format(title)
        self.addMethod(constructor)

    def addMethod(self, method):
        self.__class += str(method) + '\n'

    def addMethods(self, methods):
        for method in methods:
            self.__class += str(method) + '\n'

    def __str__(self):
        return self.__class


class Method:

    def __init__(self, title, body=['pass'], args=None, args_inits=None, tabs='\t'):
        self.__title = title
        self.__body = body
        self.__args = args
        self.__inits_args = args_inits
        self.__tabs = tabs

    def myBody(self):
        body = ''
        for row in self.__body:
            body += self.__tabs + '\t' + row + '\n'
        return body

    def __str__(self):
        if self.__args and len(self.__args) > 0:
            if self.__inits_args:
                args_string = ', '.join(args_list(self.__args, self.__inits_args))
            else:
                args_string = ', '.join(self.__args)
            __method = '{0}def {1}(self, {2}):\n'.format(self.__tabs, self.__title, args_string)
        else:
            __method = '{0}def {1}(self):\n'.format(self.__tabs, self.__title)
        __method += self.myBody()
        return __method


class RegressorElement:

    IMPORTS = 'from numpy import cos, sin, sqrt, tan, zeros, array\n\n'

    CONSTRUCTOR_ARGS = ['q', 'dq', 'ddq', 'a', 'd', 'delta']
    CONSTRUCTOR_ARGS_INITS = [tuple(0 for i in range(n)) for arg in range(len(CONSTRUCTOR_ARGS))]
    CONSTRUCTOR_BODY = ['self.{0} = {0}'.format(arg) for arg in CONSTRUCTOR_ARGS]

    OPL_X_METHOD_ARGS = ['q', 'dq', 'ddq', 'a', 'd', 'theta']
    OPL_X_METHOD_BODY = ['opL_{0} = {1}',
                         'return opL_{0}']

    GET_XI_METHOD_ARGS = ['q', 'dq', 'ddq']
    GET_XI_METHOD_BODY = ['XI = zeros({0})'.format(nL)] + \
                         ['theta = array(self.delta) - array(q)'] + \
                         ['XI[{0}] = self.opL{0}(q, dq, ddq, self.a, self.d, theta)'.format(i) for i in range(nL)] + \
                         ['return XI']

    def __init__(self, j, i):
        constructor = Method('__init__', args=RegressorElement.CONSTRUCTOR_ARGS,
                             args_inits=RegressorElement.CONSTRUCTOR_ARGS_INITS,
                             body=RegressorElement.CONSTRUCTOR_BODY)
        self.regressorElementClass = Class('Xi{0}{1}'.format(j, i), constructor,
                                           imports=RegressorElement.IMPORTS)

    def addOpL(self, k, expr):
        opL = Method('opL{0}', body=RegressorElement.OPL_X_METHOD_BODY,
                     args=RegressorElement.OPL_X_METHOD_ARGS)
        self.regressorElementClass.addMethod(str(opL).format(k, expr))

    def addGetter(self):
        getter = Method('getXi', body=RegressorElement.GET_XI_METHOD_BODY,
                        args=RegressorElement.GET_XI_METHOD_ARGS)
        self.regressorElementClass.addMethod(getter)

    def __str__(self):
        self.addGetter()
        return str(self.regressorElementClass)


class Regressor:

    IMPORTS = 'import numpy as np\n'

    CONSTRUCTOR_ARGS = ['a', 'd', 'delta']
    CONSTRUCTOR_BODY = []

    GETTER_ARGS = ['q', 'dq', 'ddq']
    GETTER_BODY = [
        'n = len(q)',
        'xi = np.empty(0)',
        'for i in range(n):',
        '\tfor j in range(n):',
        '\t\tif self._xi[i][j] != 0:',
        '\t\t\txi_ij = self._xi[i][j].getXi(q, dq, ddq, theta)',
        '\t\t\txi = np.concatenate((xi, xi_ij))',
        '\t\telse:',
        '\t\t\txi = np.concatenate((xi, np.zeros({0})))'.format(nL),
        'return xi.reshape(n, n * {0})'.format(nL)
    ]

    GET_WELL_COLS_BODY = [
        'all = set(range({0}))'.format(n * nL),
        'well = all.difference(self._linerCols)',
        'return list(well)'
    ]

    def __init__(self, path='xi/', zero_columns=[]):
        self.init(path, zero_columns)     # !!! auto-generator

        constructor = Method('__init__', Regressor.CONSTRUCTOR_BODY, args=Regressor.CONSTRUCTOR_ARGS)
        self.regressorClass = Class('Xi', constructor=constructor, imports=Regressor.IMPORTS)

        # getter = Method('getXi', args=['q', 'dq', 'ddq'], body=Regressor.GETTER_BODY)
        # self.regressorClass.addMethod(getter)

        getWellColNums = Method('getWellColNums', body=Regressor.GET_WELL_COLS_BODY)
        self.regressorClass.addMethod(getWellColNums)

    def init(self, path, zero_columns):
        # generate imports
        for i in range(n):
            for j in range(i, n):
                import_element = 'from {0}xi_{1}{2} import Xi{1}{2}\n'.format(path.replace('/', '.'), i, j)
                Regressor.IMPORTS += import_element

        # generate constructor body
        for i in range(n):
            args_string = ', '.join(args_list(Regressor.CONSTRUCTOR_ARGS, Regressor.CONSTRUCTOR_ARGS))
            row_list = ['\t\t\t\t0' for j in range(0, i)] + ['Xi{0}{1}({2})'.format(i, j, args_string) for j in range(i, n)]
            row = ', '.join(row_list)
            Regressor.CONSTRUCTOR_BODY.append('xi{0} = [{1}]'.format(i, row))
        Regressor.CONSTRUCTOR_BODY.append('self._xi = [{0}]'.format(', '.join(['xi{0}'.format(i) for i in range(n)])))
        Regressor.CONSTRUCTOR_BODY.append('self._linerCols = [{0}]'.format(', '.join(zero_columns)))

    def __str__(self):
        return str(self.regressorClass)
