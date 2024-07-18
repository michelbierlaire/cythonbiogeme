import os

import setuptools
from Cython.Build import cythonize
import numpy
import platform


def get_ext_modules():
    ext_modules = [
        setuptools.Extension(
            "cythonbiogeme.cythonbiogeme",
            sources=[
                "src/cythonbiogeme/cpp/bioFormula.cc",
                "src/cythonbiogeme/cpp/bioMemoryManagement.cc",
                "src/cythonbiogeme/cpp/biogeme.cc",
                "src/cythonbiogeme/cpp/bioDerivatives.cc",
                "src/cythonbiogeme/cpp/bioExceptions.cc",
                "src/cythonbiogeme/cpp/bioExprAnd.cc",
                "src/cythonbiogeme/cpp/bioExprBelongsTo.cc",
                "src/cythonbiogeme/cpp/bioExprConditionalSum.cc",
                "src/cythonbiogeme/cpp/bioExprCos.cc",
                "src/cythonbiogeme/cpp/bioExprDerive.cc",
                "src/cythonbiogeme/cpp/bioExprDivide.cc",
                "src/cythonbiogeme/cpp/bioExprDraws.cc",
                "src/cythonbiogeme/cpp/bioExprElem.cc",
                "src/cythonbiogeme/cpp/bioExprEqual.cc",
                "src/cythonbiogeme/cpp/bioExprExp.cc",
                "src/cythonbiogeme/cpp/bioExprFixedParameter.cc",
                "src/cythonbiogeme/cpp/bioExprFreeParameter.cc",
                "src/cythonbiogeme/cpp/bioExprGaussHermite.cc",
                "src/cythonbiogeme/cpp/bioExprGreater.cc",
                "src/cythonbiogeme/cpp/bioExprGreaterOrEqual.cc",
                "src/cythonbiogeme/cpp/bioExprIntegrate.cc",
                "src/cythonbiogeme/cpp/bioExprLess.cc",
                "src/cythonbiogeme/cpp/bioExprLessOrEqual.cc",
                "src/cythonbiogeme/cpp/bioExprLinearUtility.cc",
                "src/cythonbiogeme/cpp/bioExprLiteral.cc",
                "src/cythonbiogeme/cpp/bioExprLog.cc",
                "src/cythonbiogeme/cpp/bioExprLogLogit.cc",
                "src/cythonbiogeme/cpp/bioExprLogLogitFullChoiceSet.cc",
                "src/cythonbiogeme/cpp/bioExprLogzero.cc",
                "src/cythonbiogeme/cpp/bioExprMax.cc",
                "src/cythonbiogeme/cpp/bioExprMin.cc",
                "src/cythonbiogeme/cpp/bioExprMinus.cc",
                "src/cythonbiogeme/cpp/bioExprMontecarlo.cc",
                "src/cythonbiogeme/cpp/bioExprMultSum.cc",
                "src/cythonbiogeme/cpp/bioExprNormalCdf.cc",
                "src/cythonbiogeme/cpp/bioExprNotEqual.cc",
                "src/cythonbiogeme/cpp/bioExprNumeric.cc",
                "src/cythonbiogeme/cpp/bioExprOr.cc",
                "src/cythonbiogeme/cpp/bioExprPanelTrajectory.cc",
                "src/cythonbiogeme/cpp/bioExprPlus.cc",
                "src/cythonbiogeme/cpp/bioExprPower.cc",
                "src/cythonbiogeme/cpp/bioExprPowerConstant.cc",
                "src/cythonbiogeme/cpp/bioExprRandomVariable.cc",
                "src/cythonbiogeme/cpp/bioExprSin.cc",
                "src/cythonbiogeme/cpp/bioExprTimes.cc",
                "src/cythonbiogeme/cpp/bioExprUnaryMinus.cc",
                "src/cythonbiogeme/cpp/bioExprVariable.cc",
                "src/cythonbiogeme/cpp/bioExpression.cc",
                "src/cythonbiogeme/cpp/bioGaussHermite.cc",
                "src/cythonbiogeme/cpp/bioGhFunction.cc",
                "src/cythonbiogeme/cpp/bioNormalCdf.cc",
                "src/cythonbiogeme/cpp/bioSeveralExpressions.cc",
                "src/cythonbiogeme/cpp/bioSeveralFormulas.cc",
                "src/cythonbiogeme/cpp/bioString.cc",
                "src/cythonbiogeme/cpp/bioThreadMemory.cc",
                "src/cythonbiogeme/cpp/bioThreadMemoryOneExpression.cc",
                "src/cythonbiogeme/cpp/bioThreadMemorySimul.cc",
                "src/cythonbiogeme/cpp/bioVectorOfDerivatives.cc",
                "src/cythonbiogeme/cpp/cythonbiogeme.pyx",
                "src/cythonbiogeme/cpp/evaluateExpressions.cc",
                "src/cythonbiogeme/cpp/validity_check.cc",
            ],
            include_dirs=["src", numpy.get_include()],
            language="c++",
            extra_compile_args=["-std=c++14"],
            extra_link_args=["-std=c++14"],
        )
    ]

    if platform.system() == "Windows":
        ext_modules[0].extra_compile_args.append("/std:c++14")
        #ext_modules[0].extra_link_args.extend(
        #    [
        #        "-static",
        #        "-static-libstdc++",
        #        "-static-libgcc",
        #        "-lpthread",
        #        "-mms-bitfields",
        #        "-mwindows",
        #    ]
        #)

    return cythonize(ext_modules)


def get_numpy_include():
    import numpy

    return numpy.get_include()


setuptools.setup(
    packages=setuptools.find_packages(where="src"),
    package_dir={"": "src"},
    ext_modules=get_ext_modules(),
    include_package_data=True,
)
