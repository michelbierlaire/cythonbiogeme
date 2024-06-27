import setuptools
import os
import platform

def get_ext_modules():
    from Cython.Build import cythonize
    ext_modules = [
        setuptools.Extension(
            "cythonbiogeme",
            sources=[
                "src/cythonbiogeme/cpp/cythonbiogeme.pyx",
                # add other sources here
            ],
            include_dirs=["src", os.path.join(os.path.dirname(__file__), "src"), get_numpy_include()],
            language="c++",
            extra_compile_args=["-std=c++11"],
            extra_link_args=["-std=c++11"]
        )
    ]
    if platform.system() == "Windows":
        ext_modules[0].extra_compile_args.append("-DMS_WIN64")
        ext_modules[0].extra_link_args.extend([
            "-static", "-static-libstdc++", "-static-libgcc", "-lpthread", "-mms-bitfields", "-mwindows"
        ])
    return cythonize(ext_modules)

def get_numpy_include():
    import numpy
    return numpy.get_include()

setuptools.setup(
    ext_modules=get_ext_modules()
)
