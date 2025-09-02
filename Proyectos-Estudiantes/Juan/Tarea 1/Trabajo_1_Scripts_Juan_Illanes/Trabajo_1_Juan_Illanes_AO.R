library(reticulate)

use_python("/Library/Frameworks/Python.framework/Versions/3.13/bin/python3", required = TRUE)

sapply(c("numpy","pandas","matplotlib","scipy","networkx"), py_module_available)
py_config()  # para verificar el intérprete y la arquitectura

reticulate::py_install(packages = "matplotlib")
reticulate::py_install(packages = "numpy")
reticulate::py_install(packages = "networkx")
reticulate::py_install(packages = "pandas")
reticulate::py_install(packages = "scipy")

py_require(c('matplotlib', 'numpy', 'networkx', 'pandas', 'scipy'))

py_run_file('Script_Auxiliar.py')
