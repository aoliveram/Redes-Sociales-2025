library(reticulate)

reticulate::install_python()

reticulate::py_install(packages = "matplotlib")
reticulate::py_install(packages = "numpy")
reticulate::py_install(packages = "networkx")
reticulate::py_install(packages = "pandas")
reticulate::py_install(packages = "scipy")

py_require(c('matplotlib', 'numpy', 'networkx', 'pandas', 'scipy'))

py_run_file('Script_Auxiliar.py')
