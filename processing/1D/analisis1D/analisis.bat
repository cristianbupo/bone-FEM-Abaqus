rmdir /S /Q resultados
del result_analisis.vtu
set user_file=..\..\master\general2DElasticDiffusion.for
abaqus job=analisis user=%user_file% ask_delete=OFF
