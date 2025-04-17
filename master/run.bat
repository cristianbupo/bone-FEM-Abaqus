rmdir /S /Q resultados
set user_file=..\..\master\general2DElasticDiffusion.for
abaqus job={} user=%user_file% ask_delete=OFF
