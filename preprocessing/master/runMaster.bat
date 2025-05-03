rmdir /S /Q resultados
set project_path=C:\Users\crist\git\bone-FEM-Abaqus\processing
set version=mechanicsFirstStep
set user_file=general2DElasticDiffusion.for
abaqus job={} user=%project_path%\%version%\%user_file% ask_delete=OFF