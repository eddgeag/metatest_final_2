mkdir -p logs

bash -c '
set -e

echo "===== ESCENARIO 1 ====="
./run_pipeline_from_template.sh plantilla_FLI.tsv

echo "===== ESCENARIO 2 ====="
sleep 30
./run_pipeline_from_template.sh plantilla_sin_FLI.tsv

echo "===== ESCENARIO 3 ====="
sleep 30
./run_pipeline_from_template.sh plantilla_pipeline_3_analisis.tsv

echo "===== TODO TERMINADO ====="
' 2>&1 | tee logs/run_3_escenarios.log
