# Base template
template_script="population-dux.sh"

# List of sarcoma_normal samples
sarcoma_normals=(
    "PD41847b" "PD41848b" "PD42181b" "PD46692b" "PD46696b" "PD46697b"
    "PD46698b" "PD47706b" "PD47708b" "PD51353b" "PD54858b" "PD54859b"
)

# Generate a new script per sample
for sample in "${sarcoma_normals[@]}"; do
    new_script="script_${sample}.sh"

    sed -e "s/PD54859b/${sample}/g" \
        "$template_script" > "$new_script"

    echo "Generated $new_script"
done
