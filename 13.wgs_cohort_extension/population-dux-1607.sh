# Base template
template_script="population-dux.sh"

# List of orphan samples
orphan_tumours=(
    'PD46693b' 'PD36158b' 'PD36166b' 'PD34956b' 'PD36157b' 'PD36161b' 'PD34278b'
    'PD31013b' 'PD36156b' 'PD31012b' 'PD34279b' 'PD36160b' 'PD34957b' 'PD34958b'
    'PD36159b' 'PD34952b' 'PD34953b' 'PD36167b' 'PD34954b' 'PD34276b' 'PD34955b'
    'PD34956b' 'PD36160b' 'PD42727b'
)

# Generate a new script per sample
for sample in "${orphan_tumours[@]}"; do
    new_script="script_${sample}.sh"

    sed -e "s/PD54859b/${sample}/g" \
        "$template_script" > "$new_script"

    echo "Generated $new_script"
done
