### base of template script ###
template_script="population-dux.sh"

### list of all wilms normal samples ###
wilms=(
	"PD40713d" "PD48687b" "PD48688b" "PD48690b" "PD48691b" "PD48692b" "PD48696b"
	"PD48699b" "PD48700b" "PD48701b" "PD48705b" "PD48706b" "PD48709b" "PD48710b"
	"PD48712b" "PD48713b" "PD48714b" "PD48715b" "PD48717b" "PD48718b" "PD48719b"
	"PD48720b" "PD48721b" "PD48722b" "PD48723b" "PD48724b" "PD48726b" "PD49171b"
	"PD49173b" "PD49174b" "PD49177b" "PD49183b" "PD49184b" "PD49185b" "PD49186b"
	"PD49187b" "PD49189b" "PD49191b" "PD49193b" "PD49195b" "PD49197b" "PD49198b"
	"PD49199b" "PD49200b" "PD49204b" "PD49209b" "PD49210b" "PD49211b" "PD49213b"
	"PD49214b" "PD49215b" "PD49216b" "PD49223b" "PD49348v" "PD50589b" "PD50590b"
	"PD50593b" "PD50594f" "PD50596i" "PD50599b" "PD50600b" "PD50602b" "PD50662b"
	"PD50663b" "PD50667b" "PD50669b" "PD50675b" "PD50678b" "PD50682b" "PD50683b"
	"PD50686b" "PD50695b" "PD50696b" "PD50699b" "PD50707b" "PD50709b" "PD50715b"
	"PD50716b" "PD50717b" "PD50719b" "PD50720b" "PD50724b" "PD50726b" "PD50730b"
	"PD50733b" "PD50734b" "PD51622b" "PD52241n" "PD53640d" "PD54846m" "PD54846m")

### generate a new script using the sample ID for each replacing all instances of the PD54859b ###
for sample in "${wilms[@]}"; do
    new_script="script_${sample}.sh"

    sed -e "s/PD54859b/${sample}/g" \
        "$template_script" > "$new_script"

    echo "Generated $new_script"
done
