#!/bin/bash

##### ##### ##### ##### ##### ##### ##### ##### #####
# It's like nf-core's version, but jankier      #####
##### ##### ##### ##### ##### ##### ##### ##### #####

# this 
# 1. finds all the params
# 2. checks that they are in nextflow_schema.json
# 3. checks that their default value is in nextflow_schema.json
# 4. checks that every param and its default value is in the config file
# 5. checks that every process with its default values are in the config file

##### ##### ##### ##### ##### ##### ##### ##### #####
# finding all the params                        #####
##### ##### ##### ##### ##### ##### ##### ##### #####

params=$(grep param main.nf | awk '{ ( $1=$1 ) ; print $0 }' | grep ^param | awk '{print $1}' | sort | uniq | cut -f 2 -d '.')

for param in ${params[@]}
do
    echo "evaluating params.$param"

    # default value
    default="$(grep -w params.$param main.nf | grep "=" | head -n 1 | sed 's/.*=//g' | awk '{ ( $1=$1 ) ; print $0 }' ) "
    echo -e "default:\t$default"

    # 2. checks that they are in nextflow_schema.json
    schema_check=$(grep "\"$param\":" nextflow_schema.json | head -n 1)
    if [ -n "$schema_check" ]
    then
        # 3. checks that their default value is in nextflow_schema.json
        schema_default="$(grep "\"$param\":" -A 10 nextflow_schema.json | grep "}" -B 10 -m 1 | grep "default" | head -n 1 | sed 's/.*://g' | awk '{ ( $1=$1 ) ; print $0 }'  )"
        echo -e "schema:\t\t$schema_default"
       
    else
        echo "$param was not included in nextflow_schema.json file!!!"
    fi

    config_check=$(grep "params.$param =" configs/cecret_config_template.config | head -n 1)
    if [ -n "$config_check" ]
    then
        # 4. checks that every param and its default value is in the config file
        config_default="$(grep -w params.$param configs/cecret_config_template.config | grep -v "#" | head -n 1 | sed 's/.*=//g' | awk '{ ( $1=$1 ) ; print $0 }' )"
        echo -e "config:\t\t$config_default"
       
    else
        echo "$param was not included in configs/cecret_config_template.config file!!!"
    fi

    echo ""

done

echo "##### ##### ##### ##### #####"
echo "# Now for the processes #####"
echo "##### ##### ##### ##### #####"


processes=$(grep ^process -h modules/*nf | awk '{print $2}' | sort | uniq )
echo ${processes[@]}

for process in ${processes[@]}
do
    #echo "getting information for $process"
    #grep "process $process {" -h modules/*nf -A 200 | grep -e "when:" -e "input:" -B 200 -m 1 | grep -v "#UPHLICA" | grep -v "when:" | grep -v "input:" | grep -v "tag"

    echo -e "//\twithName:$process{"
    while read line
    do
        key=$(echo $line | awk '{print $1}' | grep -v process )
        if [ -n "$key" ]
        then
            value=$(echo $line | sed "s|$key||g" | awk '{ ( $1=$1 ) ; print $0 }' | sed 's/\"${params.outdir}\"/cecret/g')
            echo -e "//\t\t$key = \"$value\""
        fi
    done < <(grep "process $process {" -h modules/*nf -A 200 | grep -e "when:" -e "input:" -B 200 -m 1 | grep -v "#UPHLICA" | grep -v "when:" | grep -v "input:" | grep -v "tag")
    echo -e "//\t}"

    #echo ""
    #grep "//  withName:$process{" -A 20 configs/cecret_config_template.config | grep -m 1 "}" -B 20
done