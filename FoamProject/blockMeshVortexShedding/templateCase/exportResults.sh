#!/bin/bash
cd "${0%/*}"

let file_name="force_coeffs_default"
let dir_name="test_case"
let export_fields=false

#### ARGUMENT HANDLING ####

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -f|--filename)
      file_name="$2"
      shift # past argument
      shift # past value
      ;;
    -d|--directory)
      dir_name="$2"
      shift # past argument
      shift # past value
      ;;
    -a|--all)
      export_fields=true
      shift # past argument
      ;;
    --default)
      DEFAULT=YES
      shift # past argument
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

###########################
# Create directory
mkdir ../../../Data/${dir_name}

checkMesh > ../../../Data/${dir_name}/mesh.txt
grep -oP 'cells:\s*\K\d+' ../../../Data/${dir_name}/mesh.txt > ../../../Data/${dir_name}/num_cells.dat
grep -oP 'Total volume =\s*\K\d+.*\d+' ../../../Data/${dir_name}/mesh.txt > ../../../Data/${dir_name}/total_volume.dat

# Export solution fields
if [ "$export_fields" = true ]; then
  cp -r ./constant ../../../Data/${dir_name}
  cp -r ./processor* ../../../Data/${dir_name}
fi

# Uncomment header and save as csv
sed -e 's/# T/T/g' ./postProcessing/forceCoeffs_object/0/coefficient_0.dat > ../../../Data/${dir_name}/${file_name}.csv

# Change dir
cd ../../../Data/${dir_name}

# Reformat csv with format: comment=#, sep="\t"
sed -i -e 's/ //g' ${file_name}.csv
