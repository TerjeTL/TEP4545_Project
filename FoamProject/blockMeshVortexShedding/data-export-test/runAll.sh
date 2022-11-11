#!/bin/bash
cd "${0%/*}"

#### CONFIG ####

num_cpu=6
clean_run=true

################


#### DEFS ######

no_prompt=false

################


#### ARGUMENT HANDLING ####

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -y|-Y)
      no_prompt=true
      shift # past argument (no value)
      ;;
    -p|--processors)
      num_cpu=$2
      shift # past argument
      shift # past value
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

# Copy mesh file to run directory (nie - using blockmesh here)
#cp ./../../../../Meshing/flat_plate_in_channel.msh .

if [ "$clean_run" = true ]; then
  # Remove previous results
  openfoam2206 foamListTimes -rm

  # Generate mesh
  openfoam2206 blockMesh

  # Check mesh
  openfoam2206 checkMesh
  if [ "$no_prompt" = false ]; then

    read -r -p "Continue? [y/N] " response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
    then
      continue
    else
      exit 0
    fi
  fi
  
  # Split mesh for parallelization
  let num_cpu_x=num_cpu/2
  sed -i "s/^numberOfSubdomains.*$/numberOfSubdomains\t${num_cpu};/" ./system/decomposeParDict
  sed -i "s/^\tn.*$/\tn (${num_cpu_x} 2 1);/" ./system/decomposeParDict

  openfoam2206 decomposePar -force
fi


# Run icoFoam in parallel
openfoam2206 mpirun --allow-run-as-root -np "$num_cpu" icoFoam -parallel
