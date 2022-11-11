#!/bin/bash
cd "${0%/*}"

no_prompt=false

if [ "$1" = "-y" ]; then
  no_prompt=true
fi

if [ "$no_prompt" = false ]; then
  read -r -p "Remove previous time data? [y/N] " response
  if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
    then
      # Remove previous results
      foamListTimes -rm
      rm log.foamTimes
    else
      exit 0
  fi
fi


if [ "$no_prompt" = false ]; then
  read -r -p "Remove postprocessing results? [y/N]" response
  if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
    then
      # Remove previous results
      rm -r ./postProcessing
    else
      exit 0
  fi
fi


if [ "$no_prompt" = false ]; then
  read -r -p "Remove processor directories? [y/N]" response
  if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
    then
      # Remove previous results
      rm -r ./processor*
    else
      exit 0
  fi
fi
