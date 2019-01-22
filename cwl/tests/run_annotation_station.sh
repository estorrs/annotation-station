#!/bin/bash

CWL="cwl/annotation_station.cwl"
YAML="cwl/tests/annotation_station_config.yaml"

mkdir -p cwl/tests/test_results
RABIX_ARGS="--basedir cwl/tests/test_results"

rabix $RABIX_ARGS $CWL $YAML
