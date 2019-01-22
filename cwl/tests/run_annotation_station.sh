#!/bin/bash

CWL="cwl/mytool.cwl"
YAML="cwl/tests/mytool_config.yaml"

mkdir -p cwl/tests/test_results
RABIX_ARGS="--basedir cwl/tests/test_results"

rabix $RABIX_ARGS $CWL $YAML
