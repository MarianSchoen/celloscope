#!/usr/bin/env bash

kubectl delete -f ./ingress.yaml -f ./service.yaml -f ./deployment.yaml -f ./pvc.yaml
