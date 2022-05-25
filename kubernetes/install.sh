#!/usr/bin/env bash

kubectl apply -f ./pvc.yaml
kubectl apply -f ./ingress.yaml -f ./service.yaml
kubectl apply -f ./deployment.yaml


