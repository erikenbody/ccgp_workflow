#!/bin/bash
export CLUSTER_NAME="snakemake"
export ZONE="us-west1-b"
export PROJECT_ID="ccgp-ucsc"

gcloud container clusters create $CLUSTER_NAME \
    --project=${PROJECT_ID} \
    --zone=${ZONE} \
    --machine-type="n1-standard-1" \
    --num-nodes=1 \

# gcloud container node-pools create "homestar-runner" \
#     --cluster=${CLUSTER_NAME} \
#     --machine-type="e2-custom-4-10240" \
#     --node-labels=machine_type=e2-custom-4-10240 \
#     --num-nodes=0 \
#     --zone=${ZONE} \
#     --preemptible \
#     --project=${PROJECT_ID} \
#     --scopes storage-rw \
#     --image-type=UBUNTU \
#     --disk-size=200GB \
#     --enable-autoscaling \
#     --min-nodes=0 \
#     --max-nodes=30 \

# gcloud container node-pools create "homestar-runner2" \
#     --cluster=${CLUSTER_NAME} \
#     --machine-type="e2-highcpu-16" \
#     --node-labels=machine_type=e2-highcpu-16 \
#     --num-nodes=0 \
#     --zone=${ZONE} \
#     --preemptible \
#     --project=${PROJECT_ID} \
#     --scopes storage-rw \
#     --image-type=UBUNTU \
#     --disk-size=200GB \
#     --enable-autoscaling \
#     --min-nodes=0 \
    # --max-nodes=3 \

gcloud container node-pools create "homestar-runner3" \
    --cluster=${CLUSTER_NAME} \
    --machine-type="n2d-standard-32" \
    --num-nodes=0 \
    --zone=${ZONE} \
    --preemptible \
    --project=${PROJECT_ID} \
    --scopes storage-rw \
    --image-type=UBUNTU \
    --disk-size=200GB \
    --enable-autoscaling \
    --min-nodes=0 \
    --max-nodes=5 \


gcloud container clusters get-credentials --zone=$ZONE $CLUSTER_NAME
