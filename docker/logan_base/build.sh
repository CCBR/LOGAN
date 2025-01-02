
# Build image
#docker buildx create --platform linux/amd64 --use
#docker buildx use upbeat_ganguly
#docker buildx inspect upbeat_ganguly
#docker buildx build --platform linux/amd64 -f Dockerfile -t dnousome/ccbr_logan_base:v0.3.0 -t dnousome/ccbr_logan_base:latest --push .

docker build --platform linux/amd64 --tag ccbr_logan_base:v0.3.8 -f Dockerfile . 

docker tag ccbr_logan_base:v0.3.8 dnousome/ccbr_logan_base:v0.3.8
docker tag ccbr_logan_base:v0.3.8 dnousome/ccbr_logan_base:latest


docker push dnousome/ccbr_logan_base:v0.3.8
docker push dnousome/ccbr_logan_base:latest

#Pull to CCBR
cd /data/CCBR_Pipeliner/SIFS
singularity pull dnousome-ccbr_logan_base-v0.3.8.img docker://dnousome/ccbr_logan_base:v0.3.8 
